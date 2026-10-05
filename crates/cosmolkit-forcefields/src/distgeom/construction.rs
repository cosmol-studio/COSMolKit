//! Complete pinned distance-geometry force-field construction. The one private kernel
//! borrows caller-owned coordinates; no molecule or commit authority enters this owner.
use super::{ChiralSetPtr, ChiralViolationContribs, DistViolationContribs, FourthDimContribs};
use crate::{
    crystalff::{CrystalFFDetails, TorsionAngleContribs},
    geometry::Point3,
    kernel::{
        AngleConstraintContribs, DistanceConstraintContribs, ForceField, ForceFieldKernelError,
    },
    mmff::NonbondedContrib,
    uff::{InversionContribs, InversionContributionError},
};
use std::collections::BTreeMap;
const KNOWN_DIST_TOL: f64 = 0.01;
const KNOWN_DIST_FORCE_CONSTANT: f64 = 100.0;
/// Read-only bounds supplied by the existing conformer owner; this trait stores no second matrix.
pub trait DistanceBoundsRead {
    fn dimension(&self) -> usize;
    fn get_lower(&self, i: usize, j: usize) -> f64;
    fn get_upper(&self, i: usize, j: usize) -> f64;
}
trait Position3Read {
    fn point3(&self) -> Point3;
}
impl Position3Read for Vec<f64> {
    fn point3(&self) -> Point3 {
        point(self)
    }
}
impl Position3Read for &mut [f64] {
    fn point3(&self) -> Point3 {
        point(self)
    }
}
impl Position3Read for Point3 {
    fn point3(&self) -> Point3 {
        *self
    }
}
fn point(row: &[f64]) -> Point3 {
    Point3 {
        x: row[0],
        y: row[1],
        z: row[2],
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DistGeomForceFieldError {
    Kernel(ForceFieldKernelError),
    Inversion(InversionContributionError),
    IndexOverflow(usize),
}
impl From<ForceFieldKernelError> for DistGeomForceFieldError {
    fn from(e: ForceFieldKernelError) -> Self {
        Self::Kernel(e)
    }
}
impl From<InversionContributionError> for DistGeomForceFieldError {
    fn from(e: InversionContributionError) -> Self {
        Self::Inversion(e)
    }
}
fn source_index(i: usize) -> Result<u32, DistGeomForceFieldError> {
    u32::try_from(i).map_err(|_| DistGeomForceFieldError::IndexOverflow(i))
}
fn construct_distgeom_forcefield<'a>(
    mmat: &impl DistanceBoundsRead,
    positions: &'a mut [Vec<f64>],
    csets: &[ChiralSetPtr],
    weight_chiral: f64,
    weight_fourth_dim: f64,
    extra_weights: Option<&BTreeMap<(usize, usize), f64>>,
    basin_size_tol: f64,
    fixed_pts: Option<&[bool]>,
) -> Result<ForceField<'a>, DistGeomForceFieldError> {
    // BEGIN RDKIT CPP FUNCTION DistGeom::constructForceField (DistGeomUtils.cpp:184-253)
    // RDKit✔️✔️: ForceFields::ForceField *constructForceField(
    // RDKit✔️✔️:     const BoundsMatrix &mmat, RDGeom::PointPtrVect &positions,
    // RDKit✔️✔️:     const VECT_CHIRALSET &csets, double weightChiral, double weightFourthDim,
    // RDKit✔️✔️:     std::map<std::pair<int, int>, double> *extraWeights, double basinSizeTol,
    // RDKit✔️✔️:     boost::dynamic_bitset<> *fixedPts) {
    // RDKit✔️✔️:   unsigned int N = mmat.numRows();
    // RDKit✔️✔️:   CHECK_INVARIANT(N == positions.size(), "");
    // RDKit✔️✔️:   auto *field = new ForceFields::ForceField(positions[0]->dimension());
    // RDKit✔️✔️:   field->positions().insert(field->positions().begin(), positions.begin(),
    // RDKit✔️✔️:                             positions.end());
    // END RDKIT CPP FUNCTION DistGeom::constructForceField
    let n = mmat.dimension();
    assert_eq!(n, positions.len());
    assert!(!positions.is_empty());
    let dimension = positions[0].len();
    assert!((3..=4).contains(&dimension), "unsupported point dimension");
    assert!(
        positions.iter().all(|point| point.len() == dimension),
        "inconsistent point dimension"
    );
    if let Some(fixed_pts) = fixed_pts {
        assert!(fixed_pts.len() >= n, "bad fixed point bitset");
    }
    let mut field = ForceField::new(dimension as u32);
    field
        .positions_mut()
        .extend(positions.iter_mut().map(Vec::as_mut_slice));

    // RDKit✔️✔️:   auto contrib = new DistViolationContribs(field);
    // RDKit✔️✔️:   for (unsigned int i = 1; i < N; i++) {
    // RDKit✔️✔️:     for (unsigned int j = 0; j < i; j++) {
    // RDKit✔️✔️:       if (fixedPts != nullptr && (*fixedPts)[i] && (*fixedPts)[j]) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       double w = 1.0;
    // RDKit✔️✔️:       double l = mmat.getLowerBound(i, j);
    // RDKit✔️✔️:       double u = mmat.getUpperBound(i, j);
    // RDKit✔️✔️:       bool includeIt = false;
    // RDKit✔️✔️:       if (extraWeights) {
    // RDKit✔️✔️:         std::map<std::pair<int, int>, double>::const_iterator mapIt;
    // RDKit✔️✔️:         mapIt = extraWeights->find(std::make_pair(i, j));
    // RDKit✔️✔️:         if (mapIt != extraWeights->end()) {
    // RDKit✔️✔️:           w = mapIt->second;
    // RDKit✔️✔️:           includeIt = true;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (u - l <= basinSizeTol) {
    // RDKit✔️✔️:         includeIt = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (includeIt) {
    // RDKit✔️✔️:         contrib->addContrib(i, j, u, l, w);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!contrib->empty()) {
    // RDKit✔️✔️:     field->contribs().push_back(ForceFields::ContribPtr(contrib));
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     delete contrib;
    // RDKit✔️✔️:   }
    let mut dist_contrib = DistViolationContribs::new(&field);
    for i in 1..n {
        for j in 0..i {
            if fixed_pts.is_some_and(|fixed_pts| fixed_pts[i] && fixed_pts[j]) {
                continue;
            }
            let mut weight = 1.0;
            let lower = mmat.get_lower(i, j);
            let upper = mmat.get_upper(i, j);
            let mut include_it = false;
            if let Some(extra_weights) = extra_weights
                && let Some(extra_weight) = extra_weights.get(&(i, j))
            {
                weight = *extra_weight;
                include_it = true;
            }
            if upper - lower <= basin_size_tol {
                include_it = true;
            }
            if include_it {
                dist_contrib.add_contrib(i, j, upper, lower, weight);
            }
        }
    }
    if !dist_contrib.empty() {
        field.add_contribution(Box::new(dist_contrib));
    }

    // RDKit✔️✔️:   // now add chiral constraints
    // RDKit✔️✔️:   if (weightChiral > 1.0e-8) {
    // RDKit✔️✔️:     auto contrib = new ChiralViolationContribs(field);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     for (const auto &cset : csets) {
    // RDKit✔️✔️:       contrib->addContrib(cset.get(), weightChiral);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!contrib->empty()) {
    // RDKit✔️✔️:       field->contribs().push_back(ForceFields::ContribPtr(contrib));
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       delete contrib;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if weight_chiral > 1.0e-8 {
        let mut chiral_contrib = ChiralViolationContribs::new(&field);
        for cset in csets {
            chiral_contrib.add_contrib(cset, weight_chiral);
        }
        if !chiral_contrib.empty() {
            field.add_contribution(Box::new(chiral_contrib));
        }
    }

    // RDKit✔️✔️:   // finally the contribution from the fourth dimension if we need to
    // RDKit✔️✔️:   if ((field->dimension() == 4) && (weightFourthDim > 1.0e-8)) {
    // RDKit✔️✔️:     auto contrib = new FourthDimContribs(field);
    // RDKit✔️✔️:     for (unsigned int i = 0; i < N; i++) {
    // RDKit✔️✔️:       contrib->addContrib(i, weightFourthDim);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!contrib->empty()) {
    // RDKit✔️✔️:       field->contribs().push_back(ForceFields::ContribPtr(contrib));
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       delete contrib;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return field;
    // RDKit✔️✔️: }  // constructForceField
    if field.dimension() == 4 && weight_fourth_dim > 1.0e-8 {
        let mut fourth_dim_contrib = FourthDimContribs::new(&field);
        for i in 0..n {
            fourth_dim_contrib.add_contrib(i, weight_fourth_dim);
        }
        if !fourth_dim_contrib.empty() {
            field.add_contribution(Box::new(fourth_dim_contrib));
        }
    }
    Ok(field)
}

fn add_improper_torsion_terms(
    ff: &mut ForceField<'_>,
    force_scaling_factor: f64,
    improper_atoms: &[Vec<i32>],
    is_improper_constrained: &mut [bool],
) -> Result<(), DistGeomForceFieldError> {
    // BEGIN RDKIT CPP FUNCTION DistGeom::addImproperTorsionTerms (DistGeomUtils.cpp:255-307)
    // RDKit✔️✔️: void addImproperTorsionTerms(ForceFields::ForceField *ff,
    // RDKit✔️✔️:                              double forceScalingFactor,
    // RDKit✔️✔️:                              const std::vector<std::vector<int>> &improperAtoms,
    // RDKit✔️✔️:                              boost::dynamic_bitset<> &isImproperConstrained) {
    // RDKit✔️✔️:   PRECONDITION(ff, "bad force field");
    // RDKit✔️✔️:   auto inversionContribs =
    // RDKit✔️✔️:       std::make_unique<ForceFields::UFF::InversionContribs>(ff);
    // END RDKIT CPP FUNCTION DistGeom::addImproperTorsionTerms
    let mut inversion_contribs = InversionContribs::default();

    // RDKit✔️✔️:   for (const auto &improperAtom : improperAtoms) {
    // RDKit✔️✔️:     std::vector<int> n(4);
    for improper_atom in improper_atoms {
        let mut n = [0_usize; 4];
        // RDKit✔️✔️:     for (unsigned int i = 0; i < 3; ++i) {
        for i in 0..3 {
            // RDKit✔️✔️:       n[1] = 1;
            n[1] = 1;
            // RDKit✔️✔️:       switch (i) {
            match i {
                // RDKit✔️✔️:         case 0:
                // RDKit✔️✔️:           n[0] = 0;
                // RDKit✔️✔️:           n[2] = 2;
                // RDKit✔️✔️:           n[3] = 3;
                // RDKit✔️✔️:           break;
                0 => {
                    n[0] = 0;
                    n[2] = 2;
                    n[3] = 3;
                }
                // RDKit✔️✔️:         case 1:
                // RDKit✔️✔️:           n[0] = 0;
                // RDKit✔️✔️:           n[2] = 3;
                // RDKit✔️✔️:           n[3] = 2;
                // RDKit✔️✔️:           break;
                1 => {
                    n[0] = 0;
                    n[2] = 3;
                    n[3] = 2;
                }
                // RDKit✔️✔️:         case 2:
                // RDKit✔️✔️:           n[0] = 2;
                // RDKit✔️✔️:           n[2] = 3;
                // RDKit✔️✔️:           n[3] = 0;
                // RDKit✔️✔️:           break;
                2 => {
                    n[0] = 2;
                    n[2] = 3;
                    n[3] = 0;
                }
                _ => unreachable!("loop bounds guarantee 0..3"),
            }

            // RDKit✔️✔️:       inversionContribs->addContrib(
            // RDKit✔️✔️:           improperAtom[n[0]], improperAtom[n[1]], improperAtom[n[2]],
            // RDKit✔️✔️:           improperAtom[n[3]], improperAtom[4],
            // RDKit✔️✔️:           static_cast<bool>(improperAtom[5]), forceScalingFactor);
            inversion_contribs.add_contrib(
                ff.positions(),
                improper_atom[n[0]] as u32,
                improper_atom[n[1]] as u32,
                improper_atom[n[2]] as u32,
                improper_atom[n[3]] as u32,
                improper_atom[4],
                improper_atom[5] != 0,
                force_scaling_factor,
            )?;

            // RDKit✔️✔️:       isImproperConstrained[improperAtom[n[1]]] = 1;
            is_improper_constrained[improper_atom[n[1]] as usize] = true;
        }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!inversionContribs->empty()) {
    // RDKit✔️✔️:     ff->contribs().push_back(std::move(inversionContribs));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    if !inversion_contribs.empty() {
        ff.add_contribution(Box::new(inversion_contribs));
    }
    Ok(())
}

fn add_experimental_torsion_terms(
    ff: &mut ForceField<'_>,
    etkdg_details: &CrystalFFDetails,
    atom_pairs: &mut [bool],
    num_atoms: usize,
) -> Result<(), DistGeomForceFieldError> {
    // BEGIN RDKIT CPP FUNCTION DistGeom::addExperimentalTorsionTerms (DistGeomUtils.cpp:309-340)
    // RDKit✔️✔️: void addExperimentalTorsionTerms(
    // RDKit✔️✔️:     ForceFields::ForceField *ff,
    // RDKit✔️✔️:     const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails,
    // RDKit✔️✔️:     boost::dynamic_bitset<> &atomPairs, unsigned int numAtoms) {
    // RDKit✔️✔️:   PRECONDITION(ff, "bad force field");
    // RDKit✔️✔️:   auto torsionContribs =
    // RDKit✔️✔️:       std::make_unique<ForceFields::CrystalFF::TorsionAngleContribs>(ff);
    // END RDKIT CPP FUNCTION DistGeom::addExperimentalTorsionTerms
    let mut torsion_contribs = TorsionAngleContribs::new(ff);

    // RDKit✔️✔️:   for (unsigned int t = 0; t < etkdgDetails.expTorsionAtoms.size(); ++t) {
    for t in 0..etkdg_details.exp_torsion_atoms.len() {
        // RDKit✔️✔️:     int i = etkdgDetails.expTorsionAtoms[t][0];
        // RDKit✔️✔️:     int j = etkdgDetails.expTorsionAtoms[t][1];
        // RDKit✔️✔️:     int k = etkdgDetails.expTorsionAtoms[t][2];
        // RDKit✔️✔️:     int l = etkdgDetails.expTorsionAtoms[t][3];
        let i = etkdg_details.exp_torsion_atoms[t][0] as usize;
        let j = etkdg_details.exp_torsion_atoms[t][1] as usize;
        let k = etkdg_details.exp_torsion_atoms[t][2] as usize;
        let l = etkdg_details.exp_torsion_atoms[t][3] as usize;

        // RDKit✔️✔️:     if (i < l) {
        // RDKit✔️✔️:       atomPairs[i * numAtoms + l] = 1;
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       atomPairs[l * numAtoms + i] = 1;
        // RDKit✔️✔️:     }
        if i < l {
            atom_pairs[i * num_atoms + l] = true;
        } else {
            atom_pairs[l * num_atoms + i] = true;
        }

        // RDKit✔️✔️:     torsionContribs->addContrib(i, j, k, l,
        // RDKit✔️✔️:                                 etkdgDetails.expTorsionAngles[t].second,
        // RDKit✔️✔️:                                 etkdgDetails.expTorsionAngles[t].first);
        torsion_contribs.add_contrib(
            i,
            j,
            k,
            l,
            etkdg_details.exp_torsion_angles[t].1.clone(),
            etkdg_details.exp_torsion_angles[t].0.clone(),
        );
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!torsionContribs->empty()) {
    // RDKit✔️✔️:     ff->contribs().push_back(std::move(torsionContribs));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    if !torsion_contribs.is_empty() {
        ff.add_contribution(Box::new(torsion_contribs));
    }
    Ok(())
}

fn build_12_terms(
    ff: &ForceField<'_>,
    etkdg_details: &CrystalFFDetails,
    atom_pairs: &mut [bool],
    positions: &[impl Position3Read],
    force_constant: f64,
    num_atoms: usize,
) -> Result<DistanceConstraintContribs, DistGeomForceFieldError> {
    // BEGIN RDKIT CPP FUNCTION DistGeom::add12Terms (DistGeomUtils.cpp:342-383)
    // RDKit✔️✔️: void add12Terms(ForceFields::ForceField *ff,
    // RDKit✔️✔️:                 const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails,
    // RDKit✔️✔️:                 boost::dynamic_bitset<> &atomPairs,
    // RDKit✔️✔️:                 RDGeom::Point3DPtrVect &positions, double forceConstant,
    // RDKit✔️✔️:                 unsigned int numAtoms) {
    // RDKit✔️✔️:   PRECONDITION(ff, "bad force field");
    // RDKit✔️✔️:   auto distContribs =
    // RDKit✔️✔️:       std::make_unique<ForceFields::DistanceConstraintContribs>(ff);
    // END RDKIT CPP FUNCTION DistGeom::add12Terms
    let mut dist_contribs = DistanceConstraintContribs::new(ff);

    // RDKit✔️✔️:   for (const auto &bond : etkdgDetails.bonds) {
    for &(first, second) in &etkdg_details.bonds {
        // RDKit✔️✔️:     unsigned int i = bond.first;
        // RDKit✔️✔️:     unsigned int j = bond.second;
        let i = first as usize;
        let j = second as usize;

        // RDKit✔️✔️:     if (i < j) {
        // RDKit✔️✔️:       atomPairs[i * numAtoms + j] = 1;
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       atomPairs[j * numAtoms + i] = 1;
        // RDKit✔️✔️:     }
        if i < j {
            atom_pairs[i * num_atoms + j] = true;
        } else {
            atom_pairs[j * num_atoms + i] = true;
        }

        // RDKit✔️✔️:     double d = ((*positions[i]) - (*positions[j])).length();
        let d = Point3::difference(&positions[i].point3(), &positions[j].point3()).length();
        // RDKit✔️✔️:     distContribs->addContrib(i, j, d - KNOWN_DIST_TOL, d + KNOWN_DIST_TOL,
        // RDKit✔️✔️:                              forceConstant);
        dist_contribs.add_contrib(
            ff,
            source_index(i)?,
            source_index(j)?,
            d - KNOWN_DIST_TOL,
            d + KNOWN_DIST_TOL,
            force_constant,
        )?;
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!distContribs->empty()) {
    // RDKit✔️✔️:     ff->contribs().push_back(std::move(distContribs));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    Ok(dist_contribs)
}

fn build_13_terms(
    ff: &ForceField<'_>,
    etkdg_details: &CrystalFFDetails,
    atom_pairs: &mut [bool],
    positions: &[impl Position3Read],
    force_constant: f64,
    is_improper_constrained: &[bool],
    use_basic_knowledge: bool,
    mmat: &impl DistanceBoundsRead,
    num_atoms: usize,
) -> Result<(AngleConstraintContribs, DistanceConstraintContribs), DistGeomForceFieldError> {
    // BEGIN RDKIT CPP FUNCTION DistGeom::add13Terms (DistGeomUtils.cpp:385-455)
    // RDKit✔️✔️: void add13Terms(ForceFields::ForceField *ff,
    // RDKit✔️✔️:                 const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails,
    // RDKit✔️✔️:                 boost::dynamic_bitset<> &atomPairs,
    // RDKit✔️✔️:                 RDGeom::Point3DPtrVect &positions, double forceConstant,
    // RDKit✔️✔️:                 const boost::dynamic_bitset<> &isImproperConstrained,
    // RDKit✔️✔️:                 bool useBasicKnowledge, const BoundsMatrix &mmat,
    // RDKit✔️✔️:                 unsigned int numAtoms) {
    // RDKit✔️✔️:   PRECONDITION(ff, "bad force field");
    // RDKit✔️✔️:   auto distContribs =
    // RDKit✔️✔️:       std::make_unique<ForceFields::DistanceConstraintContribs>(ff);
    // RDKit✔️✔️:   auto angleContribs =
    // RDKit✔️✔️:       std::make_unique<ForceFields::AngleConstraintContribs>(ff);
    // END RDKIT CPP FUNCTION DistGeom::add13Terms
    let mut dist_contribs = DistanceConstraintContribs::new(ff);
    let mut angle_contribs = AngleConstraintContribs::new(ff);

    // RDKit✔️✔️:   for (const auto &angle : etkdgDetails.angles) {
    for angle in &etkdg_details.angles {
        // RDKit✔️✔️:     unsigned int i = angle[0];
        // RDKit✔️✔️:     unsigned int j = angle[1];
        // RDKit✔️✔️:     unsigned int k = angle[2];
        let i = angle[0] as usize;
        let j = angle[1] as usize;
        let k = angle[2] as usize;

        // RDKit✔️✔️:     if (i < k) {
        // RDKit✔️✔️:       atomPairs[i * numAtoms + k] = 1;
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       atomPairs[k * numAtoms + i] = 1;
        // RDKit✔️✔️:     }
        if i < k {
            atom_pairs[i * num_atoms + k] = true;
        } else {
            atom_pairs[k * num_atoms + i] = true;
        }

        // RDKit✔️✔️:     // check for triple bonds
        // RDKit✔️✔️:     if (useBasicKnowledge && angle[3]) {
        // RDKit✔️✔️:       angleContribs->addContrib(i, j, k, 179.0, 180.0, 1);
        if use_basic_knowledge && angle[3] != 0 {
            angle_contribs.add_contrib(
                ff,
                source_index(i)?,
                source_index(j)?,
                source_index(k)?,
                179.0,
                180.0,
                1.0,
            )?;
        // RDKit✔️✔️:     } else if (isImproperConstrained[j]) {
        // RDKit✔️✔️:       distContribs->addContrib(i, k, mmat.getLowerBound(i, k),
        // RDKit✔️✔️:                                mmat.getUpperBound(i, k), forceConstant);
        } else if is_improper_constrained[j] {
            dist_contribs.add_contrib(
                ff,
                source_index(i)?,
                source_index(k)?,
                mmat.get_lower(i, k),
                mmat.get_upper(i, k),
                force_constant,
            )?;
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       double d = ((*positions[i]) - (*positions[k])).length();
        // RDKit✔️✔️:       distContribs->addContrib(i, k, d - KNOWN_DIST_TOL, d + KNOWN_DIST_TOL,
        // RDKit✔️✔️:                                forceConstant);
        // RDKit✔️✔️:     }
        } else {
            let d = Point3::difference(&positions[i].point3(), &positions[k].point3()).length();
            dist_contribs.add_contrib(
                ff,
                source_index(i)?,
                source_index(k)?,
                d - KNOWN_DIST_TOL,
                d + KNOWN_DIST_TOL,
                force_constant,
            )?;
        }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!angleContribs->empty()) {
    // RDKit✔️✔️:     ff->contribs().push_back(std::move(angleContribs));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!distContribs->empty()) {
    // RDKit✔️✔️:     ff->contribs().push_back(std::move(distContribs));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    Ok((angle_contribs, dist_contribs))
}

fn build_long_range_distance_constraints(
    ff: &ForceField<'_>,
    etkdg_details: &CrystalFFDetails,
    atom_pairs: &[bool],
    positions: &[impl Position3Read],
    known_distance_force_constant: f64,
    mmat: &impl DistanceBoundsRead,
    num_atoms: usize,
) -> Result<DistanceConstraintContribs, DistGeomForceFieldError> {
    // BEGIN RDKIT CPP FUNCTION DistGeom::addLongRangeDistanceConstraints (DistGeomUtils.cpp:457-505)
    // RDKit✔️✔️: void addLongRangeDistanceConstraints(
    // RDKit✔️✔️:     ForceFields::ForceField *ff,
    // RDKit✔️✔️:     const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails,
    // RDKit✔️✔️:     const boost::dynamic_bitset<> &atomPairs, RDGeom::Point3DPtrVect &positions,
    // RDKit✔️✔️:     double knownDistanceForceConstant, const BoundsMatrix &mmat,
    // RDKit✔️✔️:     unsigned int numAtoms) {
    // RDKit✔️✔️:   PRECONDITION(ff, "bad force field");
    // RDKit✔️✔️:   auto distContribs =
    // RDKit✔️✔️:       std::make_unique<ForceFields::DistanceConstraintContribs>(ff);
    // RDKit✔️✔️:   double fdist = knownDistanceForceConstant;
    // END RDKIT CPP FUNCTION DistGeom::addLongRangeDistanceConstraints
    let mut dist_contribs = DistanceConstraintContribs::new(ff);

    // RDKit✔️✔️:   for (unsigned int i = 1; i < numAtoms; ++i) {
    for i in 1..num_atoms {
        // RDKit✔️✔️:     for (unsigned int j = 0; j < i; ++j) {
        for j in 0..i {
            // RDKit✔️✔️:       if (!atomPairs[j * numAtoms + i]) {
            if !atom_pairs[j * num_atoms + i] {
                // RDKit✔️✔️:         fdist = etkdgDetails.boundsMatForceScaling * 10.0;
                let mut fdist = etkdg_details.bounds_mat_force_scaling * 10.0;
                // RDKit✔️✔️:         double l = mmat.getLowerBound(i, j);
                // RDKit✔️✔️:         double u = mmat.getUpperBound(i, j);
                let mut l = mmat.get_lower(i, j);
                let mut u = mmat.get_upper(i, j);

                // RDKit✔️✔️:         if (!etkdgDetails.constrainedAtoms.empty() &&
                // RDKit✔️✔️:             etkdgDetails.constrainedAtoms[i] &&
                // RDKit✔️✔️:             etkdgDetails.constrainedAtoms[j]) {
                if !etkdg_details.constrained_atoms.is_empty()
                    && etkdg_details.constrained_atoms[i]
                    && etkdg_details.constrained_atoms[j]
                {
                    // RDKit✔️✔️:           // we're constrained, so use very tight bounds
                    // RDKit✔️✔️:           l = u = ((*positions[i]) - (*positions[j])).length();
                    let d =
                        Point3::difference(&positions[i].point3(), &positions[j].point3()).length();
                    l = d;
                    u = d;
                    // RDKit✔️✔️:           l -= KNOWN_DIST_TOL;
                    // RDKit✔️✔️:           u += KNOWN_DIST_TOL;
                    l -= KNOWN_DIST_TOL;
                    u += KNOWN_DIST_TOL;
                    // RDKit✔️✔️:           fdist = knownDistanceForceConstant;
                    fdist = known_distance_force_constant;
                    // RDKit✔️✔️:         }
                }
                // RDKit✔️✔️:         distContribs->addContrib(i, j, l, u, fdist);
                dist_contribs.add_contrib(ff, source_index(i)?, source_index(j)?, l, u, fdist)?;
            }
            // RDKit✔️✔️:       }
        }
        // RDKit✔️✔️:     }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!distContribs->empty()) {
    // RDKit✔️✔️:     ff->contribs().push_back(std::move(distContribs));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    Ok(dist_contribs)
}

fn construct_3d_forcefield<'a>(
    mmat: &impl DistanceBoundsRead,
    positions: &'a mut [Vec<f64>],
    etkdg_details: &CrystalFFDetails,
) -> Result<ForceField<'a>, DistGeomForceFieldError> {
    // BEGIN RDKIT CPP FUNCTION DistGeom::construct3DForceField (DistGeomUtils.cpp:495-524)
    // RDKit✔️✔️: ForceFields::ForceField *construct3DForceField(
    // RDKit✔️✔️:     const BoundsMatrix &mmat, RDGeom::Point3DPtrVect &positions,
    // RDKit✔️✔️:     const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails) {
    // RDKit✔️✔️:   unsigned int N = mmat.numRows();
    // RDKit✔️✔️:   CHECK_INVARIANT(N == positions.size(), "");
    // RDKit✔️✔️:   CHECK_INVARIANT(etkdgDetails.expTorsionAtoms.size() ==
    // RDKit✔️✔️:                       etkdgDetails.expTorsionAngles.size(),
    // RDKit✔️✔️:                   "");
    let n = mmat.dimension();
    assert_eq!(n, positions.len());
    assert_eq!(
        etkdg_details.exp_torsion_atoms.len(),
        etkdg_details.exp_torsion_angles.len()
    );

    // RDKit✔️✔️:   auto *field = new ForceFields::ForceField(positions[0]->dimension());
    // RDKit✔️✔️:   field->positions().insert(field->positions().begin(), positions.begin(),
    // RDKit✔️✔️:                             positions.end());
    // Rust stores RDGeom::Point3D-equivalent coordinates by value; Point3D
    // dimension is therefore fixed at three for this overload.
    assert!(!positions.is_empty());
    let mut field = ForceField::new(3);
    assert!(positions.iter().all(|point| point.len() == 3));
    field
        .positions_mut()
        .extend(positions.iter_mut().map(Vec::as_mut_slice));

    // RDKit✔️✔️:   // keep track which atoms are 1,2-, 1,3- or 1,4-restrained
    // RDKit✔️✔️:   boost::dynamic_bitset<> atomPairs(N * N);
    // RDKit✔️✔️:   // don't add 1-3 Distances constraints for angles where the
    // RDKit✔️✔️:   // central atom of the angle is the central atom of an improper torsion.
    // RDKit✔️✔️:   boost::dynamic_bitset<> isImproperConstrained(N);
    let mut atom_pairs = vec![false; n * n];
    let mut is_improper_constrained = vec![false; n];

    // RDKit✔️✔️:   addExperimentalTorsionTerms(field, etkdgDetails, atomPairs, N);
    // RDKit✔️✔️:   addImproperTorsionTerms(field, 10.0, etkdgDetails.improperAtoms,
    // RDKit✔️✔️:                           isImproperConstrained);
    // RDKit✔️✔️:   add12Terms(field, etkdgDetails, atomPairs, positions,
    // RDKit✔️✔️:              KNOWN_DIST_FORCE_CONSTANT, N);
    // RDKit✔️✔️:   add13Terms(field, etkdgDetails, atomPairs, positions,
    // RDKit✔️✔️:              KNOWN_DIST_FORCE_CONSTANT, isImproperConstrained, true, mmat, N);
    // RDKit✔️✔️:   // minimum distance for all other atom pairs that aren't constrained
    // RDKit✔️✔️:   addLongRangeDistanceConstraints(field, etkdgDetails, atomPairs, positions,
    // RDKit✔️✔️:                                   KNOWN_DIST_FORCE_CONSTANT, mmat, N);
    // RDKit✔️✔️:   return field;
    // RDKit✔️✔️: }  // construct3DForceField
    add_experimental_torsion_terms(&mut field, etkdg_details, &mut atom_pairs, n)?;
    add_improper_torsion_terms(
        &mut field,
        10.0,
        &etkdg_details.improper_atoms,
        &mut is_improper_constrained,
    )?;
    let group = build_12_terms(
        &field,
        etkdg_details,
        &mut atom_pairs,
        field.positions(),
        KNOWN_DIST_FORCE_CONSTANT,
        n,
    )?;
    if !group.empty() {
        field.add_contribution(Box::new(group));
    }
    let groups = build_13_terms(
        &field,
        etkdg_details,
        &mut atom_pairs,
        field.positions(),
        KNOWN_DIST_FORCE_CONSTANT,
        &is_improper_constrained,
        true,
        mmat,
        n,
    )?;
    if !groups.0.empty() {
        field.add_contribution(Box::new(groups.0));
    }
    if !groups.1.empty() {
        field.add_contribution(Box::new(groups.1));
    }
    let group = build_long_range_distance_constraints(
        &field,
        etkdg_details,
        &atom_pairs,
        field.positions(),
        KNOWN_DIST_FORCE_CONSTANT,
        mmat,
        n,
    )?;
    if !group.empty() {
        field.add_contribution(Box::new(group));
    }
    Ok(field)
}

fn construct_plain_3d_forcefield<'a>(
    mmat: &impl DistanceBoundsRead,
    positions: &'a mut [Vec<f64>],
    etkdg_details: &CrystalFFDetails,
) -> Result<ForceField<'a>, DistGeomForceFieldError> {
    // BEGIN RDKIT CPP FUNCTION DistGeom::constructPlain3DForceField (DistGeomUtils.cpp:545-573)
    // RDKit✔️✔️: ForceFields::ForceField *constructPlain3DForceField(
    // RDKit✔️✔️:     const BoundsMatrix &mmat, RDGeom::Point3DPtrVect &positions,
    // RDKit✔️✔️:     const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails) {
    // RDKit✔️✔️:   unsigned int N = mmat.numRows();
    // RDKit✔️✔️:   CHECK_INVARIANT(N == positions.size(), "");
    // RDKit✔️✔️:   CHECK_INVARIANT(etkdgDetails.expTorsionAtoms.size() ==
    // RDKit✔️✔️:                       etkdgDetails.expTorsionAngles.size(),
    // RDKit✔️✔️:                   "");
    let n = mmat.dimension();
    assert_eq!(n, positions.len());
    assert_eq!(
        etkdg_details.exp_torsion_atoms.len(),
        etkdg_details.exp_torsion_angles.len()
    );

    // RDKit✔️✔️:   auto *field = new ForceFields::ForceField(positions[0]->dimension());
    // RDKit✔️✔️:   field->positions().insert(field->positions().begin(), positions.begin(),
    // RDKit✔️✔️:                             positions.end());
    assert!(!positions.is_empty());
    let mut field = ForceField::new(3);
    assert!(positions.iter().all(|point| point.len() == 3));
    field
        .positions_mut()
        .extend(positions.iter_mut().map(Vec::as_mut_slice));

    // RDKit✔️✔️:   // keep track which atoms are 1,2-, 1,3- or 1,4-restrained
    // RDKit✔️✔️:   boost::dynamic_bitset<> atomPairs(N * N);
    // RDKit✔️✔️:   // don't add 1-3 Distances constraints for angles where the
    // RDKit✔️✔️:   // central atom of the angle is the central atom of an improper torsion.
    // RDKit✔️✔️:   boost::dynamic_bitset<> isImproperConstrained(N);
    let mut atom_pairs = vec![false; n * n];
    let is_improper_constrained = vec![false; n];

    // RDKit✔️✔️:   addExperimentalTorsionTerms(field, etkdgDetails, atomPairs, N);
    // RDKit✔️✔️:   add12Terms(field, etkdgDetails, atomPairs, positions,
    // RDKit✔️✔️:              KNOWN_DIST_FORCE_CONSTANT, N);
    // RDKit✔️✔️:   add13Terms(field, etkdgDetails, atomPairs, positions,
    // RDKit✔️✔️:              KNOWN_DIST_FORCE_CONSTANT, isImproperConstrained, false, mmat, N);
    // RDKit✔️✔️:   // minimum distance for all other atom pairs that aren't constrained
    // RDKit✔️✔️:   addLongRangeDistanceConstraints(field, etkdgDetails, atomPairs, positions,
    // RDKit✔️✔️:                                   KNOWN_DIST_FORCE_CONSTANT, mmat, N);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return field;
    // RDKit✔️✔️: }  // constructPlain3DForceField
    add_experimental_torsion_terms(&mut field, etkdg_details, &mut atom_pairs, n)?;
    let group = build_12_terms(
        &field,
        etkdg_details,
        &mut atom_pairs,
        field.positions(),
        KNOWN_DIST_FORCE_CONSTANT,
        n,
    )?;
    if !group.empty() {
        field.add_contribution(Box::new(group));
    }
    let groups = build_13_terms(
        &field,
        etkdg_details,
        &mut atom_pairs,
        field.positions(),
        KNOWN_DIST_FORCE_CONSTANT,
        &is_improper_constrained,
        false,
        mmat,
        n,
    )?;
    if !groups.0.empty() {
        field.add_contribution(Box::new(groups.0));
    }
    if !groups.1.empty() {
        field.add_contribution(Box::new(groups.1));
    }
    let group = build_long_range_distance_constraints(
        &field,
        etkdg_details,
        &atom_pairs,
        field.positions(),
        KNOWN_DIST_FORCE_CONSTANT,
        mmat,
        n,
    )?;
    if !group.empty() {
        field.add_contribution(Box::new(group));
    }
    Ok(field)
}

fn construct_3d_improper_forcefield_from_parts<'a>(
    mmat: &impl DistanceBoundsRead,
    positions: &'a mut [Vec<f64>],
    improper_atoms: &[Vec<i32>],
    angles: &[Vec<i32>],
    atom_nums: &[i32],
) -> Result<ForceField<'a>, DistGeomForceFieldError> {
    // BEGIN RDKIT CPP FUNCTION DistGeom::construct3DImproperForceField (DistGeomUtils.cpp:575-612)
    // RDKit✔️✔️: ForceFields::ForceField *construct3DImproperForceField(
    // RDKit✔️✔️:     const BoundsMatrix &mmat, RDGeom::Point3DPtrVect &positions,
    // RDKit✔️✔️:     const std::vector<std::vector<int>> &improperAtoms,
    // RDKit✔️✔️:     const std::vector<std::vector<int>> &angles,
    // RDKit✔️✔️:     const std::vector<int> &atomNums) {
    // RDKit✔️✔️:   RDUNUSED_PARAM(atomNums);
    let _ = atom_nums;
    // RDKit✔️✔️:   unsigned int N = mmat.numRows();
    // RDKit✔️✔️:   CHECK_INVARIANT(N == positions.size(), "");
    let n = mmat.dimension();
    assert_eq!(n, positions.len());

    // RDKit✔️✔️:   auto *field = new ForceFields::ForceField(positions[0]->dimension());
    // RDKit✔️✔️:   field->positions().insert(field->positions().begin(), positions.begin(),
    // RDKit✔️✔️:                             positions.end());
    assert!(!positions.is_empty());
    let mut field = ForceField::new(3);
    assert!(positions.iter().all(|point| point.len() == 3));
    field
        .positions_mut()
        .extend(positions.iter_mut().map(Vec::as_mut_slice));

    // RDKit✔️✔️:   // improper torsions / out-of-plane bend / inversion
    // RDKit✔️✔️:   double oobForceScalingFactor = 10.0;
    // RDKit✔️✔️:   boost::dynamic_bitset<> isImproperConstrained(N);
    // RDKit✔️✔️:   addImproperTorsionTerms(field, oobForceScalingFactor, improperAtoms,
    // RDKit✔️✔️:                           isImproperConstrained);
    let oob_force_scaling_factor = 10.0;
    let mut is_improper_constrained = vec![false; n];
    add_improper_torsion_terms(
        &mut field,
        oob_force_scaling_factor,
        improper_atoms,
        &mut is_improper_constrained,
    )?;

    // RDKit✔️✔️:   // Check that SP Centers have an angle of 180 degrees.
    // RDKit✔️✔️:   auto angleContribs =
    // RDKit✔️✔️:       std::make_unique<ForceFields::AngleConstraintContribs>(field);
    let mut angle_contribs = AngleConstraintContribs::new(&field);
    // RDKit✔️✔️:   for (const auto &angle : angles) {
    for angle in angles {
        // RDKit✔️✔️:     if (angle[3]) {
        if angle[3] != 0 {
            // RDKit✔️✔️:       angleContribs->addContrib(angle[0], angle[1], angle[2], 179.0, 180.0,
            // RDKit✔️✔️:                                 oobForceScalingFactor);
            angle_contribs.add_contrib(
                &field,
                angle[0] as u32,
                angle[1] as u32,
                angle[2] as u32,
                179.0,
                180.0,
                oob_force_scaling_factor,
            )?;
        }
        // RDKit✔️✔️:     }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!angleContribs->empty()) {
    // RDKit✔️✔️:     field->contribs().push_back(std::move(angleContribs));
    // RDKit✔️✔️:   }
    if !angle_contribs.empty() {
        field.add_contribution(Box::new(angle_contribs));
    }
    // RDKit✔️✔️:   return field;
    // RDKit✔️✔️: }  // construct3DImproperForceField
    Ok(field)
}

fn construct_3d_improper_forcefield<'a>(
    mmat: &impl DistanceBoundsRead,
    positions: &'a mut [Vec<f64>],
    etkdg_details: &CrystalFFDetails,
) -> Result<ForceField<'a>, DistGeomForceFieldError> {
    // BEGIN RDKIT CPP INLINE OVERLOAD DistGeom::construct3DImproperForceField (DistGeomUtils.h:215-220)
    // RDKit✔️✔️: inline ForceFields::ForceField *construct3DImproperForceField(
    // RDKit✔️✔️:     const BoundsMatrix &mmat, RDGeom::Point3DPtrVect &positions,
    // RDKit✔️✔️:     const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails) {
    // RDKit✔️✔️:   return construct3DImproperForceField(
    // RDKit✔️✔️:       mmat, positions, etkdgDetails.improperAtoms, etkdgDetails.angles,
    // RDKit✔️✔️:       etkdgDetails.atomNums);
    // RDKit✔️✔️: }
    construct_3d_improper_forcefield_from_parts(
        mmat,
        positions,
        &etkdg_details.improper_atoms,
        &etkdg_details.angles,
        &etkdg_details.atom_nums,
    )
}

fn construct_3d_forcefield_with_cpci<'a>(
    mmat: &impl DistanceBoundsRead,
    positions: &'a mut [Vec<f64>],
    etkdg_details: &CrystalFFDetails,
    cpci: &BTreeMap<(usize, usize), f64>,
) -> Result<ForceField<'a>, DistGeomForceFieldError> {
    // BEGIN RDKIT CPP FUNCTION DistGeom::construct3DForceField CPCI overload (DistGeomUtils.cpp:526-543)
    // RDKit✔️✔️: ForceFields::ForceField *construct3DForceField(
    // RDKit✔️✔️:     const BoundsMatrix &mmat, RDGeom::Point3DPtrVect &positions,
    // RDKit✔️✔️:     const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails,
    // RDKit✔️✔️:     const std::map<std::pair<unsigned int, unsigned int>, double> &CPCI) {
    // RDKit✔️✔️:   auto *field = construct3DForceField(mmat, positions, etkdgDetails);
    let mut field = construct_3d_forcefield(mmat, positions, etkdg_details)?;

    // RDKit✔️✔️:   bool is1_4 = false;
    // RDKit✔️✔️:   // double dielConst = 1.0;
    // RDKit✔️✔️:   boost::uint8_t dielModel = 1;
    // RDKit✔️✔️:   auto *contrib = new ForceFields::MMFF::EleContrib(field);
    let is_1_4 = false;
    let diel_model = 1;
    let mut contrib = NonbondedContrib::new(&field);

    // RDKit✔️✔️:   field->contribs().emplace_back(contrib);
    // RDKit✔️✔️:   for (const auto &charge : CPCI) {
    // RDKit✔️✔️:     contrib->addTerm(charge.first.first, charge.first.second, charge.second,
    // RDKit✔️✔️:                      dielModel, is1_4);
    // RDKit✔️✔️:   }
    for (&(idx1, idx2), &charge_term) in cpci {
        contrib.add_term(
            field.positions(),
            source_index(idx1)?,
            source_index(idx2)?,
            None,
            true,
            charge_term,
            diel_model,
            is_1_4,
        )?;
    }

    // RDKit✔️✔️:
    // RDKit✔️✔️:   return field;
    // RDKit✔️✔️: }
    field.add_contribution(Box::new(contrib));
    Ok(field)
}

/// Numeric input settings for distance-geometry contribution construction.
#[derive(Debug, Clone, Copy)]
pub struct DistanceGeometryForceFieldParams<'a> {
    pub weight_chiral: f64,
    pub weight_fourth_dimension: f64,
    pub basin_size_tolerance: f64,
    pub extra_weights: Option<&'a BTreeMap<(usize, usize), f64>>,
    pub fixed_pair_points: Option<&'a [bool]>,
    pub fixed_points: &'a [usize],
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ConformerOptimizerError {
    cause: ConformerOptimizerCause,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ConformerOptimizerCause {
    Input(&'static str),
    Construction(DistGeomForceFieldError),
}
impl std::fmt::Display for ConformerOptimizerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.cause {
            ConformerOptimizerCause::Input(message) => f.write_str(message),
            ConformerOptimizerCause::Construction(DistGeomForceFieldError::Kernel(e)) => e.fmt(f),
            ConformerOptimizerCause::Construction(DistGeomForceFieldError::Inversion(e)) => {
                e.fmt(f)
            }
            ConformerOptimizerCause::Construction(DistGeomForceFieldError::IndexOverflow(i)) => {
                write!(f, "index {i} exceeds the source unsigned index range")
            }
        }
    }
}
impl std::error::Error for ConformerOptimizerError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match &self.cause {
            ConformerOptimizerCause::Construction(DistGeomForceFieldError::Kernel(e)) => Some(e),
            ConformerOptimizerCause::Construction(DistGeomForceFieldError::Inversion(e)) => Some(e),
            _ => None,
        }
    }
}
impl From<DistGeomForceFieldError> for ConformerOptimizerError {
    fn from(e: DistGeomForceFieldError) -> Self {
        Self {
            cause: ConformerOptimizerCause::Construction(e),
        }
    }
}
impl From<ForceFieldKernelError> for ConformerOptimizerError {
    fn from(e: ForceFieldKernelError) -> Self {
        DistGeomForceFieldError::Kernel(e).into()
    }
}
fn input_error(message: &'static str) -> ConformerOptimizerError {
    ConformerOptimizerError {
        cause: ConformerOptimizerCause::Input(message),
    }
}

/// Evaluate the unique chiral-volume scalar on checked detached coordinate rows.
/// Only four stack points are projected; the private kernel and geometry types stay private.
pub fn calc_chiral_volume_rows(
    indices: [usize; 4],
    rows: &[Vec<f64>],
) -> Result<f64, ConformerOptimizerError> {
    // BEGIN RDKIT CPP FUNCTION calc_chiral_volume_rows (DistGeom/ChiralViolationContribs.cpp:38-59)
    // RDKit❗✔️: double calcChiralVolume(const unsigned int idx1, const unsigned int idx2,
    // RDKit❗✔️:                         const unsigned int idx3, const unsigned int idx4,
    // RDKit❗✔️:                         const RDGeom::PointPtrVect &pts) {
    // RDKit❗✔️:   // even if we are minimizing in higher dimension the chiral volume is
    // RDKit❗✔️:   // calculated using only the first 3 dimensions
    // RDKit❗✔️:   RDGeom::Point3D v1((*pts[idx1])[0] - (*pts[idx4])[0],
    // RDKit❗✔️:                      (*pts[idx1])[1] - (*pts[idx4])[1],
    // RDKit❗✔️:                      (*pts[idx1])[2] - (*pts[idx4])[2]);
    // RDKit❗✔️:
    // RDKit❗✔️:   RDGeom::Point3D v2((*pts[idx2])[0] - (*pts[idx4])[0],
    // RDKit❗✔️:                      (*pts[idx2])[1] - (*pts[idx4])[1],
    // RDKit❗✔️:                      (*pts[idx2])[2] - (*pts[idx4])[2]);
    // RDKit❗✔️:
    // RDKit❗✔️:   RDGeom::Point3D v3((*pts[idx3])[0] - (*pts[idx4])[0],
    // RDKit❗✔️:                      (*pts[idx3])[1] - (*pts[idx4])[1],
    // RDKit❗✔️:                      (*pts[idx3])[2] - (*pts[idx4])[2]);
    // RDKit❗✔️:
    // RDKit❗✔️:   RDGeom::Point3D v2xv3 = v2.crossProduct(v3);
    // RDKit❗✔️:
    // RDKit❗✔️:   double vol = v1.dotProduct(v2xv3);
    // RDKit❗✔️:   return vol;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION calc_chiral_volume_rows

    // RDKit❗✔️: The complete calcChiralVolume(PointPtrVect) source anchor and scalar
    // implementation remain in super::calc_chiral_volume_points; this boundary
    // only checks and projects the selected four source points without allocation.
    let mut points = [Point3 {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    }; 4];
    for (target, index) in points.iter_mut().zip(indices) {
        let row = rows
            .get(index)
            .ok_or_else(|| input_error("chiral point index out of range"))?;
        if row.len() < 3 {
            return Err(input_error("chiral point dimension less than three"));
        }
        *target = point(row);
    }
    Ok(super::calc_chiral_volume_points(0, 1, 2, 3, &points))
}

/// A conformer-only numeric optimizer borrowing detached coordinate rows.
/// The existing kernel, contributions, distance cache and positional bindings stay private.
/// Dropping this value releases the borrows; coordinates are updated directly by the one kernel.
pub struct ConformerOptimizer<'a> {
    field: ForceField<'a>,
}
impl<'a> ConformerOptimizer<'a> {
    pub fn distance_geometry(
        bounds: &impl DistanceBoundsRead,
        positions: &'a mut [Vec<f64>],
        chiral_sets: &[ChiralSetPtr],
        params: DistanceGeometryForceFieldParams<'_>,
    ) -> Result<Self, ConformerOptimizerError> {
        validate_rows(bounds, positions, None, params.fixed_points)?;
        if let Some(fixed) = params.fixed_pair_points {
            if fixed.len() < positions.len() {
                return Err(input_error("bad fixed point bitset"));
            }
        }
        if params.weight_chiral > 1.0e-8 {
            for set in chiral_sets {
                for i in [set.idx1, set.idx2, set.idx3, set.idx4] {
                    if i >= positions.len() {
                        return Err(input_error("chiral point index out of range"));
                    }
                }
            }
        }
        let field = construct_distgeom_forcefield(
            bounds,
            positions,
            chiral_sets,
            params.weight_chiral,
            params.weight_fourth_dimension,
            params.extra_weights,
            params.basin_size_tolerance,
            params.fixed_pair_points,
        )?;
        Self::initialize(field, params.fixed_points)
    }
    pub fn torsions(
        bounds: &impl DistanceBoundsRead,
        positions: &'a mut [Vec<f64>],
        details: &CrystalFFDetails,
        use_basic_knowledge: bool,
        cpci: Option<&BTreeMap<(usize, usize), f64>>,
        fixed_points: &[usize],
    ) -> Result<Self, ConformerOptimizerError> {
        validate_rows(bounds, positions, Some(3), fixed_points)?;
        validate_details(details, positions.len(), use_basic_knowledge, true)?;
        // The pinned plain ETDG branch ignores CPCI; only the basic-knowledge branch consumes it.
        let field = if use_basic_knowledge {
            if let Some(cpci) = cpci {
                construct_3d_forcefield_with_cpci(bounds, positions, details, cpci)?
            } else {
                construct_3d_forcefield(bounds, positions, details)?
            }
        } else {
            construct_plain_3d_forcefield(bounds, positions, details)?
        };
        Self::initialize(field, fixed_points)
    }
    pub fn improper(
        bounds: &impl DistanceBoundsRead,
        positions: &'a mut [Vec<f64>],
        details: &CrystalFFDetails,
        fixed_points: &[usize],
    ) -> Result<Self, ConformerOptimizerError> {
        validate_rows(bounds, positions, Some(3), fixed_points)?;
        validate_details(details, positions.len(), true, false)?;
        let field = construct_3d_improper_forcefield(bounds, positions, details)?;
        Self::initialize(field, fixed_points)
    }
    fn initialize(
        mut field: ForceField<'a>,
        fixed_points: &[usize],
    ) -> Result<Self, ConformerOptimizerError> {
        // BEGIN RDKIT CPP CALLS DGeomHelpers::EmbeddingOps initialization (Embedder.cpp)
        // RDKit✔️✔️: field->fixedPoints().push_back(v.first);
        // RDKit✔️✔️: field->initialize();
        // END RDKIT CPP CALLS DGeomHelpers::EmbeddingOps initialization
        // The source index type is unsigned; the existing kernel stores signed fixed-point indices.
        for &i in fixed_points {
            field.fixed_points_mut().push(i as i32);
        }
        field.initialize()?;
        Ok(Self { field })
    }
    pub fn energy(
        &mut self,
        contribution_energies: Option<&mut Vec<f64>>,
    ) -> Result<f64, ConformerOptimizerError> {
        self.field
            .calc_energy_current(contribution_energies)
            .map_err(Into::into)
    }
    pub fn minimize(
        &mut self,
        max_iterations: u32,
        force_tolerance: f64,
        energy_tolerance: f64,
    ) -> Result<i32, ConformerOptimizerError> {
        self.field
            .minimize(max_iterations, force_tolerance, energy_tolerance)
            .map_err(Into::into)
    }
}
fn validate_rows(
    bounds: &impl DistanceBoundsRead,
    positions: &[Vec<f64>],
    dimension: Option<usize>,
    fixed: &[usize],
) -> Result<(), ConformerOptimizerError> {
    if bounds.dimension() != positions.len() {
        return Err(input_error("Wrong size metric matrix"));
    }
    let Some(first) = positions.first() else {
        return Err(input_error("bad vector"));
    };
    let dimension = dimension.unwrap_or(first.len());
    if !(3..=4).contains(&dimension) || positions.iter().any(|p| p.len() != dimension) {
        return Err(input_error("inconsistent point dimension"));
    }
    for &i in fixed {
        if i >= positions.len() || i > i32::MAX as usize {
            return Err(input_error("bad fixed point index"));
        }
    }
    Ok(())
}
fn validate_details(
    details: &CrystalFFDetails,
    n: usize,
    basic: bool,
    all_terms: bool,
) -> Result<(), ConformerOptimizerError> {
    let check = |i: i32| {
        if i < 0 || i as usize >= n {
            Err(input_error("point index out of range"))
        } else {
            Ok(())
        }
    };
    if all_terms {
        if details.exp_torsion_atoms.len() != details.exp_torsion_angles.len() {
            return Err(input_error("torsion atom/angle size mismatch"));
        }
        for (atoms, (signs, constants)) in details
            .exp_torsion_atoms
            .iter()
            .zip(&details.exp_torsion_angles)
        {
            if atoms.len() < 4 || signs.len() < 6 || constants.len() < 6 {
                return Err(input_error("bad torsion term vector"));
            }
            for &i in &atoms[..4] {
                check(i)?;
            }
        }
        for &(i, j) in &details.bonds {
            check(i)?;
            check(j)?;
        }
        if !details.constrained_atoms.is_empty() && details.constrained_atoms.len() < n {
            return Err(input_error("bad constrained atom bitset"));
        }
    }
    if basic {
        for atoms in &details.improper_atoms {
            if atoms.len() < 6 {
                return Err(input_error("bad improper atom record"));
            }
            for &i in &atoms[..4] {
                check(i)?;
            }
        }
    }
    for angle in &details.angles {
        if angle.len() < if basic { 4 } else { 3 } {
            return Err(input_error("bad angle record"));
        }
        for &i in &angle[..3] {
            check(i)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::distgeom::ChiralSet;
    use crate::kernel::{EvaluationContext, ForceFieldContribution};
    use std::sync::Arc;
    // Literal test fixture bounds only, no production matrix or smoothing algorithm.
    struct BoundsMatrix {
        rows: Vec<Vec<f64>>,
    }
    impl BoundsMatrix {
        fn new(n: usize) -> Self {
            let mut rows = vec![vec![0.; n]; n];
            for i in 1..n {
                for j in 0..i {
                    rows[i][j] = 0.001;
                    rows[j][i] = 1000.0;
                }
            }
            Self { rows }
        }
        fn set_lower(&mut self, i: usize, j: usize, v: f64) -> Result<(), ()> {
            let (hi, lo) = (i.max(j), i.min(j));
            self.rows[hi][lo] = v;
            Ok(())
        }
        fn set_upper(&mut self, i: usize, j: usize, v: f64) -> Result<(), ()> {
            let (hi, lo) = (i.max(j), i.min(j));
            self.rows[lo][hi] = v;
            Ok(())
        }
    }
    impl DistanceBoundsRead for BoundsMatrix {
        fn dimension(&self) -> usize {
            self.rows.len()
        }
        fn get_lower(&self, i: usize, j: usize) -> f64 {
            self.rows[i.max(j)][i.min(j)]
        }
        fn get_upper(&self, i: usize, j: usize) -> f64 {
            self.rows[i.min(j)][i.max(j)]
        }
    }
    fn point3(x: f64, y: f64, z: f64) -> Point3 {
        Point3 { x, y, z }
    }
    trait FixtureRows {
        fn row(&self) -> Vec<f64>;
    }
    impl FixtureRows for Point3 {
        fn row(&self) -> Vec<f64> {
            vec![self.x, self.y, self.z]
        }
    }
    impl FixtureRows for Vec<f64> {
        fn row(&self) -> Vec<f64> {
            self.clone()
        }
    }
    fn rows_for_constructor(positions: &[impl FixtureRows]) -> Vec<Vec<f64>> {
        positions.iter().map(FixtureRows::row).collect()
    }
    fn fixture_field(rows: &mut [Vec<f64>], dimension: u32) -> ForceField<'_> {
        let mut field = ForceField::new(dimension);
        field
            .positions_mut()
            .extend(rows.iter_mut().map(Vec::as_mut_slice));
        field
    }
    fn field_points(field: &ForceField<'_>) -> Vec<Point3> {
        field.positions().iter().map(|p| point(p)).collect()
    }
    fn field_rows(field: &ForceField<'_>) -> Vec<Vec<f64>> {
        field.positions().iter().map(|p| p.to_vec()).collect()
    }
    fn evaluate_energy(contribution: &dyn ForceFieldContribution, pos: &[f64], n: usize) -> f64 {
        let mut cache = vec![-1.; n * (n + 1) / 2];
        let mut context = EvaluationContext::for_test(pos, &mut cache, n as u32);
        contribution
            .get_energy(&mut context)
            .expect("expected energy")
    }
    fn add_12_terms(
        ff: &mut ForceField<'_>,
        details: &CrystalFFDetails,
        pairs: &mut [bool],
        positions: &[impl Position3Read],
        force: f64,
        n: usize,
    ) -> Result<(), DistGeomForceFieldError> {
        let group = build_12_terms(ff, details, pairs, positions, force, n)?;
        if !group.empty() {
            ff.add_contribution(Box::new(group));
        }
        Ok(())
    }
    fn add_13_terms(
        ff: &mut ForceField<'_>,
        details: &CrystalFFDetails,
        pairs: &mut [bool],
        positions: &[impl Position3Read],
        force: f64,
        improper: &[bool],
        basic: bool,
        bounds: &impl DistanceBoundsRead,
        n: usize,
    ) -> Result<(), DistGeomForceFieldError> {
        let groups = build_13_terms(
            ff, details, pairs, positions, force, improper, basic, bounds, n,
        )?;
        if !groups.0.empty() {
            ff.add_contribution(Box::new(groups.0));
        }
        if !groups.1.empty() {
            ff.add_contribution(Box::new(groups.1));
        }
        Ok(())
    }
    fn add_long_range_distance_constraints(
        ff: &mut ForceField<'_>,
        details: &CrystalFFDetails,
        pairs: &[bool],
        positions: &[impl Position3Read],
        force: f64,
        bounds: &impl DistanceBoundsRead,
        n: usize,
    ) -> Result<(), DistGeomForceFieldError> {
        let group =
            build_long_range_distance_constraints(ff, details, pairs, positions, force, bounds, n)?;
        if !group.empty() {
            ff.add_contribution(Box::new(group));
        }
        Ok(())
    }
    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1.0e-12,
            "actual={actual} expected={expected}"
        );
    }

    #[test]
    fn construct_distgeom_forcefield_adds_distance_terms_for_basin_and_extra_weights() {
        let mut mmat = BoundsMatrix::new(3);
        mmat.set_lower(1, 0, 1.0).expect("set lower");
        mmat.set_upper(1, 0, 2.0).expect("set upper");
        mmat.set_lower(2, 0, 1.0).expect("set lower");
        mmat.set_upper(2, 0, 2.0).expect("set upper");
        mmat.set_lower(2, 1, 1.0).expect("set lower");
        mmat.set_upper(2, 1, 4.0).expect("set upper");
        let positions = vec![
            vec![0.0, 0.0, 0.0],
            vec![3.0, 0.0, 0.0],
            vec![0.0, 4.0, 0.0],
        ];
        let mut extra_weights = std::collections::BTreeMap::new();
        extra_weights.insert((2, 1), 2.0);
        let fixed_pts = vec![true, false, true];

        let mut field_bound_rows = rows_for_constructor(&positions);
        let mut field = construct_distgeom_forcefield(
            &mmat,
            &mut field_bound_rows,
            &[],
            0.0,
            0.0,
            Some(&extra_weights),
            1.0,
            Some(&fixed_pts),
        )
        .expect("source construction");
        field.initialize().expect("initialize");
        let mut contrib_energies = Vec::new();
        let energy = field
            .calc_energy_current(Some(&mut contrib_energies))
            .expect("energy");

        assert_eq!(field.dimension(), 3);
        assert_eq!(field_rows(&field), positions);
        assert_eq!(contrib_energies, vec![1.5625 + 0.6328125]);
        assert_eq!(energy, 2.1953125);
    }

    #[test]
    fn construct_distgeom_forcefield_adds_chiral_and_fourth_dimension_terms() {
        let mmat = BoundsMatrix::new(4);
        let positions = vec![
            vec![1.0, 0.0, 0.0, 1.0],
            vec![0.0, 1.0, 0.0, 2.0],
            vec![0.0, 0.0, 1.0, 3.0],
            vec![0.0, 0.0, 0.0, 4.0],
        ];
        let cset = Arc::new(ChiralSet::with_default_structure_flags(
            99, 0, 1, 2, 3, -0.1, 0.1,
        ));

        let mut field_bound_rows = rows_for_constructor(&positions);
        let mut field = construct_distgeom_forcefield(
            &mmat,
            &mut field_bound_rows,
            &[cset],
            2.0,
            0.5,
            None,
            0.0,
            None,
        )
        .expect("source construction");
        field.initialize().expect("initialize");
        let mut contrib_energies = Vec::new();
        let energy = field
            .calc_energy_current(Some(&mut contrib_energies))
            .expect("energy");
        let expected_chiral = 2.0 * (1.0_f64 - 0.1) * (1.0_f64 - 0.1);
        let expected_fourth = 0.5 * (1.0_f64 + 4.0 + 9.0 + 16.0);

        assert_eq!(field.dimension(), 4);
        assert_eq!(
            field_rows(&field),
            vec![
                vec![1.0, 0.0, 0.0, 1.0],
                vec![0.0, 1.0, 0.0, 2.0],
                vec![0.0, 0.0, 1.0, 3.0],
                vec![0.0, 0.0, 0.0, 4.0],
            ]
        );
        assert_eq!(contrib_energies, vec![expected_chiral, expected_fourth]);
        assert_eq!(energy, expected_chiral + expected_fourth);
    }

    #[test]
    fn construct_distgeom_forcefield_omits_empty_contrib_groups() {
        let mmat = BoundsMatrix::new(2);
        let positions = vec![vec![0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0]];

        let mut field_bound_rows = rows_for_constructor(&positions);
        let mut field = construct_distgeom_forcefield(
            &mmat,
            &mut field_bound_rows,
            &[],
            0.0,
            0.0,
            None,
            0.0,
            None,
        )
        .expect("source construction");
        field.initialize().expect("initialize");
        let mut contrib_energies = Vec::new();

        assert_eq!(
            field
                .calc_energy_current(Some(&mut contrib_energies))
                .expect("energy"),
            0.0
        );
        assert!(contrib_energies.is_empty());
    }

    #[test]
    #[should_panic]
    fn construct_distgeom_forcefield_rejects_size_mismatch() {
        let mmat = BoundsMatrix::new(2);
        let positions = vec![vec![0.0, 0.0, 0.0]];

        let mut __bound_rows = rows_for_constructor(&positions);
        let _ =
            construct_distgeom_forcefield(&mmat, &mut __bound_rows, &[], 0.0, 0.0, None, 0.0, None)
                .expect("source construction");
    }

    fn distgeom_improper_force_field_rows() -> Vec<Vec<f64>> {
        vec![
            vec![0.1, -0.2, 0.3],
            vec![1.0, 0.0, 0.1],
            vec![-0.2, 1.1, 0.2],
            vec![0.2, -0.1, 1.3],
        ]
    }

    fn flatten_force_field_positions(ff: &ForceField<'_>) -> Vec<f64> {
        ff.positions()
            .iter()
            .flat_map(|p| p.iter().copied())
            .collect()
    }

    #[test]
    fn distgeom_improper_torsion_terms_add_three_source_permutations_and_mark_center() {
        let mut ff_rows = distgeom_improper_force_field_rows();
        let mut ff = fixture_field(&mut ff_rows, 3);
        let improper_atoms = vec![vec![0, 1, 2, 3, 6, 1]];
        let mut is_improper_constrained = vec![false; 4];

        add_improper_torsion_terms(&mut ff, 2.0, &improper_atoms, &mut is_improper_constrained)
            .expect("source construction");
        ff.initialize().expect("initialize");
        let mut contrib_energies = Vec::new();
        let energy = ff
            .calc_energy_current(Some(&mut contrib_energies))
            .expect("energy");

        let mut expected_ff_rows = distgeom_improper_force_field_rows();
        let expected_ff = fixture_field(&mut expected_ff_rows, 3);
        let mut expected = InversionContribs::default();
        expected
            .add_contrib(expected_ff.positions(), 0, 1, 2, 3, 6, true, 2.0)
            .expect("expected inversion");
        expected
            .add_contrib(expected_ff.positions(), 0, 1, 3, 2, 6, true, 2.0)
            .expect("expected inversion");
        expected
            .add_contrib(expected_ff.positions(), 2, 1, 3, 0, 6, true, 2.0)
            .expect("expected inversion");
        let expected_energy = evaluate_energy(
            &expected,
            &flatten_force_field_positions(&expected_ff),
            expected_ff.positions().len(),
        );

        assert_eq!(is_improper_constrained, vec![false, true, false, false]);
        assert_eq!(contrib_energies.len(), 1);
        assert_close(contrib_energies[0], expected_energy);
        assert_close(energy, expected_energy);
    }

    #[test]
    fn distgeom_improper_torsion_terms_accumulate_multiple_centers_in_one_contrib_group() {
        let mut ff_rows = distgeom_improper_force_field_rows();
        let mut ff = fixture_field(&mut ff_rows, 3);
        let improper_atoms = vec![vec![0, 1, 2, 3, 6, 0], vec![3, 2, 1, 0, 15, 0]];
        let mut is_improper_constrained = vec![false; 4];

        add_improper_torsion_terms(&mut ff, 0.5, &improper_atoms, &mut is_improper_constrained)
            .expect("source construction");
        ff.initialize().expect("initialize");
        let mut contrib_energies = Vec::new();
        let energy = ff
            .calc_energy_current(Some(&mut contrib_energies))
            .expect("energy");

        let mut expected_ff_rows = distgeom_improper_force_field_rows();
        let expected_ff = fixture_field(&mut expected_ff_rows, 3);
        let mut expected = InversionContribs::default();
        expected
            .add_contrib(expected_ff.positions(), 0, 1, 2, 3, 6, false, 0.5)
            .expect("expected inversion");
        expected
            .add_contrib(expected_ff.positions(), 0, 1, 3, 2, 6, false, 0.5)
            .expect("expected inversion");
        expected
            .add_contrib(expected_ff.positions(), 2, 1, 3, 0, 6, false, 0.5)
            .expect("expected inversion");
        expected
            .add_contrib(expected_ff.positions(), 3, 2, 1, 0, 15, false, 0.5)
            .expect("expected inversion");
        expected
            .add_contrib(expected_ff.positions(), 3, 2, 0, 1, 15, false, 0.5)
            .expect("expected inversion");
        expected
            .add_contrib(expected_ff.positions(), 1, 2, 0, 3, 15, false, 0.5)
            .expect("expected inversion");
        let expected_energy = evaluate_energy(
            &expected,
            &flatten_force_field_positions(&expected_ff),
            expected_ff.positions().len(),
        );

        assert_eq!(is_improper_constrained, vec![false, true, true, false]);
        assert_eq!(contrib_energies.len(), 1);
        assert_close(contrib_energies[0], expected_energy);
        assert_close(energy, expected_energy);
    }

    #[test]
    fn distgeom_improper_torsion_terms_leave_force_field_unchanged_for_empty_input() {
        let mut ff_rows = distgeom_improper_force_field_rows();
        let mut ff = fixture_field(&mut ff_rows, 3);
        let mut is_improper_constrained = vec![false; 4];

        add_improper_torsion_terms(&mut ff, 1.0, &[], &mut is_improper_constrained)
            .expect("source construction");
        ff.initialize().expect("initialize");
        let mut contrib_energies = Vec::new();

        assert_eq!(
            ff.calc_energy_current(Some(&mut contrib_energies))
                .expect("energy"),
            0.0
        );
        assert!(contrib_energies.is_empty());
        assert_eq!(is_improper_constrained, vec![false; 4]);
    }

    #[test]
    #[should_panic]
    fn distgeom_improper_torsion_terms_reject_short_improper_atom_record() {
        let mut ff_rows = distgeom_improper_force_field_rows();
        let mut ff = fixture_field(&mut ff_rows, 3);
        let improper_atoms = vec![vec![0, 1, 2, 3, 6]];
        let mut is_improper_constrained = vec![false; 4];

        add_improper_torsion_terms(&mut ff, 1.0, &improper_atoms, &mut is_improper_constrained)
            .expect("source construction");
    }

    fn distgeom_experimental_torsion_force_field_rows() -> Vec<Vec<f64>> {
        vec![
            vec![0.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0],
            vec![1.0, 1.0, 0.0],
            vec![1.0, 1.0, 1.0],
            vec![2.0, 1.0, 1.0],
        ]
    }

    #[test]
    fn distgeom_experimental_torsion_terms_mark_ordered_endpoint_pairs_and_add_contribs() {
        let mut ff_rows = distgeom_experimental_torsion_force_field_rows();
        let mut ff = fixture_field(&mut ff_rows, 3);
        let details = CrystalFFDetails {
            exp_torsion_atoms: vec![vec![0, 1, 2, 3], vec![4, 3, 2, 1]],
            exp_torsion_angles: vec![
                (
                    vec![1, -1, 1, -1, 1, -1],
                    vec![0.10, 0.20, 0.30, 0.40, 0.50, 0.60],
                ),
                (
                    vec![-1, 1, -1, 1, -1, 1],
                    vec![0.05, 0.15, 0.25, 0.35, 0.45, 0.55],
                ),
            ],
            ..CrystalFFDetails::default()
        };
        let mut atom_pairs = vec![false; 25];

        add_experimental_torsion_terms(&mut ff, &details, &mut atom_pairs, 5)
            .expect("source construction");
        ff.initialize().expect("initialize");
        let mut contrib_energies = Vec::new();
        let energy = ff
            .calc_energy_current(Some(&mut contrib_energies))
            .expect("energy");

        let mut expected_ff_rows = distgeom_experimental_torsion_force_field_rows();
        let expected_ff = fixture_field(&mut expected_ff_rows, 3);
        let mut expected = TorsionAngleContribs::new(&expected_ff);
        expected.add_contrib(
            0,
            1,
            2,
            3,
            vec![0.10, 0.20, 0.30, 0.40, 0.50, 0.60],
            vec![1, -1, 1, -1, 1, -1],
        );
        expected.add_contrib(
            4,
            3,
            2,
            1,
            vec![0.05, 0.15, 0.25, 0.35, 0.45, 0.55],
            vec![-1, 1, -1, 1, -1, 1],
        );
        let expected_energy = evaluate_energy(
            &expected,
            &flatten_force_field_positions(&expected_ff),
            expected_ff.positions().len(),
        );

        assert!(atom_pairs[3]);
        assert!(atom_pairs[9]);
        assert_eq!(atom_pairs.iter().filter(|&&set| set).count(), 2);
        assert_eq!(contrib_energies.len(), 1);
        assert_close(contrib_energies[0], expected_energy);
        assert_close(energy, expected_energy);
    }

    #[test]
    fn distgeom_experimental_torsion_terms_leave_force_field_unchanged_for_empty_input() {
        let mut ff_rows = distgeom_experimental_torsion_force_field_rows();
        let mut ff = fixture_field(&mut ff_rows, 3);
        let details = CrystalFFDetails::default();
        let mut atom_pairs = vec![false; 25];

        add_experimental_torsion_terms(&mut ff, &details, &mut atom_pairs, 5)
            .expect("source construction");
        ff.initialize().expect("initialize");
        let mut contrib_energies = Vec::new();

        assert_eq!(
            ff.calc_energy_current(Some(&mut contrib_energies))
                .expect("energy"),
            0.0
        );
        assert!(contrib_energies.is_empty());
        assert_eq!(atom_pairs, vec![false; 25]);
    }

    #[test]
    #[should_panic]
    fn distgeom_experimental_torsion_terms_reject_short_torsion_atom_record() {
        let mut ff_rows = distgeom_experimental_torsion_force_field_rows();
        let mut ff = fixture_field(&mut ff_rows, 3);
        let details = CrystalFFDetails {
            exp_torsion_atoms: vec![vec![0, 1, 2]],
            exp_torsion_angles: vec![(vec![1; 6], vec![0.1; 6])],
            ..CrystalFFDetails::default()
        };
        let mut atom_pairs = vec![false; 25];

        add_experimental_torsion_terms(&mut ff, &details, &mut atom_pairs, 5)
            .expect("source construction");
    }

    #[test]
    fn distgeom_add12_terms_mark_ordered_bond_pairs_and_add_distance_contribs() {
        let mut ff_rows = distgeom_experimental_torsion_force_field_rows();
        let mut ff = fixture_field(&mut ff_rows, 3);
        let positions = vec![
            point3(0.0, 0.0, 0.0),
            point3(2.0, 0.0, 0.0),
            point3(1.0, 1.0, 0.0),
            point3(1.0, 1.0, 1.0),
            point3(2.0, 3.0, 0.0),
        ];
        let details = CrystalFFDetails {
            bonds: vec![(0, 1), (4, 1)],
            ..CrystalFFDetails::default()
        };
        let mut atom_pairs = vec![false; 25];

        add_12_terms(&mut ff, &details, &mut atom_pairs, &positions, 7.0, 5)
            .expect("source construction");
        ff.initialize().expect("initialize");
        let mut contrib_energies = Vec::new();
        let energy = ff
            .calc_energy_current(Some(&mut contrib_energies))
            .expect("energy");

        let mut expected_ff_rows = distgeom_experimental_torsion_force_field_rows();
        let mut expected_ff = fixture_field(&mut expected_ff_rows, 3);
        expected_ff.initialize().expect("initialize");
        let mut expected = DistanceConstraintContribs::new(&expected_ff);
        expected
            .add_contrib(&expected_ff, 0, 1, 1.99, 2.01, 7.0)
            .expect("expected constraint");
        let d41 = Point3::difference(&positions[4], &positions[1]).length();
        expected
            .add_contrib(
                &expected_ff,
                4,
                1,
                d41 - KNOWN_DIST_TOL,
                d41 + KNOWN_DIST_TOL,
                7.0,
            )
            .expect("expected constraint");
        let expected_energy = evaluate_energy(
            &expected,
            &flatten_force_field_positions(&expected_ff),
            expected_ff.positions().len(),
        );

        assert!(atom_pairs[1]);
        assert!(atom_pairs[9]);
        assert_eq!(atom_pairs.iter().filter(|&&set| set).count(), 2);
        assert_eq!(contrib_energies.len(), 1);
        assert_close(contrib_energies[0], expected_energy);
        assert_close(energy, expected_energy);
    }

    #[test]
    fn distgeom_add12_terms_leave_force_field_unchanged_for_empty_input() {
        let mut ff_rows = distgeom_experimental_torsion_force_field_rows();
        let mut ff = fixture_field(&mut ff_rows, 3);
        let positions = field_points(&ff);
        let details = CrystalFFDetails::default();
        let mut atom_pairs = vec![false; 25];

        add_12_terms(&mut ff, &details, &mut atom_pairs, &positions, 3.0, 5)
            .expect("source construction");
        ff.initialize().expect("initialize");
        let mut contrib_energies = Vec::new();

        assert_eq!(
            ff.calc_energy_current(Some(&mut contrib_energies))
                .expect("energy"),
            0.0
        );
        assert!(contrib_energies.is_empty());
        assert_eq!(atom_pairs, vec![false; 25]);
    }

    #[test]
    #[should_panic]
    fn distgeom_add12_terms_reject_out_of_range_position_index() {
        let mut ff_rows = distgeom_experimental_torsion_force_field_rows();
        let mut ff = fixture_field(&mut ff_rows, 3);
        let positions = field_points(&ff);
        let details = CrystalFFDetails {
            bonds: vec![(0, 5)],
            ..CrystalFFDetails::default()
        };
        let mut atom_pairs = vec![false; 30];

        add_12_terms(&mut ff, &details, &mut atom_pairs, &positions, 3.0, 6)
            .expect("source construction");
    }

    #[test]
    fn distgeom_add13_terms_add_angle_improper_bounds_and_current_distance_contribs() {
        let mut ff_rows = distgeom_experimental_torsion_force_field_rows();
        let mut ff = fixture_field(&mut ff_rows, 3);
        let positions = vec![
            point3(0.0, 0.0, 0.0),
            point3(2.0, 0.0, 0.0),
            point3(1.0, 1.0, 0.0),
            point3(1.0, 1.0, 1.0),
            point3(2.0, 3.0, 0.0),
        ];
        let details = CrystalFFDetails {
            angles: vec![vec![0, 1, 2, 1], vec![0, 1, 3, 0], vec![4, 3, 1, 0]],
            ..CrystalFFDetails::default()
        };
        let mut atom_pairs = vec![false; 25];
        let mut mmat = BoundsMatrix::new(5);
        mmat.set_lower(0, 3, 1.5).expect("set lower");
        mmat.set_upper(0, 3, 1.6).expect("set upper");
        let is_improper_constrained = vec![false, true, false, false, false];

        add_13_terms(
            &mut ff,
            &details,
            &mut atom_pairs,
            &positions,
            4.0,
            &is_improper_constrained,
            true,
            &mmat,
            5,
        )
        .expect("source construction");
        ff.initialize().expect("initialize");
        let mut contrib_energies = Vec::new();
        let energy = ff
            .calc_energy_current(Some(&mut contrib_energies))
            .expect("energy");

        let mut expected_ff_rows = distgeom_experimental_torsion_force_field_rows();
        let mut expected_ff = fixture_field(&mut expected_ff_rows, 3);
        expected_ff.initialize().expect("initialize");
        let mut expected_angle = AngleConstraintContribs::new(&expected_ff);
        expected_angle
            .add_contrib(&expected_ff, 0, 1, 2, 179.0, 180.0, 1.0)
            .expect("expected constraint");
        let expected_angle_energy = evaluate_energy(
            &expected_angle,
            &flatten_force_field_positions(&expected_ff),
            expected_ff.positions().len(),
        );

        let mut expected_dist = DistanceConstraintContribs::new(&expected_ff);
        expected_dist
            .add_contrib(&expected_ff, 0, 3, 1.5, 1.6, 4.0)
            .expect("expected constraint");
        let d41 = Point3::difference(&positions[4], &positions[1]).length();
        expected_dist
            .add_contrib(
                &expected_ff,
                4,
                1,
                d41 - KNOWN_DIST_TOL,
                d41 + KNOWN_DIST_TOL,
                4.0,
            )
            .expect("expected constraint");
        let expected_dist_energy = evaluate_energy(
            &expected_dist,
            &flatten_force_field_positions(&expected_ff),
            expected_ff.positions().len(),
        );

        assert!(atom_pairs[2]);
        assert!(atom_pairs[3]);
        assert!(atom_pairs[9]);
        assert_eq!(atom_pairs.iter().filter(|&&set| set).count(), 3);
        assert_eq!(contrib_energies.len(), 2);
        assert_close(contrib_energies[0], expected_angle_energy);
        assert_close(contrib_energies[1], expected_dist_energy);
        assert_close(energy, expected_angle_energy + expected_dist_energy);
    }

    #[test]
    fn distgeom_add13_terms_leave_force_field_unchanged_for_empty_input() {
        let mut ff_rows = distgeom_experimental_torsion_force_field_rows();
        let mut ff = fixture_field(&mut ff_rows, 3);
        let positions = field_points(&ff);
        let details = CrystalFFDetails::default();
        let mut atom_pairs = vec![false; 25];
        let mmat = BoundsMatrix::new(5);
        let is_improper_constrained = vec![false; 5];

        add_13_terms(
            &mut ff,
            &details,
            &mut atom_pairs,
            &positions,
            4.0,
            &is_improper_constrained,
            true,
            &mmat,
            5,
        )
        .expect("source construction");
        ff.initialize().expect("initialize");
        let mut contrib_energies = Vec::new();

        assert_eq!(
            ff.calc_energy_current(Some(&mut contrib_energies))
                .expect("energy"),
            0.0
        );
        assert!(contrib_energies.is_empty());
        assert_eq!(atom_pairs, vec![false; 25]);
    }

    #[test]
    #[should_panic]
    fn distgeom_add13_terms_reject_short_angle_record() {
        let mut ff_rows = distgeom_experimental_torsion_force_field_rows();
        let mut ff = fixture_field(&mut ff_rows, 3);
        let positions = field_points(&ff);
        let details = CrystalFFDetails {
            angles: vec![vec![0, 1, 2]],
            ..CrystalFFDetails::default()
        };
        let mut atom_pairs = vec![false; 25];
        let mmat = BoundsMatrix::new(5);
        let is_improper_constrained = vec![false; 5];

        add_13_terms(
            &mut ff,
            &details,
            &mut atom_pairs,
            &positions,
            4.0,
            &is_improper_constrained,
            true,
            &mmat,
            5,
        )
        .expect("source construction");
    }

    #[test]
    fn distgeom_long_range_distance_constraints_skip_atom_pairs_and_use_bounds_or_tight_constraints()
     {
        let mut ff_rows = distgeom_experimental_torsion_force_field_rows();
        let mut ff = fixture_field(&mut ff_rows, 3);
        let positions = vec![
            point3(0.0, 0.0, 0.0),
            point3(2.0, 0.0, 0.0),
            point3(1.0, 1.0, 0.0),
            point3(1.0, 1.0, 1.0),
        ];
        let mut mmat = BoundsMatrix::new(4);
        for i in 1..4 {
            for j in 0..i {
                mmat.set_lower(i, j, 0.5 + i as f64 + j as f64 * 0.1)
                    .expect("set lower");
                mmat.set_upper(i, j, 2.5 + i as f64 + j as f64 * 0.1)
                    .expect("set upper");
            }
        }
        let details = CrystalFFDetails {
            bounds_mat_force_scaling: 2.5,
            constrained_atoms: vec![true, true, false, true],
            ..CrystalFFDetails::default()
        };
        let mut atom_pairs = vec![false; 16];
        atom_pairs[2] = true;

        add_long_range_distance_constraints(
            &mut ff,
            &details,
            &atom_pairs,
            &positions,
            8.0,
            &mmat,
            4,
        )
        .expect("source construction");
        ff.initialize().expect("initialize");
        let mut contrib_energies = Vec::new();
        let energy = ff
            .calc_energy_current(Some(&mut contrib_energies))
            .expect("energy");

        let mut expected_ff_rows = distgeom_experimental_torsion_force_field_rows();
        let mut expected_ff = fixture_field(&mut expected_ff_rows, 3);
        expected_ff.initialize().expect("initialize");
        let mut expected = DistanceConstraintContribs::new(&expected_ff);
        let d10 = Point3::difference(&positions[1], &positions[0]).length();
        expected
            .add_contrib(
                &expected_ff,
                1,
                0,
                d10 - KNOWN_DIST_TOL,
                d10 + KNOWN_DIST_TOL,
                8.0,
            )
            .expect("expected constraint");
        expected
            .add_contrib(
                &expected_ff,
                2,
                1,
                mmat.get_lower(2, 1),
                mmat.get_upper(2, 1),
                25.0,
            )
            .expect("expected constraint");
        let d30 = Point3::difference(&positions[3], &positions[0]).length();
        expected
            .add_contrib(
                &expected_ff,
                3,
                0,
                d30 - KNOWN_DIST_TOL,
                d30 + KNOWN_DIST_TOL,
                8.0,
            )
            .expect("expected constraint");
        let d31 = Point3::difference(&positions[3], &positions[1]).length();
        expected
            .add_contrib(
                &expected_ff,
                3,
                1,
                d31 - KNOWN_DIST_TOL,
                d31 + KNOWN_DIST_TOL,
                8.0,
            )
            .expect("expected constraint");
        expected
            .add_contrib(
                &expected_ff,
                3,
                2,
                mmat.get_lower(3, 2),
                mmat.get_upper(3, 2),
                25.0,
            )
            .expect("expected constraint");
        let expected_energy = evaluate_energy(
            &expected,
            &flatten_force_field_positions(&expected_ff),
            expected_ff.positions().len(),
        );

        assert_eq!(contrib_energies.len(), 1);
        assert_close(contrib_energies[0], expected_energy);
        assert_close(energy, expected_energy);
    }

    #[test]
    fn distgeom_long_range_distance_constraints_leave_force_field_unchanged_when_all_pairs_present()
    {
        let mut ff_rows = distgeom_experimental_torsion_force_field_rows();
        let mut ff = fixture_field(&mut ff_rows, 3);
        let positions = field_points(&ff);
        let details = CrystalFFDetails {
            bounds_mat_force_scaling: 2.5,
            ..CrystalFFDetails::default()
        };
        let atom_pairs = vec![true; 25];
        let mmat = BoundsMatrix::new(5);

        add_long_range_distance_constraints(
            &mut ff,
            &details,
            &atom_pairs,
            &positions,
            8.0,
            &mmat,
            5,
        )
        .expect("source construction");
        ff.initialize().expect("initialize");
        let mut contrib_energies = Vec::new();

        assert_eq!(
            ff.calc_energy_current(Some(&mut contrib_energies))
                .expect("energy"),
            0.0
        );
        assert!(contrib_energies.is_empty());
    }

    #[test]
    #[should_panic]
    fn distgeom_long_range_distance_constraints_reject_short_atom_pair_bitset() {
        let mut ff_rows = distgeom_experimental_torsion_force_field_rows();
        let mut ff = fixture_field(&mut ff_rows, 3);
        let positions = field_points(&ff);
        let details = CrystalFFDetails::default();
        let atom_pairs = vec![false; 3];
        let mmat = BoundsMatrix::new(5);

        add_long_range_distance_constraints(
            &mut ff,
            &details,
            &atom_pairs,
            &positions,
            8.0,
            &mmat,
            5,
        )
        .expect("source construction");
    }

    fn construct_3d_forcefield_positions() -> Vec<Point3> {
        vec![
            point3(0.0, 0.0, 0.0),
            point3(1.2, 0.0, 0.0),
            point3(1.2, 1.1, 0.0),
            point3(1.2, 1.1, 1.3),
            point3(2.4, 1.1, 1.3),
        ]
    }

    fn construct_3d_forcefield_bounds_matrix(num_atoms: usize) -> BoundsMatrix {
        let mut mmat = BoundsMatrix::new(num_atoms);
        for i in 1..num_atoms {
            for j in 0..i {
                mmat.set_lower(i, j, 0.75 + i as f64 * 0.1 + j as f64 * 0.01)
                    .expect("set lower");
                mmat.set_upper(i, j, 3.25 + i as f64 * 0.1 + j as f64 * 0.01)
                    .expect("set upper");
            }
        }
        mmat
    }

    fn construct_3d_forcefield_details() -> CrystalFFDetails {
        CrystalFFDetails {
            exp_torsion_atoms: vec![vec![0, 1, 2, 3]],
            exp_torsion_angles: vec![(
                vec![1, -1, 1, -1, 1, -1],
                vec![0.10, 0.20, 0.30, 0.40, 0.50, 0.60],
            )],
            improper_atoms: vec![vec![1, 2, 3, 4, 6, 0]],
            bonds: vec![(0, 1), (3, 4)],
            angles: vec![vec![0, 1, 2, 1], vec![1, 2, 4, 0]],
            bounds_mat_force_scaling: 1.7,
            constrained_atoms: vec![true, true, false, true, true],
            ..CrystalFFDetails::default()
        }
    }

    #[test]
    fn construct_3d_forcefield_reproduces_source_helper_sequence_and_positions() {
        let positions = construct_3d_forcefield_positions();
        let mmat = construct_3d_forcefield_bounds_matrix(positions.len());
        let details = construct_3d_forcefield_details();

        let mut field_bound_rows = rows_for_constructor(&positions);
        let mut field = construct_3d_forcefield(&mmat, &mut field_bound_rows, &details)
            .expect("source construction");
        field.initialize().expect("initialize");
        let mut contrib_energies = Vec::new();
        let energy = field
            .calc_energy_current(Some(&mut contrib_energies))
            .expect("energy");

        let mut expected_bound_rows = rows_for_constructor(&positions);
        let mut expected = fixture_field(&mut expected_bound_rows, 3);
        let mut atom_pairs = vec![false; positions.len() * positions.len()];
        let mut is_improper_constrained = vec![false; positions.len()];
        add_experimental_torsion_terms(&mut expected, &details, &mut atom_pairs, positions.len())
            .expect("source construction");
        add_improper_torsion_terms(
            &mut expected,
            10.0,
            &details.improper_atoms,
            &mut is_improper_constrained,
        )
        .expect("source construction");
        add_12_terms(
            &mut expected,
            &details,
            &mut atom_pairs,
            &positions,
            KNOWN_DIST_FORCE_CONSTANT,
            positions.len(),
        )
        .expect("source construction");
        add_13_terms(
            &mut expected,
            &details,
            &mut atom_pairs,
            &positions,
            KNOWN_DIST_FORCE_CONSTANT,
            &is_improper_constrained,
            true,
            &mmat,
            positions.len(),
        )
        .expect("source construction");
        add_long_range_distance_constraints(
            &mut expected,
            &details,
            &atom_pairs,
            &positions,
            KNOWN_DIST_FORCE_CONSTANT,
            &mmat,
            positions.len(),
        )
        .expect("source construction");
        expected.initialize().expect("initialize");
        let mut expected_contrib_energies = Vec::new();
        let expected_energy = expected
            .calc_energy_current(Some(&mut expected_contrib_energies))
            .expect("energy");

        assert_eq!(field.dimension(), 3);
        assert_eq!(field_points(&field), positions);
        assert_eq!(contrib_energies.len(), expected_contrib_energies.len());
        assert_eq!(contrib_energies.len(), 6);
        for (observed, expected) in contrib_energies
            .iter()
            .zip(expected_contrib_energies.iter())
        {
            assert_close(*observed, *expected);
        }
        assert_close(energy, expected_energy);
    }

    #[test]
    #[should_panic]
    fn construct_3d_forcefield_rejects_bounds_position_size_mismatch() {
        let positions = construct_3d_forcefield_positions();
        let mmat = BoundsMatrix::new(positions.len() + 1);
        let details = CrystalFFDetails::default();

        let mut __bound_rows = rows_for_constructor(&positions);
        let _ = construct_3d_forcefield(&mmat, &mut __bound_rows, &details)
            .expect("source construction");
    }

    #[test]
    #[should_panic]
    fn construct_3d_forcefield_rejects_torsion_atom_angle_size_mismatch() {
        let positions = construct_3d_forcefield_positions();
        let mmat = BoundsMatrix::new(positions.len());
        let details = CrystalFFDetails {
            exp_torsion_atoms: vec![vec![0, 1, 2, 3]],
            exp_torsion_angles: Vec::new(),
            ..CrystalFFDetails::default()
        };

        let mut __bound_rows = rows_for_constructor(&positions);
        let _ = construct_3d_forcefield(&mmat, &mut __bound_rows, &details)
            .expect("source construction");
    }

    #[test]
    fn construct_3d_forcefield_with_cpci_appends_electrostatic_terms_after_base_field() {
        let positions = vec![
            point3(0.0, 0.0, 0.0),
            point3(1.0, 0.0, 0.0),
            point3(0.0, 2.0, 0.0),
        ];
        let mmat = construct_3d_forcefield_bounds_matrix(positions.len());
        let details = CrystalFFDetails::default();
        let mut base_bound_rows = rows_for_constructor(&positions);
        let base = construct_3d_forcefield(&mmat, &mut base_bound_rows, &details)
            .expect("source construction");
        let cpci = std::collections::BTreeMap::from([((0, 1), 0.5), ((1, 2), -0.25)]);

        let mut field_bound_rows = rows_for_constructor(&positions);
        let mut field =
            construct_3d_forcefield_with_cpci(&mmat, &mut field_bound_rows, &details, &cpci)
                .expect("source construction");
        field.initialize().expect("initialize");
        let mut contrib_energies = Vec::new();
        let energy = field
            .calc_energy_current(Some(&mut contrib_energies))
            .expect("energy");

        let mut expected_cpci_field_bound_rows = rows_for_constructor(&positions);
        let mut expected_cpci_field = fixture_field(&mut expected_cpci_field_bound_rows, 3);
        expected_cpci_field.initialize().expect("initialize");
        let mut expected_cpci = NonbondedContrib::new(&expected_cpci_field);
        expected_cpci
            .add_term(
                expected_cpci_field.positions(),
                0,
                1,
                None,
                true,
                0.5,
                1,
                false,
            )
            .expect("expected CPCI");
        expected_cpci
            .add_term(
                expected_cpci_field.positions(),
                1,
                2,
                None,
                true,
                -0.25,
                1,
                false,
            )
            .expect("expected CPCI");
        let expected_cpci_energy = evaluate_energy(
            &expected_cpci,
            &flatten_force_field_positions(&expected_cpci_field),
            expected_cpci_field.positions().len(),
        );

        let mut base = base;
        base.initialize().expect("initialize");
        let mut base_contrib_energies = Vec::new();
        let base_energy = base
            .calc_energy_current(Some(&mut base_contrib_energies))
            .expect("energy");

        assert_eq!(contrib_energies.len(), base_contrib_energies.len() + 1);
        for (observed, expected) in contrib_energies
            .iter()
            .take(base_contrib_energies.len())
            .zip(base_contrib_energies.iter())
        {
            assert_close(*observed, *expected);
        }
        assert_close(
            *contrib_energies.last().expect("CPCI contribution"),
            expected_cpci_energy,
        );
        assert_close(energy, base_energy + expected_cpci_energy);
    }

    #[test]
    fn construct_plain_3d_forcefield_reproduces_source_helper_sequence_without_improper_knowledge()
    {
        let positions = construct_3d_forcefield_positions();
        let mmat = construct_3d_forcefield_bounds_matrix(positions.len());
        let details = construct_3d_forcefield_details();

        let mut field_bound_rows = rows_for_constructor(&positions);
        let mut field = construct_plain_3d_forcefield(&mmat, &mut field_bound_rows, &details)
            .expect("source construction");
        field.initialize().expect("initialize");
        let mut contrib_energies = Vec::new();
        let energy = field
            .calc_energy_current(Some(&mut contrib_energies))
            .expect("energy");

        let mut expected_bound_rows = rows_for_constructor(&positions);
        let mut expected = fixture_field(&mut expected_bound_rows, 3);
        let mut atom_pairs = vec![false; positions.len() * positions.len()];
        let is_improper_constrained = vec![false; positions.len()];
        add_experimental_torsion_terms(&mut expected, &details, &mut atom_pairs, positions.len())
            .expect("source construction");
        add_12_terms(
            &mut expected,
            &details,
            &mut atom_pairs,
            &positions,
            KNOWN_DIST_FORCE_CONSTANT,
            positions.len(),
        )
        .expect("source construction");
        add_13_terms(
            &mut expected,
            &details,
            &mut atom_pairs,
            &positions,
            KNOWN_DIST_FORCE_CONSTANT,
            &is_improper_constrained,
            false,
            &mmat,
            positions.len(),
        )
        .expect("source construction");
        add_long_range_distance_constraints(
            &mut expected,
            &details,
            &atom_pairs,
            &positions,
            KNOWN_DIST_FORCE_CONSTANT,
            &mmat,
            positions.len(),
        )
        .expect("source construction");
        expected.initialize().expect("initialize");
        let mut expected_contrib_energies = Vec::new();
        let expected_energy = expected
            .calc_energy_current(Some(&mut expected_contrib_energies))
            .expect("energy");

        assert_eq!(field.dimension(), 3);
        assert_eq!(field_points(&field), positions);
        assert_eq!(contrib_energies.len(), expected_contrib_energies.len());
        assert_eq!(contrib_energies.len(), 4);
        for (observed, expected) in contrib_energies
            .iter()
            .zip(expected_contrib_energies.iter())
        {
            assert_close(*observed, *expected);
        }
        assert_close(energy, expected_energy);
    }

    #[test]
    #[should_panic]
    fn construct_plain_3d_forcefield_rejects_bounds_position_size_mismatch() {
        let positions = construct_3d_forcefield_positions();
        let mmat = BoundsMatrix::new(positions.len() + 1);
        let details = CrystalFFDetails::default();

        let mut __bound_rows = rows_for_constructor(&positions);
        let _ = construct_plain_3d_forcefield(&mmat, &mut __bound_rows, &details)
            .expect("source construction");
    }

    #[test]
    #[should_panic]
    fn construct_plain_3d_forcefield_rejects_torsion_atom_angle_size_mismatch() {
        let positions = construct_3d_forcefield_positions();
        let mmat = BoundsMatrix::new(positions.len());
        let details = CrystalFFDetails {
            exp_torsion_atoms: vec![vec![0, 1, 2, 3]],
            exp_torsion_angles: Vec::new(),
            ..CrystalFFDetails::default()
        };

        let mut __bound_rows = rows_for_constructor(&positions);
        let _ = construct_plain_3d_forcefield(&mmat, &mut __bound_rows, &details)
            .expect("source construction");
    }

    #[test]
    fn construct_3d_improper_forcefield_reproduces_parts_overload_terms() {
        let positions = construct_3d_forcefield_positions();
        let mmat = construct_3d_forcefield_bounds_matrix(positions.len());
        let improper_atoms = vec![vec![1, 2, 3, 4, 6, 0]];
        let angles = vec![vec![0, 1, 2, 1], vec![1, 2, 4, 0], vec![4, 3, 2, 1]];
        let atom_nums = vec![6, 7, 8, 9, 16];

        let mut field_bound_rows = rows_for_constructor(&positions);
        let mut field = construct_3d_improper_forcefield_from_parts(
            &mmat,
            &mut field_bound_rows,
            &improper_atoms,
            &angles,
            &atom_nums,
        )
        .expect("source construction");
        field.initialize().expect("initialize");
        let mut contrib_energies = Vec::new();
        let energy = field
            .calc_energy_current(Some(&mut contrib_energies))
            .expect("energy");

        let mut expected_bound_rows = rows_for_constructor(&positions);
        let mut expected = fixture_field(&mut expected_bound_rows, 3);
        let mut is_improper_constrained = vec![false; positions.len()];
        add_improper_torsion_terms(
            &mut expected,
            10.0,
            &improper_atoms,
            &mut is_improper_constrained,
        )
        .expect("source construction");
        let mut expected_angles = AngleConstraintContribs::new(&expected);
        expected_angles
            .add_contrib(&expected, 0, 1, 2, 179.0, 180.0, 10.0)
            .expect("expected constraint");
        expected_angles
            .add_contrib(&expected, 4, 3, 2, 179.0, 180.0, 10.0)
            .expect("expected constraint");
        expected.add_contribution(Box::new(expected_angles));
        expected.initialize().expect("initialize");
        let mut expected_contrib_energies = Vec::new();
        let expected_energy = expected
            .calc_energy_current(Some(&mut expected_contrib_energies))
            .expect("energy");

        assert_eq!(field.dimension(), 3);
        assert_eq!(field_points(&field), positions);
        assert_eq!(contrib_energies.len(), expected_contrib_energies.len());
        assert_eq!(contrib_energies.len(), 2);
        for (observed, expected) in contrib_energies
            .iter()
            .zip(expected_contrib_energies.iter())
        {
            assert_close(*observed, *expected);
        }
        assert_close(energy, expected_energy);
    }

    #[test]
    fn construct_3d_improper_forcefield_details_overload_delegates_to_parts_and_ignores_atom_nums()
    {
        let positions = construct_3d_forcefield_positions();
        let mmat = construct_3d_forcefield_bounds_matrix(positions.len());
        let details = CrystalFFDetails {
            improper_atoms: vec![vec![1, 2, 3, 4, 6, 0]],
            angles: vec![vec![0, 1, 2, 1], vec![1, 2, 4, 0], vec![4, 3, 2, 1]],
            atom_nums: vec![6, 7, 8, 9, 16],
            ..CrystalFFDetails::default()
        };

        let mut field_bound_rows = rows_for_constructor(&positions);
        let mut field = construct_3d_improper_forcefield(&mmat, &mut field_bound_rows, &details)
            .expect("source construction");
        let mut expected_bound_rows = rows_for_constructor(&positions);
        let mut expected = construct_3d_improper_forcefield_from_parts(
            &mmat,
            &mut expected_bound_rows,
            &details.improper_atoms,
            &details.angles,
            &[999, 998],
        )
        .expect("source construction");
        field.initialize().expect("initialize");
        expected.initialize().expect("initialize");
        let mut contrib_energies = Vec::new();
        let mut expected_contrib_energies = Vec::new();

        let energy = field
            .calc_energy_current(Some(&mut contrib_energies))
            .expect("energy");
        let expected_energy = expected
            .calc_energy_current(Some(&mut expected_contrib_energies))
            .expect("energy");

        assert_eq!(contrib_energies.len(), expected_contrib_energies.len());
        for (observed, expected) in contrib_energies
            .iter()
            .zip(expected_contrib_energies.iter())
        {
            assert_close(*observed, *expected);
        }
        assert_close(energy, expected_energy);
    }

    #[test]
    #[should_panic]
    fn construct_3d_improper_forcefield_rejects_bounds_position_size_mismatch() {
        let positions = construct_3d_forcefield_positions();
        let mmat = BoundsMatrix::new(positions.len() + 1);
        let details = CrystalFFDetails::default();

        let mut __bound_rows = rows_for_constructor(&positions);
        let _ = construct_3d_improper_forcefield(&mmat, &mut __bound_rows, &details)
            .expect("source construction");
    }

    #[test]
    #[should_panic]
    fn construct_3d_improper_forcefield_rejects_short_angle_record() {
        let positions = construct_3d_forcefield_positions();
        let mmat = BoundsMatrix::new(positions.len());
        let details = CrystalFFDetails {
            angles: vec![vec![0, 1, 2]],
            ..CrystalFFDetails::default()
        };

        let mut __bound_rows = rows_for_constructor(&positions);
        let _ = construct_3d_improper_forcefield(&mmat, &mut __bound_rows, &details)
            .expect("source construction");
    }

    #[test]
    fn source_first_minimization_exact_two_point_coordinates_via_checked_optimizer() {
        let mut mmat = BoundsMatrix::new(2);
        mmat.set_lower(0, 1, 1.0).expect("set lower");
        mmat.set_upper(0, 1, 1.0).expect("set upper");
        let mut positions = vec![vec![0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0]];
        {
            let mut optimizer = ConformerOptimizer::distance_geometry(
                &mmat,
                &mut positions,
                &[],
                DistanceGeometryForceFieldParams {
                    weight_chiral: 1.0,
                    weight_fourth_dimension: 0.1,
                    basin_size_tolerance: 5.0,
                    extra_weights: None,
                    fixed_pair_points: None,
                    fixed_points: &[],
                },
            )
            .expect("source optimizer");
            assert_eq!(optimizer.energy(None).expect("energy"), 0.0);
        }
        assert_eq!(positions, vec![vec![0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0]]);
    }
    #[test]
    fn source_fourth_dimension_random_coordinate_map_stays_fixed_via_checked_optimizer() {
        let mut mmat = BoundsMatrix::new(2);
        mmat.set_lower(0, 1, 1.0).expect("set lower");
        mmat.set_upper(0, 1, 1.0).expect("set upper");
        let mut positions = vec![vec![7.0, 8.0, 9.0, 0.0], vec![10.0, 8.0, 9.0, 3.0]];
        {
            let mut optimizer = ConformerOptimizer::distance_geometry(
                &mmat,
                &mut positions,
                &[],
                DistanceGeometryForceFieldParams {
                    weight_chiral: 0.2,
                    weight_fourth_dimension: 1.0,
                    basin_size_tolerance: 5.0,
                    extra_weights: None,
                    fixed_pair_points: None,
                    fixed_points: &[0],
                },
            )
            .expect("source optimizer");
            if optimizer.energy(None).expect("energy") > 0.00001 {
                while optimizer.minimize(200, 0.001, 1.0e-6).expect("minimize") != 0 {}
            }
        }
        assert_eq!(positions[0], vec![7.0, 8.0, 9.0, 0.0]);
    }
    #[test]
    fn source_fourth_dimension_exact_satisfied_two_point_bounds_via_checked_optimizer() {
        let mut mmat = BoundsMatrix::new(2);
        mmat.set_lower(0, 1, 1.0).expect("set lower");
        mmat.set_upper(0, 1, 1.0).expect("set upper");
        let mut positions = vec![vec![0.0, 0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0, 0.0]];
        {
            let mut optimizer = ConformerOptimizer::distance_geometry(
                &mmat,
                &mut positions,
                &[],
                DistanceGeometryForceFieldParams {
                    weight_chiral: 0.2,
                    weight_fourth_dimension: 1.0,
                    basin_size_tolerance: 5.0,
                    extra_weights: None,
                    fixed_pair_points: None,
                    fixed_points: &[],
                },
            )
            .expect("source optimizer");
            assert_eq!(optimizer.energy(None).expect("energy"), 0.0);
        }
        assert_eq!(
            positions,
            vec![vec![0.0, 0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0, 0.0]]
        );
    }
}

#[cfg(test)]
mod shared_chiral_rows_original_conditions {
    use super::*;
    fn fixture_flat_volume(a: usize, b: usize, c: usize, d: usize, pos: &[f64], dim: usize) -> f64 {
        let rows: Vec<Vec<f64>> = pos.chunks_exact(dim).map(<[f64]>::to_vec).collect();
        calc_chiral_volume_rows([a, b, c, d], &rows).expect("original checked flat fixture")
    }
    #[test]
    fn chiral_volume_flat_returns_signed_triple_product() {
        let pos = [
            1.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 0.0, 1.0, //
            0.0, 0.0, 0.0,
        ];

        assert_eq!(fixture_flat_volume(0, 1, 2, 3, &pos, 3), 1.0);
        assert_eq!(fixture_flat_volume(0, 2, 1, 3, &pos, 3), -1.0);
    }

    #[test]
    fn chiral_volume_flat_uses_idx4_as_reference_point() {
        let pos = [
            2.0, 2.0, 2.0, //
            3.0, 2.0, 2.0, //
            2.0, 3.0, 2.0, //
            2.0, 2.0, 3.0,
        ];

        assert_eq!(fixture_flat_volume(1, 2, 3, 0, &pos, 3), 1.0);
    }

    #[test]
    fn chiral_volume_flat_ignores_dimensions_after_first_three() {
        let pos = [
            1.0, 0.0, 0.0, 100.0, //
            0.0, 1.0, 0.0, 200.0, //
            0.0, 0.0, 1.0, 300.0, //
            0.0, 0.0, 0.0, 400.0,
        ];

        assert_eq!(fixture_flat_volume(0, 1, 2, 3, &pos, 4), 1.0);
    }

    #[test]
    fn chiral_volume_points_returns_signed_triple_product() {
        let pts = [
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
            vec![0.0, 0.0, 0.0],
        ];

        assert_eq!(
            calc_chiral_volume_rows([0, 1, 2, 3], &pts).expect("original checked four points"),
            1.0
        );
        assert_eq!(
            calc_chiral_volume_rows([0, 2, 1, 3], &pts).expect("original checked four points"),
            -1.0
        );
    }
}
