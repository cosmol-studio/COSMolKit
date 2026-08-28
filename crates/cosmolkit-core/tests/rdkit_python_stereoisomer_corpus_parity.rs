use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use cosmolkit_core::{
    BondDirection, BondStereo, ChiralTag, ControllingAtom, Molecule, PotentialStereoOptions,
    StereoCenter, StereoDescriptor, StereoInfo, StereoSpecified, StereoType, StereoisomerOptions,
};
use num_bigint::BigInt;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

mod common;
use common::parity_data;

const OUTPUT_NAME: &str = "python_stereoisomer_corpus.jsonl";
const POTENTIAL_PROFILE_IDS: [&str; 2] = ["preserve_possible", "clean_possible"];
const ENUMERATION_PROFILE_IDS: [&str; 4] = [
    "default_bounded",
    "all_assigned_bounded",
    "non_unique_bounded",
    "seeded_three",
];

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
struct MoleculeState {
    canonical_smiles: String,
    atom_chiral_tags: Vec<String>,
    bonds: Vec<BondState>,
    conformer_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
struct BondState {
    direction: String,
    stereo: String,
    stereo_atoms: Vec<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
struct StereoInfoState {
    stereo_type: String,
    specified: String,
    center_kind: String,
    center_index: Option<usize>,
    descriptor: String,
    permutation: u32,
    controlling_atoms: Vec<Option<usize>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
struct PotentialProfile {
    id: String,
    clean_it: bool,
    flag_possible: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
struct PotentialRun {
    profile: PotentialProfile,
    records: Vec<StereoInfoState>,
    state: MoleculeState,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
struct EnumerationProfile {
    id: String,
    only_unassigned: bool,
    only_stereo_groups: bool,
    max_isomers: usize,
    rand: Option<i64>,
    unique: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
struct EnumerationRun {
    profile: EnumerationProfile,
    theoretical_count: String,
    bounded_out: bool,
    outputs: Vec<MoleculeState>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
struct OracleRecord {
    schema_version: u32,
    row: usize,
    smiles: String,
    parse_status: String,
    error_type: Option<String>,
    error_message: Option<String>,
    source_state: Option<MoleculeState>,
    potential_stereo: Vec<PotentialRun>,
    enumeration: Vec<EnumerationRun>,
}

fn read_records(profile: &str) -> Vec<OracleRecord> {
    let path = parity_data::expected_path_for_profile("stereo", "rdkit", profile, OUTPUT_NAME);
    let file = File::open(&path)
        .unwrap_or_else(|error| panic!("failed to open {}: {error}", path.display()));
    BufReader::new(file)
        .lines()
        .enumerate()
        .map(|(line_index, line)| {
            let line = line.unwrap_or_else(|error| {
                panic!(
                    "failed to read {} line {}: {error}",
                    path.display(),
                    line_index + 1
                )
            });
            serde_json::from_str(&line).unwrap_or_else(|error| {
                panic!(
                    "failed to parse {} line {}: {error}",
                    path.display(),
                    line_index + 1
                )
            })
        })
        .collect()
}

fn read_smiles(path: &Path) -> Vec<String> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()))
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(str::to_owned)
        .collect()
}

fn chiral_tag_name(tag: ChiralTag) -> &'static str {
    match tag {
        ChiralTag::Unspecified => "Unspecified",
        ChiralTag::TetrahedralCw => "TetrahedralCw",
        ChiralTag::TetrahedralCcw => "TetrahedralCcw",
        ChiralTag::Other => "Other",
        ChiralTag::Tetrahedral => "Tetrahedral",
        ChiralTag::Allene => "Allene",
        ChiralTag::SquarePlanar => "SquarePlanar",
        ChiralTag::TrigonalBipyramidal => "TrigonalBipyramidal",
        ChiralTag::Octahedral => "Octahedral",
    }
}

fn bond_direction_name(direction: BondDirection) -> &'static str {
    match direction {
        BondDirection::None => "None",
        BondDirection::BeginWedge => "BeginWedge",
        BondDirection::BeginDash => "BeginDash",
        BondDirection::EndDownRight => "EndDownRight",
        BondDirection::EndUpRight => "EndUpRight",
        BondDirection::EitherDouble => "EitherDouble",
        BondDirection::Unknown => "Unknown",
    }
}

fn bond_stereo_name(stereo: BondStereo) -> &'static str {
    match stereo {
        BondStereo::None => "None",
        BondStereo::Any => "Any",
        BondStereo::Z => "Z",
        BondStereo::E => "E",
        BondStereo::Cis => "Cis",
        BondStereo::Trans => "Trans",
        BondStereo::AtropCw => "AtropCw",
        BondStereo::AtropCcw => "AtropCcw",
    }
}

fn enumeration_direction_gauge(molecule: &Molecule) -> Vec<Option<&'static str>> {
    let mut parent = (0..molecule.num_bonds()).collect::<Vec<_>>();
    let mut directed = vec![false; molecule.num_bonds()];
    let mut participating = vec![false; molecule.num_bonds()];
    let mut incident = vec![Vec::new(); molecule.num_atoms()];
    for bond in molecule.bonds() {
        if matches!(
            bond.direction(),
            BondDirection::EndDownRight | BondDirection::EndUpRight
        ) {
            directed[bond.id().index()] = true;
            incident[bond.begin().index()].push(bond.id().index());
            incident[bond.end().index()].push(bond.id().index());
        }
    }

    fn find(parent: &mut [usize], mut index: usize) -> usize {
        while parent[index] != index {
            parent[index] = parent[parent[index]];
            index = parent[index];
        }
        index
    }
    fn union(parent: &mut [usize], left: usize, right: usize) {
        let left_root = find(parent, left);
        let right_root = find(parent, right);
        if left_root != right_root {
            parent[right_root] = left_root;
        }
    }

    for bond in molecule.bonds() {
        if !matches!(
            bond.stereo(),
            BondStereo::Cis | BondStereo::Trans | BondStereo::Z | BondStereo::E
        ) {
            continue;
        }
        let component = incident[bond.begin().index()]
            .iter()
            .chain(&incident[bond.end().index()])
            .copied()
            .collect::<Vec<_>>();
        if let Some(first) = component.first().copied() {
            for index in &component {
                participating[*index] = true;
            }
            for other in component.iter().copied().skip(1) {
                union(&mut parent, first, other);
            }
        }
    }

    let mut minimum = vec![usize::MAX; molecule.num_bonds()];
    for index in 0..molecule.num_bonds() {
        if participating[index] {
            let root = find(&mut parent, index);
            minimum[root] = minimum[root].min(index);
        }
    }
    let mut flip = vec![false; molecule.num_bonds()];
    for (root, first) in minimum.iter().copied().enumerate() {
        if first != usize::MAX && molecule.bonds()[first].direction() == BondDirection::EndDownRight
        {
            flip[root] = true;
        }
    }

    molecule
        .bonds()
        .iter()
        .map(|bond| {
            if !participating[bond.id().index()] {
                return None;
            }
            let root = find(&mut parent, bond.id().index());
            Some(match (bond.direction(), flip[root]) {
                (BondDirection::EndDownRight, false) | (BondDirection::EndUpRight, true) => {
                    "EndDownRight"
                }
                (BondDirection::EndUpRight, false) | (BondDirection::EndDownRight, true) => {
                    "EndUpRight"
                }
                _ => unreachable!("direction gauge contains only slash directions"),
            })
        })
        .collect()
}

fn molecule_state(
    molecule: &Molecule,
    normalize_enumeration_directions: bool,
) -> Result<MoleculeState, String> {
    let direction_gauge =
        normalize_enumeration_directions.then(|| enumeration_direction_gauge(molecule));
    Ok(MoleculeState {
        canonical_smiles: molecule
            .to_smiles(true)
            .map_err(|error| error.to_string())?,
        atom_chiral_tags: molecule
            .atoms()
            .iter()
            .map(|atom| chiral_tag_name(atom.chiral_tag()).to_owned())
            .collect(),
        bonds: molecule
            .bonds()
            .iter()
            .map(|bond| BondState {
                direction: direction_gauge
                    .as_ref()
                    .and_then(|directions| directions[bond.id().index()])
                    .unwrap_or_else(|| bond_direction_name(bond.direction()))
                    .to_owned(),
                stereo: bond_stereo_name(bond.stereo()).to_owned(),
                stereo_atoms: bond.stereo_atoms().map_or_else(Vec::new, |atoms| {
                    atoms.into_iter().map(|atom| atom.index()).collect()
                }),
            })
            .collect(),
        conformer_count: molecule.conformers_2d().len() + molecule.conformers_3d().len(),
    })
}

fn stereo_type_name(stereo_type: StereoType) -> &'static str {
    match stereo_type {
        StereoType::Unspecified => "unspecified",
        StereoType::AtomTetrahedral => "atom_tetrahedral",
        StereoType::AtomSquarePlanar => "atom_square_planar",
        StereoType::AtomTrigonalBipyramidal => "atom_trigonal_bipyramidal",
        StereoType::AtomOctahedral => "atom_octahedral",
        StereoType::BondDouble => "bond_double",
        StereoType::BondEvenCumulene => "bond_even_cumulene",
        StereoType::BondAtropisomer => "bond_atropisomer",
    }
}

fn specified_name(specified: StereoSpecified) -> &'static str {
    match specified {
        StereoSpecified::Unspecified => "unspecified",
        StereoSpecified::Specified => "specified",
        StereoSpecified::Unknown => "unknown",
    }
}

fn descriptor_name(descriptor: StereoDescriptor) -> &'static str {
    match descriptor {
        StereoDescriptor::None => "none",
        StereoDescriptor::TetrahedralClockwise => "tetrahedral_clockwise",
        StereoDescriptor::TetrahedralCounterclockwise => "tetrahedral_counterclockwise",
        StereoDescriptor::BondCis => "bond_cis",
        StereoDescriptor::BondTrans => "bond_trans",
        StereoDescriptor::BondAtropClockwise => "bond_atrop_clockwise",
        StereoDescriptor::BondAtropCounterclockwise => "bond_atrop_counterclockwise",
    }
}

fn stereo_info_state(info: &StereoInfo) -> StereoInfoState {
    let stereo_type = info.stereo_type();
    let (center_kind, center_index) = match info.center() {
        StereoCenter::Atom(atom) => ("atom", Some(atom.index())),
        StereoCenter::Bond(bond) => ("bond", Some(bond.index())),
        StereoCenter::Missing => (
            if stereo_type.is_atom_centered() {
                "atom"
            } else {
                "bond"
            },
            None,
        ),
    };
    StereoInfoState {
        stereo_type: stereo_type_name(stereo_type).to_owned(),
        specified: specified_name(info.specified()).to_owned(),
        center_kind: center_kind.to_owned(),
        center_index,
        descriptor: descriptor_name(info.descriptor()).to_owned(),
        permutation: info.permutation(),
        controlling_atoms: info
            .controlling_atoms()
            .iter()
            .map(|atom| match atom {
                ControllingAtom::Missing => None,
                ControllingAtom::Atom(atom) => Some(atom.index()),
            })
            .collect(),
    }
}

fn potential_profiles() -> Vec<PotentialProfile> {
    vec![
        PotentialProfile {
            id: "preserve_possible".to_owned(),
            clean_it: false,
            flag_possible: true,
        },
        PotentialProfile {
            id: "clean_possible".to_owned(),
            clean_it: true,
            flag_possible: true,
        },
    ]
}

fn enumeration_profiles() -> Vec<EnumerationProfile> {
    vec![
        EnumerationProfile {
            id: "default_bounded".to_owned(),
            only_unassigned: true,
            only_stereo_groups: false,
            max_isomers: 8,
            rand: None,
            unique: true,
        },
        EnumerationProfile {
            id: "all_assigned_bounded".to_owned(),
            only_unassigned: false,
            only_stereo_groups: false,
            max_isomers: 8,
            rand: None,
            unique: true,
        },
        EnumerationProfile {
            id: "non_unique_bounded".to_owned(),
            only_unassigned: true,
            only_stereo_groups: false,
            max_isomers: 8,
            rand: None,
            unique: false,
        },
        EnumerationProfile {
            id: "seeded_three".to_owned(),
            only_unassigned: true,
            only_stereo_groups: false,
            max_isomers: 3,
            rand: Some(61_453),
            unique: true,
        },
    ]
}

fn actual_record(row: usize, smiles: &str) -> Result<OracleRecord, String> {
    let source = match Molecule::from_smiles(smiles) {
        Ok(source) => source,
        Err(_) => {
            return Ok(OracleRecord {
                schema_version: 2,
                row,
                smiles: smiles.to_owned(),
                parse_status: "none".to_owned(),
                error_type: None,
                error_message: None,
                source_state: None,
                potential_stereo: Vec::new(),
                enumeration: Vec::new(),
            });
        }
    };
    let source_state = molecule_state(&source, false)?;

    let mut potential_stereo = Vec::new();
    for profile in potential_profiles() {
        let analysis = source
            .analyze_potential_stereo(PotentialStereoOptions {
                clean: profile.clean_it,
                flag_possible: profile.flag_possible,
            })
            .map_err(|error| format!("potential profile {}: {error}", profile.id))?;
        potential_stereo.push(PotentialRun {
            profile,
            records: analysis.stereo_info.iter().map(stereo_info_state).collect(),
            state: molecule_state(&analysis.molecule, false)?,
        });
    }

    let mut enumeration = Vec::new();
    for profile in enumeration_profiles() {
        let options = StereoisomerOptions {
            try_embedding: false,
            only_unassigned: profile.only_unassigned,
            only_stereo_groups: profile.only_stereo_groups,
            max_isomers: profile.max_isomers,
            random_seed: profile.rand.map(BigInt::from),
            unique: profile.unique,
        };
        let theoretical_count = source
            .stereoisomer_count(&options)
            .map_err(|error| format!("count profile {}: {error}", profile.id))?;
        let bounded_out = profile.id != "seeded_three" && theoretical_count > 8_u8.into();
        let outputs = if bounded_out {
            Vec::new()
        } else {
            source
                .stereoisomers(options)
                .map_err(|error| format!("iterator profile {}: {error}", profile.id))?
                .map(|result| {
                    result
                        .map_err(|error| format!("output profile {}: {error}", profile.id))
                        .and_then(|molecule| molecule_state(&molecule, true))
                })
                .collect::<Result<Vec<_>, _>>()?
        };
        enumeration.push(EnumerationRun {
            profile,
            theoretical_count: theoretical_count.to_string(),
            bounded_out,
            outputs,
        });
    }

    if molecule_state(&source, false)? != source_state {
        return Err("potential-stereo or enumeration mutated the source molecule".to_owned());
    }
    Ok(OracleRecord {
        schema_version: 2,
        row,
        smiles: smiles.to_owned(),
        parse_status: "ok".to_owned(),
        error_type: None,
        error_message: None,
        source_state: Some(source_state),
        potential_stereo,
        enumeration,
    })
}

fn validate_record_identity(records: &[OracleRecord], corpus: &[String]) -> Result<(), String> {
    if records.len() != corpus.len() {
        return Err(format!(
            "record count {} does not match corpus count {}",
            records.len(),
            corpus.len()
        ));
    }
    for (row, (record, smiles)) in records.iter().zip(corpus).enumerate() {
        if record.schema_version != 2 {
            return Err(format!(
                "row {row} has schema version {}",
                record.schema_version
            ));
        }
        if record.row != row {
            return Err(format!("row identity {} does not match {row}", record.row));
        }
        if record.smiles != *smiles {
            return Err(format!("row {row} input identity differs"));
        }
        match record.parse_status.as_str() {
            "ok" => {
                let potential_ids = record
                    .potential_stereo
                    .iter()
                    .map(|run| run.profile.id.as_str())
                    .collect::<Vec<_>>();
                if potential_ids != POTENTIAL_PROFILE_IDS {
                    return Err(format!("row {row} potential profile order differs"));
                }
                let enumeration_ids = record
                    .enumeration
                    .iter()
                    .map(|run| run.profile.id.as_str())
                    .collect::<Vec<_>>();
                if enumeration_ids != ENUMERATION_PROFILE_IDS {
                    return Err(format!("row {row} enumeration profile order differs"));
                }
            }
            "none" | "error" => {
                if record.source_state.is_some()
                    || !record.potential_stereo.is_empty()
                    || !record.enumeration.is_empty()
                {
                    return Err(format!("row {row} parse failure contains chemistry output"));
                }
            }
            other => return Err(format!("row {row} has unknown parse status {other:?}")),
        }
    }
    Ok(())
}

fn leaf_count(value: &serde_json::Value) -> usize {
    match value {
        serde_json::Value::Array(values) => values.iter().map(leaf_count).sum(),
        serde_json::Value::Object(values) => values.values().map(leaf_count).sum(),
        serde_json::Value::Null
        | serde_json::Value::Bool(_)
        | serde_json::Value::Number(_)
        | serde_json::Value::String(_) => 1,
    }
}

fn enumeration_difference(expected: &[EnumerationRun], actual: &[EnumerationRun]) -> String {
    if expected.len() != actual.len() {
        return format!(
            "enumeration run count differs: expected {}, got {}",
            expected.len(),
            actual.len()
        );
    }
    for (expected_run, actual_run) in expected.iter().zip(actual) {
        if expected_run.profile != actual_run.profile {
            return format!(
                "enumeration profile differs: expected {:?}, got {:?}",
                expected_run.profile, actual_run.profile
            );
        }
        let profile = &expected_run.profile.id;
        if expected_run.theoretical_count != actual_run.theoretical_count {
            return format!(
                "profile {profile} count differs: expected {}, got {}",
                expected_run.theoretical_count, actual_run.theoretical_count
            );
        }
        if expected_run.bounded_out != actual_run.bounded_out {
            return format!(
                "profile {profile} bounded flag differs: expected {}, got {}",
                expected_run.bounded_out, actual_run.bounded_out
            );
        }
        if expected_run.outputs.len() != actual_run.outputs.len() {
            return format!(
                "profile {profile} output count differs: expected {}, got {}",
                expected_run.outputs.len(),
                actual_run.outputs.len()
            );
        }
        for (output_index, (expected_output, actual_output)) in expected_run
            .outputs
            .iter()
            .zip(&actual_run.outputs)
            .enumerate()
        {
            if expected_output.canonical_smiles != actual_output.canonical_smiles {
                return format!(
                    "profile {profile} output {output_index} canonical SMILES differs: expected {:?}, got {:?}",
                    expected_output.canonical_smiles, actual_output.canonical_smiles
                );
            }
            if expected_output.atom_chiral_tags != actual_output.atom_chiral_tags {
                return format!(
                    "profile {profile} output {output_index} atom chiral tags differ: expected {:?}, got {:?}",
                    expected_output.atom_chiral_tags, actual_output.atom_chiral_tags
                );
            }
            for (bond_index, (expected_bond, actual_bond)) in expected_output
                .bonds
                .iter()
                .zip(&actual_output.bonds)
                .enumerate()
            {
                if expected_bond != actual_bond {
                    return format!(
                        "profile {profile} output {output_index} bond {bond_index} differs: expected {expected_bond:?}, got {actual_bond:?}"
                    );
                }
            }
            if expected_output.bonds.len() != actual_output.bonds.len() {
                return format!(
                    "profile {profile} output {output_index} bond count differs: expected {}, got {}",
                    expected_output.bonds.len(),
                    actual_output.bonds.len()
                );
            }
            if expected_output.conformer_count != actual_output.conformer_count {
                return format!(
                    "profile {profile} output {output_index} conformer count differs: expected {}, got {}",
                    expected_output.conformer_count, actual_output.conformer_count
                );
            }
        }
    }
    "enumeration metadata differs".to_owned()
}

fn compare_record(expected: &OracleRecord) -> Result<usize, String> {
    let actual = actual_record(expected.row, &expected.smiles)?;
    if actual != *expected {
        if actual.parse_status != expected.parse_status {
            return Err(format!(
                "parse status differs: expected {}, got {}",
                expected.parse_status, actual.parse_status
            ));
        }
        if actual.source_state != expected.source_state {
            return Err("source state differs".to_owned());
        }
        if actual.potential_stereo != expected.potential_stereo {
            return Err("potential-stereo records or cleaned state differ".to_owned());
        }
        if actual.enumeration != expected.enumeration {
            return Err(enumeration_difference(
                &expected.enumeration,
                &actual.enumeration,
            ));
        }
        return Err("record metadata differs".to_owned());
    }
    let value = serde_json::to_value(expected).map_err(|error| error.to_string())?;
    Ok(leaf_count(&value))
}

fn assert_profile(
    profile: &str,
    corpus_path: &Path,
    expected_records: usize,
    expected_comparisons: usize,
) {
    let records = read_records(profile);
    let corpus = read_smiles(corpus_path);
    assert_eq!(records.len(), expected_records, "{profile}: record count");
    validate_record_identity(&records, &corpus)
        .unwrap_or_else(|error| panic!("{profile}: {error}"));
    let results = records
        .par_iter()
        .map(|record| {
            compare_record(record)
                .map_err(|error| format!("row {} ({}): {error}", record.row, record.smiles))
        })
        .collect::<Vec<_>>();
    let failures = results
        .iter()
        .filter_map(|result| result.as_ref().err())
        .take(32)
        .cloned()
        .collect::<Vec<_>>();
    assert!(
        failures.is_empty(),
        "{profile} parity failures: {failures:?}"
    );
    let comparisons = results.into_iter().map(Result::unwrap).sum::<usize>();
    assert_eq!(
        comparisons, expected_comparisons,
        "{profile}: exact comparison total changed"
    );
}

#[test]
fn manifests_profiles_and_rejection_guards_are_enforced() {
    let small_path = parity_data::repo_root().join("testdata/smiles/corpus/smiles_small.smi");
    let corpus = read_smiles(&small_path);
    let records = read_records("python_stereoisomer_small");
    validate_record_identity(&records, &corpus).expect("committed small profile identity");

    assert!(validate_record_identity(&records[..records.len() - 1], &corpus).is_err());
    let filtered = records
        .iter()
        .filter(|record| record.parse_status == "ok")
        .cloned()
        .collect::<Vec<_>>();
    assert!(validate_record_identity(&filtered, &corpus).is_err());

    let mut mismatching = records.clone();
    mismatching[0].smiles.push('C');
    assert!(validate_record_identity(&mismatching, &corpus).is_err());

    let mut missing_profile = records.clone();
    let successful = missing_profile
        .iter_mut()
        .find(|record| record.parse_status == "ok")
        .expect("small corpus success row");
    successful.enumeration.pop();
    assert!(validate_record_identity(&missing_profile, &corpus).is_err());

    let mut chemistry_mismatch = records
        .iter()
        .find(|record| record.parse_status == "ok")
        .expect("small corpus success row")
        .clone();
    chemistry_mismatch
        .source_state
        .as_mut()
        .expect("source state")
        .canonical_smiles
        .push('C');
    assert!(compare_record(&chemistry_mismatch).is_err());
}

#[test]
fn small_corpus_matches_pinned_rdkit_exactly() {
    assert_profile(
        "python_stereoisomer_small",
        &parity_data::repo_root().join("testdata/smiles/corpus/smiles_small.smi"),
        152,
        89_047,
    );
}

#[test]
fn corpus_5000_matches_pinned_rdkit_exactly() {
    assert_profile(
        "python_stereoisomer_5000",
        &parity_data::repo_root().join("testdata/smiles/corpus/smiles_5000.smi"),
        5_000,
        4_680_490,
    );
}
