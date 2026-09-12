//! Core-owned numerical periodic-row data consumed by foundational chemistry.

use std::sync::OnceLock;

use cosmolkit_types::{Element, ElementInfo};
use thiserror::Error;

/// Error returned by a checked numerical periodic-table lookup.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Error)]
#[non_exhaustive]
pub enum PeriodicTableError {
    /// The element exists, but the requested isotope is absent from the
    /// pinned source table.
    #[error("unknown isotope {isotope} for element {element}")]
    UnknownIsotope { element: Element, isotope: u16 },
}

#[derive(Debug)]
struct IsotopeRow {
    isotope: u16,
    mass: f64,
    abundance: f64,
}

#[derive(Debug)]
struct PeriodicRow {
    period: u8,
    covalent_radius: f64,
    rb0: f64,
    van_der_waals_radius: f64,
    atomic_weight: f64,
    outer_electrons: i32,
    most_common_isotope: u16,
    most_common_isotope_mass: f64,
    valences: Box<[i32]>,
    isotopes: Box<[IsotopeRow]>,
}

static PERIODIC_ROWS: OnceLock<Box<[PeriodicRow]>> = OnceLock::new();

fn rows() -> &'static [PeriodicRow] {
    PERIODIC_ROWS.get_or_init(|| {
        // BEGIN RDKIT CPP FUNCTION PeriodicTable::PeriodicTable
        // RDKit✔️✔️: PeriodicTable::PeriodicTable() {
        // RDKit✔️✔️:   // it is assumed that the atomic atomData string constains atoms
        // RDKit✔️✔️:   // in sequence and no atoms are missing in between
        // RDKit✔️✔️:   byanum.clear();
        // RDKit✔️✔️:   byname.clear();
        // RDKit✔️✔️:   boost::char_separator<char> eolSep("\n");
        // RDKit✔️✔️:   tokenizer tokens(periodicTableAtomData, eolSep);
        // RDKit✔️✔️:   for (tokenizer::iterator token = tokens.begin(); token != tokens.end();
        // RDKit✔️✔️:        ++token) {
        // RDKit✔️✔️:     if (*token != " ") {
        // RDKit✔️✔️:       atomicData adata(*token);
        // RDKit✔️✔️:       std::string enam = adata.Symbol();
        // RDKit✔️✔️:       byname[enam] = adata.AtomicNum();
        // RDKit✔️✔️:       // there are, for backwards compatibility reasons, some duplicate rows for
        // RDKit✔️✔️:       // atomic numbers in the atomic_data data structure. It's ok to have
        // RDKit✔️✔️:       // multiple symbols map to the same atomic number (above), but we need to
        // RDKit✔️✔️:       // be sure that we only store one entry per atomic number.
        // RDKit✔️✔️:       // Note that this only works because the first atom in the adata list is
        // RDKit✔️✔️:       // the dummy atom (atomic number 0). This was #2784
        // RDKit✔️✔️:       if (rdcast<size_t>(adata.AtomicNum()) == byanum.size()) {
        // RDKit✔️✔️:         byanum.push_back(adata);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️🔝:
        // RDKit✔️🔝:   unsigned int lidx = 0;
        // RDKit✔️🔝:   std::istringstream istr;
        // RDKit✔️🔝:   istr.imbue(std::locale("C"));
        // RDKit✔️🔝:   while (isotopesAtomData[lidx] != "" && isotopesAtomData[lidx] != "EOS") {
        // RDKit✔️🔝:     tokenizer lines(isotopesAtomData[lidx++], eolSep);
        // RDKit✔️🔝:     boost::char_separator<char> spaceSep(" \t");
        // RDKit✔️🔝:     for (tokenizer::iterator line = lines.begin(); line != lines.end();
        // RDKit✔️🔝:          ++line) {
        // RDKit✔️🔝:       if (*line != " ") {
        // RDKit✔️🔝:         tokenizer tokens(*line, spaceSep);
        // RDKit✔️🔝:         tokenizer::iterator token = tokens.begin();
        // RDKit✔️🔝:         int anum;
        // RDKit✔️🔝:         istr.clear();
        // RDKit✔️🔝:         istr.str(*token);
        // RDKit✔️🔝:         istr >> anum;
        // RDKit✔️🔝:         atomicData &adata = byanum[anum];
        // RDKit✔️🔝:         ++token;
        // RDKit✔️🔝:         if (token == tokens.end()) {
        // RDKit✔️🔝:           continue;
        // RDKit✔️🔝:         }
        // RDKit✔️🔝:         ++token;
        // RDKit✔️🔝:         if (token == tokens.end()) {
        // RDKit✔️🔝:           continue;
        // RDKit✔️🔝:         }
        // RDKit✔️🔝:         unsigned int isotope;
        // RDKit✔️🔝:         istr.clear();
        // RDKit✔️🔝:         istr.str(*token);
        // RDKit✔️🔝:         istr >> isotope;
        // RDKit✔️🔝:         ++token;
        // RDKit✔️🔝:         if (token == tokens.end()) {
        // RDKit✔️🔝:           continue;
        // RDKit✔️🔝:         }
        // RDKit✔️🔝:         double mass;
        // RDKit✔️🔝:         istr.clear();
        // RDKit✔️🔝:         istr.str(*token);
        // RDKit✔️🔝:         istr >> mass;
        // RDKit✔️🔝:         ++token;
        // RDKit✔️🔝:         if (token == tokens.end()) {
        // RDKit✔️🔝:           continue;
        // RDKit✔️🔝:         }
        // RDKit✔️🔝:         double abundance;
        // RDKit✔️🔝:         istr.clear();
        // RDKit✔️🔝:         istr.str(*token);
        // RDKit✔️🔝:         istr >> abundance;
        // RDKit✔️🔝:         adata.d_isotopeInfoMap[isotope] = std::make_pair(mass, abundance);
        // RDKit✔️🔝:       }
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:   }
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION PeriodicTable::PeriodicTable
        // The normalized data removes the source symbol column because the
        // source loader itself skips it. Appending already-sorted isotope rows
        // to contiguous slices is linear and allocation-coalesced, improving
        // on one tree-node allocation and O(log n) insertion per source row;
        // immutable binary-search lookup remains O(log n) with identical keys.
        // BEGIN RDKIT CPP FUNCTION atomicData::atomicData
        // RDKit✔️✔️: atomicData::atomicData(const std::string &dataLine) {
        // RDKit✔️✔️:   boost::char_separator<char> spaceSep(" \t");
        // RDKit✔️✔️:   tokenizer tokens(dataLine, spaceSep);
        // RDKit✔️✔️:   tokenizer::iterator token = tokens.begin();
        // RDKit✔️✔️:   std::istringstream istr;
        // RDKit✔️✔️:   istr.imbue(std::locale("C"));
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // atomic number first
        // RDKit✔️✔️:   istr.clear();
        // RDKit✔️✔️:   istr.str(*token);
        // RDKit✔️✔️:   istr >> anum;
        // RDKit✔️✔️:   ++token;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // use atomic number to set element name
        // RDKit✔️✔️:   name = elementNames[anum];
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // element symbol
        // RDKit✔️✔️:   symb = *token;
        // RDKit✔️✔️:   ++token;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // covalent radius
        // RDKit✔️✔️:   istr.clear();
        // RDKit✔️✔️:   istr.str(*token);
        // RDKit✔️✔️:   istr >> row;
        // RDKit✔️✔️:   ++token;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // covalent radius
        // RDKit✔️✔️:   istr.clear();
        // RDKit✔️✔️:   istr.str(*token);
        // RDKit✔️✔️:   istr >> rCov;
        // RDKit✔️✔️:   ++token;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // rB0
        // RDKit✔️✔️:   istr.clear();
        // RDKit✔️✔️:   istr.str(*token);
        // RDKit✔️✔️:   istr >> rB0;
        // RDKit✔️✔️:   ++token;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   //  Van derWaal radius
        // RDKit✔️✔️:   istr.clear();
        // RDKit✔️✔️:   istr.str(*token);
        // RDKit✔️✔️:   istr >> rVdw;
        // RDKit✔️✔️:   ++token;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // atomic mass
        // RDKit✔️✔️:   istr.clear();
        // RDKit✔️✔️:   istr.str(*token);
        // RDKit✔️✔️:   istr >> mass;
        // RDKit✔️✔️:   ++token;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // number of outshell electrons
        // RDKit✔️✔️:   istr.clear();
        // RDKit✔️✔️:   istr.str(*token);
        // RDKit✔️✔️:   istr >> nVal;
        // RDKit✔️✔️:   ++token;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // most common isotope
        // RDKit✔️✔️:   istr.clear();
        // RDKit✔️✔️:   istr.str(*token);
        // RDKit✔️✔️:   istr >> commonIsotope;
        // RDKit✔️✔️:   ++token;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // most common isotopic mass
        // RDKit✔️✔️:   istr.clear();
        // RDKit✔️✔️:   istr.str(*token);
        // RDKit✔️✔️:   istr >> commonIsotopeMass;
        // RDKit✔️✔️:   ++token;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // now the valences
        // RDKit✔️✔️:   valence.clear();
        // RDKit✔️✔️:   while (token != tokens.end()) {
        // RDKit✔️✔️:     istr.clear();
        // RDKit✔️✔️:     istr.str(*token);
        // RDKit✔️✔️:     int tval;
        // RDKit✔️✔️:     istr >> tval;
        // RDKit✔️✔️:     valence.push_back(tval);
        // RDKit✔️✔️:     ++token;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION atomicData::atomicData
        let mut parsed = Vec::with_capacity(119);
        for line in PERIODIC_ROW_DATA.lines().filter(|line| !line.is_empty()) {
            let mut fields = line.split_whitespace();
            let atomic_number: usize = fields.next().expect("atomic number").parse().expect("u8");
            let source_symbol = fields.next().expect("symbol");
            let period = fields.next().expect("period").parse().expect("u8");
            let covalent_radius = fields
                .next()
                .expect("covalent radius")
                .parse()
                .expect("f64");
            let rb0 = fields.next().expect("Rb0").parse().expect("f64");
            let van_der_waals_radius = fields
                .next()
                .expect("van der Waals radius")
                .parse()
                .expect("f64");
            let atomic_weight = fields.next().expect("atomic weight").parse().expect("f64");
            let outer_electrons = fields
                .next()
                .expect("outer electrons")
                .parse()
                .expect("i32");
            let most_common_isotope = fields.next().expect("common isotope").parse().expect("u16");
            let most_common_isotope_mass = fields
                .next()
                .expect("common isotope mass")
                .parse()
                .expect("f64");
            let valences: Vec<i32> = fields
                .map(|field| field.parse().expect("valence"))
                .collect();
            assert!(!valences.is_empty(), "periodic row must contain a valence");
            assert_eq!(
                atomic_number,
                parsed.len(),
                "periodic rows must be contiguous"
            );
            assert_eq!(
                Element::from_atomic_number(atomic_number as u8)
                    .expect("checked source atomic number")
                    .symbol(),
                source_symbol,
                "canonical row symbol must match the vocabulary"
            );
            parsed.push(PeriodicRow {
                period,
                covalent_radius,
                rb0,
                van_der_waals_radius,
                atomic_weight,
                outer_electrons,
                most_common_isotope,
                most_common_isotope_mass,
                valences: valences.into_boxed_slice(),
                isotopes: Box::new([]),
            });
        }
        assert_eq!(
            parsed.len(),
            119,
            "periodic table must contain rows 0 through 118"
        );

        let mut isotopes_by_element: Vec<Vec<IsotopeRow>> =
            (0..parsed.len()).map(|_| Vec::new()).collect();
        let mut isotope_count = 0_usize;
        for line in ISOTOPE_ROW_DATA.lines().filter(|line| !line.is_empty()) {
            let mut fields = line.split_whitespace();
            let atomic_number: usize = fields
                .next()
                .expect("isotope atomic number")
                .parse()
                .expect("usize");
            let isotope = fields
                .next()
                .expect("isotope mass number")
                .parse()
                .expect("u16");
            let mass = fields.next().expect("isotope mass").parse().expect("f64");
            let abundance = fields
                .next()
                .expect("isotope abundance")
                .parse()
                .expect("f64");
            assert!(fields.next().is_none(), "unexpected isotope field");
            let element_rows = isotopes_by_element
                .get_mut(atomic_number)
                .expect("isotope atomic number must have a periodic row");
            assert!(
                element_rows
                    .last()
                    .is_none_or(|previous| previous.isotope < isotope),
                "source isotope rows must be unique and ordered within each element"
            );
            element_rows.push(IsotopeRow {
                isotope,
                mass,
                abundance,
            });
            isotope_count += 1;
        }
        assert_eq!(isotope_count, 3_111, "complete source isotope corpus");
        for (row, isotopes) in parsed.iter_mut().zip(isotopes_by_element) {
            row.isotopes = isotopes.into_boxed_slice();
        }
        parsed.into_boxed_slice()
    })
}

pub(crate) fn valences(atomic_number: u8) -> Option<&'static [i32]> {
    rows()
        .get(usize::from(atomic_number))
        .map(|row| row.valences.as_ref())
}

pub(crate) fn period(atomic_number: u8) -> Option<u8> {
    rows().get(usize::from(atomic_number)).map(|row| row.period)
}

pub(crate) fn outer_electrons(atomic_number: u8) -> Option<i32> {
    rows()
        .get(usize::from(atomic_number))
        .map(|row| row.outer_electrons)
}

pub(crate) fn rb0(atomic_number: u8) -> f64 {
    rows()
        .get(usize::from(atomic_number))
        .map_or(0.0, |row| row.rb0)
}

pub(crate) fn covalent_radius(atomic_number: u8) -> Option<f64> {
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::getRcovalent
    // RDKit✔️✔️: double getRcovalent(UINT atomicNumber) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   return byanum[atomicNumber].Rcov();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION PeriodicTable::getRcovalent
    rows()
        .get(usize::from(atomic_number))
        .map(|row| row.covalent_radius)
}

pub(crate) fn van_der_waals_radius(atomic_number: u8) -> Option<f64> {
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::getRvdw
    // RDKit✔️✔️: double getRvdw(UINT atomicNumber) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   return byanum[atomicNumber].Rvdw();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION PeriodicTable::getRvdw
    rows()
        .get(usize::from(atomic_number))
        .map(|row| row.van_der_waals_radius)
}

pub(crate) fn atomic_number_from_symbol(symbol: &str) -> Option<u8> {
    Element::from_symbol(symbol).map(Element::atomic_number)
}

pub(crate) fn symbol(atomic_number: u8) -> Option<&'static str> {
    Element::from_atomic_number(atomic_number).map(Element::symbol)
}

/// Return the complete dependency-light public metadata record for an element.
#[must_use]
pub fn element_info(element: Element) -> ElementInfo {
    let row = &rows()[usize::from(element.atomic_number())];
    ElementInfo {
        element,
        symbol: element.symbol(),
        atomic_number: element.atomic_number(),
        period: row.period,
        outer_electrons: row.outer_electrons,
        valences: row.valences.as_ref(),
        rb0: row.rb0,
        atomic_weight: row.atomic_weight,
    }
}

fn isotope_row(element: Element, isotope: u16) -> Option<&'static IsotopeRow> {
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::getMassForIsotope /
    // PeriodicTable::getAbundanceForIsotope
    // RDKit✔️✔️: const std::map<unsigned int, std::pair<double, double>> &m =
    // RDKit✔️✔️:     byanum[atomicNumber].d_isotopeInfoMap;
    // RDKit✔️✔️: std::map<unsigned int, std::pair<double, double>>::const_iterator item =
    // RDKit✔️✔️:     m.find(isotope);
    // RDKit✔️✔️: if (item == m.end()) {
    // RDKit✔️✔️:   return 0.0;
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:   return item->second.first;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION PeriodicTable::getMassForIsotope /
    // PeriodicTable::getAbundanceForIsotope
    let isotope_rows = &rows()[usize::from(element.atomic_number())].isotopes;
    isotope_rows
        .binary_search_by_key(&isotope, |row| row.isotope)
        .ok()
        .map(|index| &isotope_rows[index])
}

/// Return an element's atomic weight or one exact isotope mass.
pub fn atomic_mass(element: Element, isotope: Option<u16>) -> Result<f64, PeriodicTableError> {
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::getAtomicWeight
    // RDKit✔️✔️: double getAtomicWeight(UINT atomicNumber) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   double mass = byanum[atomicNumber].Mass();
    // RDKit✔️✔️:   return mass;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION PeriodicTable::getAtomicWeight
    match isotope {
        Some(isotope) => isotope_mass(element, isotope),
        None => Ok(rows()[usize::from(element.atomic_number())].atomic_weight),
    }
}

/// Return the mass number of the source-defined most common isotope.
#[must_use]
pub fn most_common_isotope(element: Element) -> u16 {
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::getMostCommonIsotope
    // RDKit✔️✔️: int getMostCommonIsotope(UINT atomicNumber) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   return byanum[atomicNumber].MostCommonIsotope();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION PeriodicTable::getMostCommonIsotope
    rows()[usize::from(element.atomic_number())].most_common_isotope
}

/// Return the source-defined mass of an element's most common isotope.
#[must_use]
pub fn most_common_isotope_mass(element: Element) -> f64 {
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::getMostCommonIsotopeMass
    // RDKit✔️✔️: double getMostCommonIsotopeMass(UINT atomicNumber) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   return byanum[atomicNumber].MostCommonIsotopeMass();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION PeriodicTable::getMostCommonIsotopeMass
    rows()[usize::from(element.atomic_number())].most_common_isotope_mass
}

/// Return an exact isotope mass, rejecting a missing source row.
pub fn isotope_mass(element: Element, isotope: u16) -> Result<f64, PeriodicTableError> {
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::getMassForIsotope
    // RDKit❗✔️: double getMassForIsotope(UINT atomicNumber, UINT isotope) const {
    // RDKit❗✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit❗✔️:   const std::map<unsigned int, std::pair<double, double>> &m =
    // RDKit❗✔️:       byanum[atomicNumber].d_isotopeInfoMap;
    // RDKit❗✔️:   std::map<unsigned int, std::pair<double, double>>::const_iterator item =
    // RDKit❗✔️:       m.find(isotope);
    // RDKit❗✔️:   if (item == m.end()) {
    // RDKit❗✔️:     return 0.0;
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     return item->second.first;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION PeriodicTable::getMassForIsotope
    // COSMolKit intentionally replaces the ambiguous missing-row 0.0 sentinel
    // with a structured error; present source rows retain their exact mass.
    isotope_row(element, isotope)
        .map(|row| row.mass)
        .ok_or(PeriodicTableError::UnknownIsotope { element, isotope })
}

/// Return an exact isotope abundance, preserving present zero-abundance rows.
pub fn isotope_abundance(element: Element, isotope: u16) -> Result<f64, PeriodicTableError> {
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::getAbundanceForIsotope
    // RDKit❗✔️: double getAbundanceForIsotope(UINT atomicNumber, UINT isotope) const {
    // RDKit❗✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit❗✔️:   const std::map<unsigned int, std::pair<double, double>> &m =
    // RDKit❗✔️:       byanum[atomicNumber].d_isotopeInfoMap;
    // RDKit❗✔️:   std::map<unsigned int, std::pair<double, double>>::const_iterator item =
    // RDKit❗✔️:       m.find(isotope);
    // RDKit❗✔️:   if (item == m.end()) {
    // RDKit❗✔️:     return 0.0;
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     return item->second.second;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION PeriodicTable::getAbundanceForIsotope
    // COSMolKit intentionally distinguishes a missing source row from a
    // present isotope whose source abundance is exactly zero.
    isotope_row(element, isotope)
        .map(|row| row.abundance)
        .ok_or(PeriodicTableError::UnknownIsotope { element, isotope })
}

pub(crate) fn more_electronegative(first: u8, second: u8) -> Option<bool> {
    let first_row = rows().get(usize::from(first))?;
    let second_row = rows().get(usize::from(second))?;
    Some(
        first_row.outer_electrons > second_row.outer_electrons
            || (first_row.outer_electrons == second_row.outer_electrons && first < second),
    )
}

#[cfg(test)]
fn isotope_count() -> usize {
    rows().iter().map(|row| row.isotopes.len()).sum()
}

// Verbatim canonical rows from RDKit's generated `periodicTableAtomData`.
// Uut/Uup duplicate-symbol rows are handled by the vocabulary alias lookup;
// RDKit also retains only the first numerical row for each atomic number.
const PERIODIC_ROW_DATA: &str = r#"0 * 0 0 0 0 0 0 0 0 -1
1 H 1 0.31 0.33 1.2 1.008 1 1 1.007825032 1
2 He 1 0.28 0.7 1.4 4.003 2 4 4.002603254 0
3 Li 2 1.28 1.23 2.2 6.941 1 7 7.01600455 1 -1
4 Be 2 0.96 0.9 1.9 9.012 2 9 9.0121822 2
5 B 2 0.84 0.82 1.8 10.812 3 11 11.0093054 3
6 C 2 0.76 0.77 1.7 12.011 4 12 12 4
7 N 2 0.71 0.7 1.6 14.007 5 14 14.003074 3
8 O 2 0.66 0.66 1.55 15.999 6 16 15.99491462 2
9 F 2 0.57 0.611 1.5 18.998 7 19 18.99840322 1
10 Ne 2 0.58 0.7 1.54 20.18 8 20 19.99244018 0
11 Na 3 1.66 1.54 2.4 22.99 1 23 22.98976928 1 -1
12 Mg 3 1.41 1.36 2.2 24.305 2 24 23.9850417 2 -1
13 Al 3 1.21 1.18 2.1 26.982 3 27 26.98153863 3
14 Si 3 1.11 0.937 2.1 28.086 4 28 27.97692653 4
15 P 3 1.07 0.89 1.95 30.974 5 31 30.97376163 3 5
16 S 3 1.05 1.04 1.8 32.067 6 32 31.972071 2 4 6
17 Cl 3 1.02 0.997 1.8 35.453 7 35 34.96885268 1
18 Ar 3 1.06 1.74 1.88 39.948 8 40 39.96238312 0
19 K 4 2.03 2.03 2.8 39.098 1 39 38.96370668 1 -1
20 Ca 4 1.76 1.74 2.4 40.078 2 40 39.96259098 2 -1
21 Sc 4 1.70 1.44 2.3 44.956 3 45 44.9559119 -1
22 Ti 4 1.60 1.32 2.15 47.867 4 48 47.9479463 -1
23 V 4 1.52 1.22 2.05 50.944 5 51 50.9439595 -1
24 Cr 4 1.39 1.18 2.05 51.996 6 52 51.9405075 -1
25 Mn 4 1.39 1.17 2.05 54.938 7 55 54.9380451 -1
26 Fe 4 1.32 1.17 2.05 55.845 8 56 55.9349375 -1
27 Co 4 1.26 1.16 2.0 58.933 9 59 58.933195 -1
28 Ni 4 1.24 1.15 2.0 58.693 10 58 57.9353429 -1
29 Cu 4 1.32 1.17 2.0 63.546 11 63 62.9295975 -1
30 Zn 4 1.22 1.25 2.1 65.39 2 64 63.9291422 -1
31 Ga 4 1.22 1.26 2.1 69.723 3 69 68.9255736 3
32 Ge 4 1.20 1.188 2.1 72.61 4 74 73.9211778 4
33 As 4 1.19 1.2 2.05 74.922 5 75 74.9215965 3 5
34 Se 4 1.20 1.17 1.9 78.96 6 80 79.9165213 2 4 6
35 Br 4 1.20 1.167 1.9 79.904 7 79 78.9183371 1
36 Kr 4 1.16 1.91 2.02 83.8 8 84 83.911507 0
37 Rb 5 2.20 2.16 2.9 85.468 1 85 84.91178974 1 -1
38 Sr 5 1.95 1.91 2.55 87.62 2 88 87.9056121 2 -1
39 Y 5 1.90 1.62 2.4 88.906 3 89 88.9058483 -1
40 Zr 5 1.75 1.45 2.3 91.224 4 90 89.9047044 -1
41 Nb 5 1.64 1.34 2.15 92.906 5 93 92.9063781 -1
42 Mo 5 1.54 1.3 2.1 95.94 6 98 97.9054082 -1
43 Tc 5 1.47 1.27 2.05 98 7 97 96.906365 -1
44 Ru 5 1.46 1.25 2.05 101.07 8 102 101.9043493 -1
45 Rh 5 1.42 1.25 2.0 102.906 9 103 102.905504 -1
46 Pd 5 1.39 1.28 2.05 106.42 10 106 105.903486 -1
47 Ag 5 1.45 1.34 2.1 107.868 11 107 106.905097 -1
48 Cd 5 1.44 1.48 2.2 112.412 2 114 113.9033585 -1
49 In 5 1.42 1.44 2.2 114.818 3 115 114.903878 3
50 Sn 5 1.39 1.385 2.25 118.711 4 120 119.9021947 2 4
51 Sb 5 1.39 1.4 2.2 121.76 5 121 120.9038157 3 5
52 Te 5 1.38 1.378 2.1 127.6 6 130 129.9062244 2 4 6
53 I 5 1.39 1.387 2.1 126.904 7 127 126.904473 1 3 5
54 Xe 5 1.40 1.98 2.16 131.29 8 132 131.9041535 0 2 4 6
55 Cs 6 2.44 2.35 3.0 132.905 1 133 132.9054519 1
56 Ba 6 2.15 1.98 2.7 137.328 2 138 137.9052472 2 -1
57 La 6 2.07 1.69 2.5 138.906 3 139 138.9063533 -1
58 Ce 6 2.04 1.83 2.48 140.116 4 140 139.9054387 -1
59 Pr 6 2.03 1.82 2.47 140.908 3 141 140.9076528 -1
60 Nd 6 2.01 1.81 2.45 144.24 4 142 141.9077233 -1
61 Pm 6 1.99 1.8 2.43 145 5 145 144.912749 -1
62 Sm 6 1.98 1.8 2.42 150.36 6 152 151.9197324 -1
63 Eu 6 1.98 1.99 2.4 151.964 7 153 152.9212303 -1
64 Gd 6 1.96 1.79 2.38 157.25 8 158 157.9241039 -1
65 Tb 6 1.94 1.76 2.37 158.925 9 159 158.9253468 -1
66 Dy 6 1.92 1.75 2.35 162.5 10 164 163.9291748 -1
67 Ho 6 1.92 1.74 2.33 164.93 11 165 164.9303221 -1
68 Er 6 1.89 1.73 2.32 167.26 12 166 165.9302931 -1
69 Tm 6 1.90 1.72 2.3 168.934 13 169 168.9342133 -1
70 Yb 6 1.87 1.94 2.28 173.04 14 174 173.9388621 -1
71 Lu 6 1.87 1.72 2.27 174.967 15 175 174.9407718 -1
72 Hf 6 1.75 1.44 2.25 178.49 4 180 179.94655 -1
73 Ta 6 1.70 1.34 2.2 180.948 5 181 180.9479958 -1
74 W 6 1.62 1.3 2.1 183.84 6 184 183.9509312 -1
75 Re 6 1.51 1.28 2.05 186.207 7 187 186.9557531 -1
76 Os 6 1.44 1.26 2.0 190.23 8 192 191.9614807 -1
77 Ir 6 1.41 1.27 2.0 192.217 9 193 192.9629264 -1
78 Pt 6 1.36 1.3 2.05 195.078 10 195 194.9647911 -1
79 Au 6 1.36 1.34 2.1 196.967 11 197 196.9665687 -1
80 Hg 6 1.32 1.49 2.05 200.59 2 202 201.970643 -1
81 Tl 6 1.45 1.48 2.2 204.383 3 205 204.9744275 -1
82 Pb 6 1.46 1.48 2.3 207.2 4 208 207.9766521 2 4
83 Bi 6 1.48 1.45 2.3 208.98 5 209 208.9803987 3 5
84 Po 6 1.40 1.46 2.0 209 6 209 208.9824304 2 4 6
85 At 6 1.50 1.45 2.0 210 7 210 209.987148 1 3 5
86 Rn 6 1.50 2.4 2.0 222 8 222 222.0175706 0
87 Fr 7 2.6 2 2.0 223 1 223 223.0197359 1
88 Ra 7 2.2 1.9 2.0 226 2 226 226.0254026 2 -1
89 Ac 7 2.15 1.88 2.0 227 3 227 227.0277521 -1
90 Th 7 2.06 1.79 2.4 232.038 4 232 232.0380553 -1
91 Pa 7 2.00 1.61 2.0 231.036 3 231 231.035884 -1
92 U 7 1.96 1.58 2.3 238.029 4 238 238.0507882 -1
93 Np 7 1.90 1.55 2.0 237 5 236 236.04657 -1
94 Pu 7 1.87 1.53 2.0 244 6 238 238.0495599 -1
95 Am 7 1.80 1.07 2.0 243 7 241 241.0568291 -1
96 Cm 7 1.69 0 2.0 247 8 243 243.0613891 -1
97 Bk 7 1.9 0 2.0 247 9 247 247.070307 -1
98 Cf 7 1.9 0 2.0 251 10 249 249.0748535 -1
99 Es 7 1.9 0 2.0 252 11 252 252.08298 -1
100 Fm 7 1.9 0 2.0 257 12 257 257.095105 -1
101 Md 7 1.9 0 2.0 258 13 258 258.098431 -1
102 No 7 1.9 0 2.0 259 14 259 259.10103 -1
103 Lr 7 1.9 0 2.0 262 15 262 262.10963 -1
104 Rf 7 1.9 0 2.0 267 2 267 267.12153 -1
105 Db 7 1.9 0 2.0 268 2 268 268.12545 -1
106 Sg 7 1.9 0 2.0 269 2 271 271.13347 -1
107 Bh 7 1.9 0 2.0 270 2 270 270.13362 -1
108 Hs 7 1.9 0 2.0 269 2 269 269.13406 -1
109 Mt 7 1.9 0 2.0 278 2 278 278.15481 -1
110 Ds 7 1.9 0 2.0 281 2 281 281.16206 -1
111 Rg 7 1.9 0 2.0 281 2 281 281.16537 -1
112 Cn 7 1.9 0 2.0 285 2 285 285.17411 -1
113 Nh 7 1.36 0 2.0 284 2 284 284.17873 -1
114 Fl 7 1.43 0 2.0 289 2 289 289.19042 -1
115 Mc 7 1.62 0 2.0 288 2 288 288.19274 -1
116 Lv 7 1.75 0 2.0 293 2 293 293.20449 -1
117 Ts 7 1.65 0 2.0 292 2 292 292.20746 -1
118 Og 7 1.57 0 2.0 294 2 294 294.21392 -1"#;

// Generated mechanically from the complete `isotopesAtomData` array in the
// pinned RDKit `GraphMol/atomic_data.cpp`. The source loader ignores the
// symbol column, so the normalized rows retain exactly the four consumed
// fields: atomic number, isotope number, mass, and abundance.
const ISOTOPE_ROW_DATA: &str = include_str!("periodic_table_isotopes.tsv");

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_rows_are_contiguous_complete_and_self_consistent() {
        for atomic_number in 0_u8..=118 {
            let row = rows()
                .get(usize::from(atomic_number))
                .expect("canonical row");
            let element = Element::from_atomic_number(atomic_number).expect("canonical element");
            assert_eq!(
                atomic_number_from_symbol(element.symbol()),
                Some(atomic_number)
            );
            assert_eq!(symbol(atomic_number), Some(element.symbol()));
            assert!(!row.valences.is_empty());
            assert!(row.period <= 7);
            assert!(row.rb0 >= 0.0);
            assert!(row.atomic_weight >= 0.0);
        }
    }

    #[test]
    fn source_rows_preserve_valence_period_outer_electron_radius_and_weight_fields() {
        assert_eq!(valences(0), Some(&[-1][..]));
        assert_eq!(valences(6), Some(&[4][..]));
        assert_eq!(valences(16), Some(&[2, 4, 6][..]));
        assert_eq!(valences(54), Some(&[0, 2, 4, 6][..]));
        assert_eq!(period(6), Some(2));
        assert_eq!(period(118), Some(7));
        assert_eq!(outer_electrons(6), Some(4));
        assert_eq!(outer_electrons(71), Some(15));
        assert_eq!(rb0(6), 0.77);
        assert_eq!(covalent_radius(6), Some(0.76));
        assert_eq!(van_der_waals_radius(6), Some(1.7));
        assert_eq!(element_info(Element::C).atomic_weight, 12.011);
    }

    #[test]
    fn out_of_range_rows_fail_closed() {
        assert_eq!(valences(119), None);
        assert_eq!(period(119), None);
        assert_eq!(outer_electrons(119), None);
        assert_eq!(symbol(119), None);
        assert_eq!(rb0(119), 0.0);
        assert_eq!(covalent_radius(119), None);
        assert_eq!(van_der_waals_radius(119), None);
    }

    #[test]
    fn legacy_symbols_share_canonical_numeric_rows() {
        assert_eq!(atomic_number_from_symbol("Uut"), Some(113));
        assert_eq!(atomic_number_from_symbol("Uup"), Some(115));
        assert_eq!(symbol(113), Some("Nh"));
        assert_eq!(symbol(115), Some("Mc"));
    }

    #[test]
    fn table_construction_is_shared_after_first_lookup() {
        let first = rows().as_ptr();
        let second = rows().as_ptr();
        assert_eq!(first, second);
    }

    #[test]
    fn complete_isotope_corpus_is_loaded_in_source_order() {
        assert_eq!(isotope_count(), 3_111);
        for row in rows() {
            assert!(
                row.isotopes
                    .windows(2)
                    .all(|pair| pair[0].isotope < pair[1].isotope)
            );
        }
    }

    #[test]
    fn internal_radius_and_electronegativity_branches_match_source_rows() {
        assert_eq!(covalent_radius(1), Some(0.31));
        assert_eq!(van_der_waals_radius(1), Some(1.2));
        assert_eq!(covalent_radius(118), Some(1.57));
        assert_eq!(van_der_waals_radius(118), Some(2.0));

        assert_eq!(more_electronegative(8, 6), Some(true));
        assert_eq!(more_electronegative(6, 8), Some(false));
        assert_eq!(more_electronegative(6, 14), Some(true));
        assert_eq!(more_electronegative(14, 6), Some(false));
        assert_eq!(more_electronegative(6, 6), Some(false));
        assert_eq!(more_electronegative(119, 6), None);
        assert_eq!(more_electronegative(6, 119), None);
    }
}
