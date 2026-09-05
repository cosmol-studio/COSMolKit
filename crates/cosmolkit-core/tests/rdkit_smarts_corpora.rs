use std::io::{BufRead, BufReader};

use cosmolkit_core::{
    BondDirection, ChiralTag, Molecule, SmartsParseParams, SmilesWriteParams, SubstructMatchParams, get_atom_smarts,
    get_bond_smarts, mol_fragment_to_cx_smarts, mol_fragment_to_smarts, mol_from_smarts, mol_to_cx_smarts,
    mol_to_smarts, try_get_substruct_matches_with_params,
};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct GoldenRecord {
    case_id: String,
    kind: String,
    source: String,
    smarts: String,
    target_smiles: Option<String>,
    parse_ok: bool,
    num_atoms: Option<usize>,
    num_bonds: Option<usize>,
    written_smarts: Option<String>,
    atom_mappings: Option<Vec<Vec<usize>>>,
}

fn records() -> impl Iterator<Item = GoldenRecord> {
    let path = cosmolkit_test_support::expected_path_for_profile("smarts", "rdkit", "smarts_source", "smarts.jsonl");
    let file = std::fs::File::open(&path).unwrap_or_else(|error| panic!("failed to open {}: {error}", path.display()));
    BufReader::new(file).lines().enumerate().map(move |(index, line)| {
        let line = line.unwrap_or_else(|error| panic!("failed to read {}:{}: {error}", path.display(), index + 1));
        serde_json::from_str(&line)
            .unwrap_or_else(|error| panic!("failed to parse {}:{}: {error}", path.display(), index + 1))
    })
}

#[test]
fn rdkit_smarts_parser_and_writer_source_corpus_matches_exactly() {
    for record in records() {
        let parsed = mol_from_smarts(&record.smarts, &SmartsParseParams::default());
        assert_eq!(
            parsed.is_ok(),
            record.parse_ok,
            "{} from {}: {}",
            record.case_id,
            record.source,
            record.smarts
        );
        let Ok(query) = parsed else {
            continue;
        };
        assert_eq!(query.num_atoms(), record.num_atoms.unwrap(), "{}", record.case_id);
        assert_eq!(query.num_bonds(), record.num_bonds.unwrap(), "{}", record.case_id);
        assert_eq!(
            mol_to_smarts(&query, &SmilesWriteParams::default()).unwrap(),
            record.written_smarts.unwrap(),
            "{}",
            record.case_id
        );
    }
}

#[test]
fn rdkit_smarts_substructure_source_corpus_matches_ordered_mappings() {
    for record in records().filter(|record| record.kind == "match") {
        let target = Molecule::from_smiles(record.target_smiles.as_deref().unwrap())
            .unwrap_or_else(|error| panic!("{} target: {error}", record.case_id));
        let query = mol_from_smarts(&record.smarts, &SmartsParseParams::default())
            .unwrap_or_else(|error| panic!("{} query: {error}", record.case_id));
        let actual = try_get_substruct_matches_with_params(&target, &query, &SubstructMatchParams::default())
            .unwrap_or_else(|error| panic!("{} match: {error}", record.case_id))
            .into_iter()
            .map(|matched| matched.atom_mapping)
            .collect::<Vec<_>>();
        assert_eq!(actual, record.atom_mappings.unwrap(), "{}", record.case_id);
    }
}

#[test]
fn smarts_fragment_and_cx_apis_compose_with_source_corpus() {
    for record in records().filter(|record| record.parse_ok) {
        let first = mol_from_smarts(&record.smarts, &SmartsParseParams::default()).unwrap();
        let text = mol_to_smarts(&first, &SmilesWriteParams::default()).unwrap();
        let second = mol_from_smarts(&text, &SmartsParseParams::default()).unwrap();
        assert_eq!(first.num_atoms(), second.num_atoms(), "{}", record.case_id);
        assert_eq!(first.num_bonds(), second.num_bonds(), "{}", record.case_id);
        assert_eq!(
            mol_to_smarts(&second, &SmilesWriteParams::default()).unwrap(),
            text,
            "{}",
            record.case_id
        );

        for atom in first.atoms() {
            let atom_text = get_atom_smarts(atom, &SmilesWriteParams::default()).unwrap();
            let atom_query = mol_from_smarts(&atom_text, &SmartsParseParams::default()).unwrap();
            assert_eq!(atom_query.num_atoms(), 1, "{}", record.case_id);
        }
        for bond in first.bonds() {
            let _ = get_bond_smarts(bond, &SmilesWriteParams::default(), None).unwrap();
        }
        let atoms = first.atoms().iter().map(|atom| atom.id()).collect::<Vec<_>>();
        let bonds = first.bonds().iter().map(|bond| bond.id()).collect::<Vec<_>>();
        let bonds_to_use = (!bonds.is_empty()).then_some(bonds.as_slice());
        let fragment = mol_fragment_to_smarts(&first, &SmilesWriteParams::default(), &atoms, bonds_to_use).unwrap();
        assert_eq!(
            mol_from_smarts(&fragment, &SmartsParseParams::default())
                .unwrap()
                .num_atoms(),
            first.num_atoms(),
            "{}",
            record.case_id
        );
        let cx = mol_to_cx_smarts(&first, &SmilesWriteParams::default()).unwrap();
        let fragment_cx =
            mol_fragment_to_cx_smarts(&first, &SmilesWriteParams::default(), &atoms, bonds_to_use).unwrap();
        mol_from_smarts(&cx, &SmartsParseParams::default())
            .unwrap_or_else(|error| panic!("{} CXSMARTS {cx:?}: {error}", record.case_id));
        mol_from_smarts(&fragment_cx, &SmartsParseParams::default())
            .unwrap_or_else(|error| panic!("{} fragment CXSMARTS {fragment_cx:?}: {error}", record.case_id));
    }
}

#[test]
fn smarts_directional_bond_state_matches_rdkit_grammar_paths() {
    for (smarts, bond_index, expected) in [
        ("C/C", 0, BondDirection::EndUpRight),
        (r"C\C", 0, BondDirection::EndDownRight),
        ("C/;@C", 0, BondDirection::EndUpRight),
        ("C1CC/1", 2, BondDirection::EndUpRight),
        (r"C1CC\1", 2, BondDirection::EndDownRight),
        ("C/1CC1", 2, BondDirection::EndUpRight),
        (r"C\1CC1", 2, BondDirection::EndDownRight),
    ] {
        let molecule =
            mol_from_smarts(smarts, &SmartsParseParams::default()).unwrap_or_else(|error| panic!("{smarts}: {error}"));
        assert_eq!(molecule.bonds()[bond_index].direction(), expected, "{smarts}");
    }
}

#[test]
fn smarts_atom_expression_keeps_the_first_chirality_specification() {
    for (smarts, written, tag, permutation) in [
        ("[C@,C@@]", "[C@,C]", ChiralTag::TetrahedralCcw, None),
        ("[C@@,C@]", "[C@@,C]", ChiralTag::TetrahedralCw, None),
        ("[C@SP1,C@SP2]", "[C,C]", ChiralTag::SquarePlanar, Some(1)),
    ] {
        let molecule =
            mol_from_smarts(smarts, &SmartsParseParams::default()).unwrap_or_else(|error| panic!("{smarts}: {error}"));
        let atom = &molecule.atoms()[0];
        assert_eq!(atom.chiral_tag(), tag, "{smarts}");
        assert_eq!(atom.chiral_permutation(), permutation, "{smarts}");
        assert_eq!(
            mol_to_smarts(&molecule, &SmilesWriteParams::default()).unwrap(),
            written,
            "{smarts}"
        );
    }
}
