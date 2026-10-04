//! Private fixed native literals, shared by the public and storage regressions.
//! Provenance/qualification/license: prepared_property_presence.md/.LICENSE.

use cosmolkit_model::{MoleculeProperties, PropertyValue, TopologyBlock};
use std::collections::{BTreeMap, BTreeSet};

pub const CASES: [(usize, &str); 3] = [
    (
        2508,
        "CC(C)(c1ccc(OCC[C@H](O)Cl)cc1)c1ccc(OCC[C@@H](O)Cl)cc1",
    ),
    (
        227,
        "CCC/C=C/C=C/C(O)=N[C@@H](Cc1c[nH]c2cc(Cl)ccc12)C(O)=N[C@H](CC(=N)O)C(O)=N[C@@H](CC(=O)O)C(O)=N[C@@H]1C(O)=NCC(O)=N[C@@H](CCCN)C(O)=N[C@@H](CC(=O)O)C(O)=N[C@H](C)C(O)=N[C@@H](CC(=O)O)C(O)=NCC(O)=N[C@H](C)C(O)=N[C@@H]([C@H](C)CC(=O)O)C(O)=N[C@@H](CC(=O)c2ccc(Cl)cc2N)C(=O)O[C@@H]1C",
    ),
    (
        155,
        "CCC[C@H]1C[C@@H]2CC[C@@H](O2)[C@H](C)C(=O)O[C@@H]([C@H](CC)[C@H]2CC[C@@H](C[C@@H](CCC)N(C)C)O2)[C@H](C)[C@@H]2CC[C@@H](O2)[C@@H](CC)C(=O)O1",
    ),
];

fn unhex(text: &str) -> String {
    String::from_utf8(
        text.as_bytes()
            .chunks_exact(2)
            .map(|pair| u8::from_str_radix(std::str::from_utf8(pair).unwrap(), 16).unwrap())
            .collect(),
    )
    .unwrap()
}

fn block<'a>(fixed: &'a str, key: &str) -> &'a str {
    let marker = format!("BLOCK\t{key}\n");
    fixed
        .split_once(&marker)
        .unwrap()
        .1
        .split_once("END\n")
        .unwrap()
        .0
}

pub fn check(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    xy: &[[f64; 2]],
    line: usize,
    orientation: bool,
) -> Vec<String> {
    let fixed = include_str!(
        "../../../../testdata/depict_2d/expected/rdkit/prepared_property_presence.tsv"
    );
    let label = format!("line:{line}:P2:O{}", u8::from(orientation));
    let cell = fixed
        .lines()
        .find(|row| row.starts_with(&format!("CELL\t{label}\t")))
        .unwrap();
    let keys = cell.split('\t').collect::<Vec<_>>();
    let expected_topology = block(fixed, keys[2]);
    let expected_properties = block(fixed, keys[3]);
    let expected_xy = block(fixed, keys[6]);
    let mut errors = Vec::new();
    let mut actual_rows = Vec::new();
    for atom in &topology.atoms {
        actual_rows.push(format!(
            "ATOM\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            atom.id().index(),
            atom.element().atomic_number(),
            atom.isotope().unwrap_or(0),
            atom.formal_charge(),
            atom.explicit_hydrogens(),
            u8::from(atom.no_implicit()),
            atom.radical_electrons(),
            u8::from(atom.is_aromatic()),
            atom.hybridization().rdkit_code(),
            atom.chiral_tag().rdkit_code(),
            atom.atom_map().unwrap_or(0)
        ));
    }
    for bond in &topology.bonds {
        let refs = bond
            .stereo_atoms()
            .map(|ids| format!("{},{}", ids[0].index(), ids[1].index()))
            .unwrap_or_default();
        actual_rows.push(format!(
            "BOND\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{refs}",
            bond.id().index(),
            bond.begin().index(),
            bond.end().index(),
            bond.order().rdkit_code(),
            u8::from(bond.is_aromatic()),
            u8::from(bond.is_conjugated()),
            bond.direction().rdkit_code(),
            bond.stereo().rdkit_code()
        ));
    }
    let native_rows = expected_topology
        .lines()
        .filter(|row| row.starts_with("ATOM\t") || row.starts_with("BOND\t"))
        .collect::<Vec<_>>();
    if actual_rows.iter().map(String::as_str).collect::<Vec<_>>() != native_rows {
        errors.push(format!("{label}: complete ordered topology differs\nactual={actual_rows:?}\nexpected={native_rows:?}"));
    }
    for (owner, count) in [
        ("ATOM", topology.atoms.len()),
        ("BOND", topology.bonds.len()),
    ] {
        for index in 0..count {
            let mut expected = BTreeMap::new();
            let mut computed = BTreeSet::new();
            for row in expected_topology
                .lines()
                .filter(|row| row.starts_with("PROP\t"))
            {
                let fields = row.split('\t').collect::<Vec<_>>();
                if fields[1] != owner || fields[2].parse::<usize>().unwrap() != index {
                    continue;
                }
                let key = unhex(fields[4]);
                if key == "__computedProps" {
                    assert_eq!(fields[5], "12");
                    continue;
                }
                let text = unhex(fields[7]);
                let value = match fields[5] {
                    "1" | "6" => PropertyValue::Int(text.parse().unwrap()),
                    "3" => PropertyValue::String(text),
                    other => panic!("unmodeled fixed native tag {other}"),
                };
                if fields[6] == "1" {
                    computed.insert(key.clone());
                }
                expected.insert(key, value);
            }
            let (actual, markers) = if owner == "ATOM" {
                (
                    topology.atoms[index].props(),
                    topology.atoms[index].computed_prop_names(),
                )
            } else {
                (
                    topology.bonds[index].props(),
                    topology.bonds[index].computed_prop_names(),
                )
            };
            if actual != &expected || markers != &computed {
                errors.push(format!("{label}: {owner}/{index} complete scalar properties/computed differ actual={actual:?}/{markers:?} expected={expected:?}/{computed:?}"));
            }
        }
    }
    let mut expected_mol = BTreeMap::new();
    let mut computed_mol = BTreeSet::new();
    for row in expected_properties.lines() {
        let fields = row.split('\t').collect::<Vec<_>>();
        let key = unhex(fields[4]);
        if key == "__computedProps" {
            assert_eq!(fields[5], "12");
            continue;
        }
        assert!(matches!(fields[5], "1" | "3"));
        if fields[6] == "1" {
            computed_mol.insert(key.clone());
        }
        expected_mol.insert(key, unhex(fields[7]));
    }
    if properties.props() != &expected_mol || properties.computed_prop_names() != &computed_mol {
        errors.push(format!("{label}: modeled MOL scalar properties/computed differ actual={properties:?} expected={expected_mol:?}/{computed_mol:?}"));
    }
    let rows = expected_xy.lines().collect::<Vec<_>>();
    if xy.len() != rows.len() {
        errors.push(format!("{label}: XY lengths {}/{}", xy.len(), rows.len()));
    }
    for (index, row) in rows.iter().enumerate() {
        let fields = row.split('\t').collect::<Vec<_>>();
        assert_eq!(fields[1].parse::<usize>().unwrap(), index);
        if let Some(actual) = xy.get(index) {
            for axis in 0..2 {
                let expected = f64::from_bits(fields[axis + 2].parse().unwrap());
                if (actual[axis] - expected).abs() > 1e-8 {
                    errors.push(format!(
                        "{label}: XY {index}/{axis} {} expected {expected}",
                        actual[axis]
                    ));
                }
            }
        }
    }
    errors
}
