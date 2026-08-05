use cosmolkit_core::{
    ChiralTag, Mol2ReadParams, Molecule, RdkitPdbMolProfile, SdfReadParams,
    mol_from_mol2_block_like_rdkit,
};

#[derive(Debug, Clone, PartialEq, Eq)]
struct AtomChiralState {
    tag: ChiralTag,
    permutation: Option<u32>,
    non_explicit_3d: Option<String>,
}

fn atom_chiral_state(molecule: &Molecule) -> Vec<AtomChiralState> {
    molecule
        .atoms()
        .iter()
        .map(|atom| AtomChiralState {
            tag: atom.chiral_tag(),
            permutation: atom.chiral_permutation(),
            non_explicit_3d: atom.prop("_NonExplicit3DChirality").map(str::to_owned),
        })
        .collect()
}

fn assert_tetrahedral_parser_result(molecule: &Molecule, expected_coordinates: &[[f64; 3]]) {
    let mut expected_state = vec![
        AtomChiralState {
            tag: ChiralTag::TetrahedralCcw,
            permutation: None,
            non_explicit_3d: Some("1".to_string()),
        },
        AtomChiralState {
            tag: ChiralTag::Unspecified,
            permutation: None,
            non_explicit_3d: None,
        },
        AtomChiralState {
            tag: ChiralTag::Unspecified,
            permutation: None,
            non_explicit_3d: None,
        },
        AtomChiralState {
            tag: ChiralTag::Unspecified,
            permutation: None,
            non_explicit_3d: None,
        },
    ];
    if molecule.num_atoms() == 5 {
        expected_state.push(AtomChiralState {
            tag: ChiralTag::Unspecified,
            permutation: None,
            non_explicit_3d: None,
        });
    }

    assert_eq!(atom_chiral_state(molecule), expected_state);
    let conformer = molecule
        .conformers_3d()
        .first()
        .expect("parser path must retain its 3D conformer");
    assert!(conformer.is_3d());
    assert_eq!(conformer.coordinates(), expected_coordinates);
}

#[test]
fn sdf_3d_path_uses_completed_chiral_assignment_kernel() {
    let mol_block = "\
parser_path_sdf
  COSMolKit

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    1.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -1.0000   -1.0000   -1.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
";
    let molecule = Molecule::from_mol_block_with_params(
        mol_block,
        SdfReadParams {
            sanitize: false,
            remove_hs: false,
            ..SdfReadParams::default()
        },
    )
    .expect("3D V2000 parser path must succeed");

    assert_tetrahedral_parser_result(
        &molecule,
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, -1.0, -1.0],
        ],
    );
}

#[test]
fn mol2_3d_path_assigns_before_sanitize_or_hydrogen_removal() {
    let mol2 = "\
@<TRIPOS>MOLECULE
parser_path_mol2
5 4
SMALL
NO_CHARGES
@<TRIPOS>ATOM
1 C1 0 0 0 C.3
2 F1 1 0 0 F
3 Cl1 0 1 0 Cl
4 Br1 0 0 1 Br
5 H1 -1 -1 -1 H
@<TRIPOS>BOND
1 1 2 1
2 1 3 1
3 1 4 1
4 1 5 1
";
    let record = mol_from_mol2_block_like_rdkit(
        mol2,
        Mol2ReadParams {
            sanitize: false,
            remove_hs: false,
            ..Mol2ReadParams::default()
        },
    )
    .expect("MOL2 parser path must succeed")
    .expect("MOL2 input must contain one molecule");

    assert_tetrahedral_parser_result(
        &record.molecule,
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, -1.0, -1.0],
        ],
    );
}

#[test]
fn pdb_molecule_3d_path_uses_completed_chiral_assignment_kernel() {
    let pdb = "\
HETATM    1  C1  UNK A   1       0.000   0.000   0.000  1.00  0.00           C
HETATM    2  F1  UNK A   1       1.000   0.000   0.000  1.00  0.00           F
HETATM    3 CL1  UNK A   1       0.000   1.000   0.000  1.00  0.00          CL
HETATM    4 BR1  UNK A   1       0.000   0.000   1.000  1.00  0.00          BR
HETATM    5  H1  UNK A   1      -1.000  -1.000  -1.000  1.00  0.00           H
CONECT    1    2    3    4    5
END
";
    let molecule = Molecule::from_pdb_block_with_options(
        pdb,
        RdkitPdbMolProfile {
            sanitize: false,
            remove_hs: false,
            flavor: 0,
            proximity_bonding: false,
        },
    )
    .expect("PDB molecule conversion path must succeed");

    assert_tetrahedral_parser_result(
        &molecule,
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, -1.0, -1.0],
        ],
    );
}

#[test]
fn cxsmiles_3d_path_uses_completed_chiral_assignment_kernel() {
    let molecule =
        Molecule::from_smiles_with_sanitize("C(F)(Cl)Br |(0,0,0;1,0,0;0,1,0;0,0,1)|", false)
            .expect("3D CXSMILES parser path must succeed");

    assert_tetrahedral_parser_result(
        &molecule,
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
    );
}
