# COSMolKit

`cosmolkit` is the public Rust API and runtime owner for COSMolKit, a Rust-native cheminformatics and structural biology toolkit. It is the sole supported Rust entrypoint for `Molecule`, operation contracts, and the feature APIs being migrated from the legacy implementation. Detached model values and source-backed algorithms are being moved into dedicated workspace crates behind this facade.

For a concise Rust-native cheminformatics overview, see <https://tools.cosmol.org/rust-cheminformatics>.

## Cargo features

The default `full` feature preserves the complete COSMolKit API. Applications
that need a smaller dependency surface can disable defaults and enable only the
capabilities they use:

```toml
cosmolkit = { version = "0.3.0", default-features = false, features = ["fingerprints"] }
```

| Feature | Capability | Implied features |
|---|---|---|
| `io` | MOL/SDF, MOL2, XYZ, PDB, and mmCIF IO | — |
| `inchi` | InChI and InChIKey conversion | — |
| `fingerprints` | Morgan, AtomPair, Pattern, Topological, MACCS, Layered, and Avalon fingerprints | `io` |
| `descriptors` | Molecular descriptors | — |
| `depict` | SVG and PNG molecule depiction | `io` |
| `batch` | Ordered parallel batch APIs | `fingerprints`, `depict`, `io` |
| `serialization` | Compact binary molecule serialization | — |
| `stereoisomers` | Stereoisomer enumeration | — |
| `hashing` | Molecular hashes and scaffolds | `fingerprints`, `io` |
| `confseq` | ConfSeq decoding | — |

Feature selection only controls compile-time API and dependency composition;
it does not change the behavior of an enabled operation. Python wheels always
build the complete `full` capability set.

## Documentation

- Rust API documentation: <https://docs.rs/cosmolkit/latest/cosmolkit/>
- Project validation: [`VALIDATION.md`](https://github.com/cosmol-studio/COSMolKit/blob/main/VALIDATION.md)
- Validation scope and evidence: [`VALIDATION.md`](https://github.com/cosmol-studio/COSMolKit/blob/main/VALIDATION.md) · [COSMolKit Web Tools](https://tools.cosmol.org/tools)
- Target crate architecture: [`dev/final_target_architecture.md`](https://github.com/cosmol-studio/COSMolKit/blob/main/dev/final_target_architecture.md)
- Project overview and Python package notes: [`README.md`](https://github.com/cosmol-studio/COSMolKit/blob/main/README.md)

## Validation Status

COSMolKit treats parity as **source-backed semantic equivalence within explicitly documented boundaries**, not as statistical agreement of final outputs. Compatibility-critical chemistry is implemented as a line-by-line, source-backed port with explicit operation contracts and traceable correspondence to pinned upstream code. Validation corpora verify that port; they are not used to iteratively tune heuristic reimplementations until outputs happen to agree.

The comparison boundary therefore extends well beyond final strings. Covered surfaces compare exact bytes, bits, return status, complete atom and bond state, stereochemistry, derived state and invariants, **RNG state, seed handling, and random draw sequences where stochastic behavior is part of the contract**, every matrix entry, coordinates, energies, and every gradient component where applicable. Discrete results must match exactly; declared numerical tolerances reach `1e-8` for matrix entries and `1e-6` for coordinates, energies, and gradients. **99% or 99.9% agreement remains unfinished when any covered mismatch exists.**

This boundary is stress-tested against a complete ChEMBL 37 profile: 2,897,819 source records, 2,897,804 of them mutually parseable, across 34 repository-defined sharded phases against pinned RDKit `2026.03.1`. The profile performs billions of comparisons, expands parameter spaces into matrices of up to 768 branches, repeats complete matrices to expose instability, permutes operation order, and checks scalar, one-thread, multi-thread, batch, and shared-object concurrent paths.

Every discovered mismatch is traced back to the corresponding upstream logic, corrected at the source-port level, and permanently retained as a focused regression rather than hidden by corpus-specific adjustments. This discipline limits **semantic debt** by preventing convenient local fixes from accumulating into undocumented chemistry behavior.

The parity suite uses three complementary validation layers. The complete ChEMBL 37 profile provides large-scale stress coverage; the maintained 5,000-record corpus runs exhaustive parameter matrices not yet practical across the full ChEMBL profile; and the 152-record project corpus keeps focused regressions fast enough for daily testing.

See [`VALIDATION.md`](https://github.com/cosmol-studio/COSMolKit/blob/main/VALIDATION.md) for exact corpus eligibility, comparison counts, tolerances, per-feature boundaries, focused regressions, and upstream surfaces outside the current claim.

## Installation

```toml
cargo add cosmolkit
```

## Quick Start

```rust
use cosmolkit::{Molecule, SmilesWriteParams};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mol = Molecule::from_smiles("CCO")?;
    let mol = mol.with_2d_coordinates()?;

    let smiles = mol.to_smiles_with_params(&SmilesWriteParams::default())?;
    let svg = mol.to_svg(300, 300)?;

    println!("{smiles}");
    println!("{}", svg.len());
    Ok(())
}
```

## Molecule Operations

Normal `Molecule` operations return new values and leave the receiver
unchanged:

```rust
let mol = Molecule::from_smiles("CCO")?;
let with_h = mol.with_hydrogens()?;
assert_ne!(mol.num_atoms(), with_h.num_atoms());
```

In-place operations are explicit and always end with `_`:

```rust
let mut mol = Molecule::from_smiles("CCO")?;
mol.add_hydrogens_()?;
mol.sanitize_()?;
```

The trailing underscore is reserved for in-place mutation on public `Molecule`
methods; it has no other meaning. In-place operations prioritize avoiding the
operation-system working-copy clone when molecule blocks are uniquely owned. If
an in-place operation returns an error, the receiver is not guaranteed to equal
its pre-call value and may retain partial changes, while its internal storage
remains complete. Use the non-mutating operation when failure-preserving value
semantics are required.

Stable molecule operations include assigning atom chiral tags from a selected
3D conformer through `with_chiral_tags_from_structure()` and its explicit
in-place counterpart `assign_chiral_tags_from_structure_()`.

## Protein Structures

```rust
use cosmolkit::Protein;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let protein = Protein::from_pdb("1crn.pdb")?;
    let summary = protein.selection_summary();

    println!("chains: {}", summary.chains);
    println!("residues: {}", summary.residues);
    println!("atoms: {}", summary.atoms);
    Ok(())
}
```

## Batch Workflows

```rust
use cosmolkit::{BatchErrorMode, MoleculeBatch};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let smiles = vec![
        "CCO".to_string(),
        "c1ccccc1".to_string(),
        "CC(=O)O".to_string(),
    ];

    let batch = MoleculeBatch::from_smiles_list(&smiles)
        .with_parallel_jobs(Some(8))
        .with_2d_coordinates(BatchErrorMode::Strict)?;

    let out = batch.to_smiles_list(BatchErrorMode::Strict)?;
    println!("{out:?}");
    Ok(())
}
```

## InChI

The Rust facade exposes the four audited scalar InChI APIs directly:

```rust
use cosmolkit::{Molecule, inchi_to_inchi_key, mol_from_inchi, mol_to_inchi};

let molecule = Molecule::from_smiles("C")?;
let generated = mol_to_inchi(&molecule, None)?;
assert_eq!(generated.inchi, b"InChI=1S/CH4/h1H4");

let key = inchi_to_inchi_key(&generated.inchi)?;
assert_eq!(key.key, b"VNWKTOKETHGBQD-UHFFFAOYSA-N");

let parsed = mol_from_inchi(&generated.inchi, false, false)?;
assert!(parsed.molecule.is_some());
```

Pinned official InChI v1.07.5 and RDKit 2026.03.1 establish exact parity for
source-defined behavior in this boundary. Official-C undefined behavior on the
audited `NormalizeAndCompare` initial-allocation path is mapped to a
deterministic structured allocation error. MolBlock, SDF/V3000, IXA, AuxInfo,
INCHIGEN, version-query, and extended-polymer InChI APIs are not exposed.

## Molecular Descriptors

The facade re-exports the source-backed descriptor functions from
`cosmolkit-core`:

```rust
use cosmolkit::{
    Molecule, calc_chi_0, calc_mol_formula, calc_mol_wt, calc_mqns,
    calc_num_aromatic_rings,
};

let molecule = Molecule::from_smiles("c1ccccc1O")?;
assert_eq!(calc_mol_formula(&molecule, false, true)?, "C6H6O");
assert!(calc_mol_wt(&molecule, false)? > 94.0);
assert_eq!(calc_num_aromatic_rings(&molecule)?, 1);
assert!(calc_chi_0(&molecule) > 0.0);
assert_eq!(calc_mqns(&molecule)?.len(), 42);
```

The documented descriptor surface includes molecular properties, connectivity
and shape indices, Lipinski and ring/stereo counts, MQN, Labute ASA, and
SlogP/SMR VSA. Supported rows and parameter combinations are checked
field-by-field against pinned RDKit golden data; unmodeled source states return
an explicit descriptor error.

## Fingerprints

The Rust facade exposes source-backed Morgan, AtomPair, Topological Torsion,
MACCS, RDKit topological, Avalon, and Layered fingerprints. ``TopologicalTorsion*`` is
the ordered atom-path torsion family; ``TopologicalFingerprint*`` remains
RDKit's distinct path/subgraph ``RDKFingerprintMol`` family. The applicable
families can also return typed provenance:

```rust
use cosmolkit::{
    AtomPairFingerprintParams, AvalonFingerprintParams, LayeredFingerprintLayers,
    LayeredFingerprintParams, Molecule,
    TopologicalFingerprintOutputRequest,
    TopologicalFingerprintParams, TopologicalTorsionFingerprintOutputRequest,
    TopologicalTorsionFingerprintParams, TopologicalTorsionFingerprintVector,
    topological_torsion_fingerprint, topological_torsion_fingerprint_with_output,
    topological_torsion_sparse_count_fingerprint,
};

let molecule = Molecule::from_smiles("c1ccccc1O")?;
let topological = molecule.topological_fingerprint(
    &TopologicalFingerprintParams::default(),
)?;
let provenance = molecule.topological_fingerprint_with_output(
    &TopologicalFingerprintParams::default(),
    TopologicalFingerprintOutputRequest {
        atom_bits: true,
        bit_info: true,
    },
)?;
let avalon = molecule.avalon_fingerprint(&AvalonFingerprintParams::default())?;
let atom_pair = molecule.atom_pair_fingerprint(&AtomPairFingerprintParams::default())?;
let layered = molecule.layered_fingerprint(&LayeredFingerprintParams {
    layers: LayeredFingerprintLayers::SUBSTRUCTURE,
    ..Default::default()
})?;
let torsion_params = TopologicalTorsionFingerprintParams::default();
let torsion_ids = topological_torsion_sparse_count_fingerprint(&molecule, &torsion_params)?;
let torsion_bits = topological_torsion_fingerprint(&molecule, &torsion_params)?;
let torsion_provenance = topological_torsion_fingerprint_with_output(
    &molecule,
    &torsion_params,
    TopologicalTorsionFingerprintOutputRequest {
        vector: TopologicalTorsionFingerprintVector::Bit,
        bit_paths: true,
        ..Default::default()
    },
)?;

assert_eq!(topological.n_bits(), 2048);
assert!(provenance.output.atom_bits.is_some());
assert_eq!(avalon.n_bits(), 512);
assert_eq!(atom_pair.n_bits(), 2048);
assert_eq!(layered.n_bits(), 2048);
assert!(!torsion_ids.nonzero_elements().is_empty());
assert_eq!(torsion_bits.n_bits(), 2048);
assert!(torsion_provenance.additional_output.is_some());
```

Topological Torsion also exposes sparse-bit and folded-count forms, ordered
``MoleculeBatch`` conveniences, shared ``AdditionalOutput`` provenance, and
three explicitly typed legacy adapters. Invalid arguments return
``FingerprintError``; batch calculation errors retain the original record
index. Exact parity is continuously checked against pinned RDKit 2026.03.1 on
focused branch fixtures and every row of a 5,000-molecule, nine-profile
matrix. The complete ChEMBL 37 audit additionally covers all 2,897,804
mutually parseable records through 127,503,376 exact vector and provenance
comparisons. Legacy adapters preserve their historical unfolded-size and
``n_bits_per_entry`` threshold differences while delegating to the same
chemistry and vector-assembly core.

The documented topological and Avalon profiles are checked against pinned
RDKit across all 2,897,804 mutually parseable ChEMBL 37 molecules. The
full-corpus audit completed 113,014,356 exact comparisons over 14 topological
vectors, 23 Avalon vectors, and two complete topological provenance outputs
with zero mismatches. The committed 5,000-row matrices remain the continuous
regression gates for these profiles.

AtomPair is additionally checked across all 2,897,804 mutually parseable
ChEMBL 37 molecules, covering 118,809,964 comparisons over 40 vectors and one
complete provenance output per molecule with zero mismatches.

Layered exposes the six source layers, arbitrary retained source flags,
inclusive path bounds, rooted linear or branched enumeration, exact-width bit
masks, and seeded atom counts through one read-only core while preserving the
upstream ``0.7.0`` compatibility metadata. ``None`` roots mean whole-molecule
enumeration; an explicitly empty root vector enumerates no paths. Invalid
bounds, widths, count lengths, masks, and roots return ``FingerprintError``.
The complete ChEMBL 37 audit covers all 2,897,804 mutually parseable records
across 18 profiles and 52,160,472 exact comparisons with zero mismatches.
Pinned RDKit's unrooted linear branch can consume atom indices as bond indices
and terminate the process; COSMolKit deliberately uses the documented
bond-path semantics instead of reproducing that crash.

## Conformer Generation And Force Field Applications

Native conformer generation uses RDKit-aligned distance-geometry parameters.
The default value-style molecule operation uses ETKDGv3 and returns a new
molecule value. Multi-conformer generation supports deterministic seeded runs,
RMS pruning, and sequential seed expansion:

```rust
use cosmolkit::{EmbedParameters, Molecule};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let molecule = Molecule::from_smiles("CC(=O)NC")?.with_hydrogens()?;

    let embedded = molecule.with_3d_conformer()?;
    println!("{}", embedded.conformers_3d().len());

    let mut params = EmbedParameters::etkdg();
    params.random_seed = 123;
    params.num_threads = 1;
    params.prune_rms_thresh = 0.5;

    let pruned = molecule.with_3d_conformers_with_params(5, params)?;
    println!("{}", pruned.conformers_3d().len());
    Ok(())
}
```

Force-field APIs operate on molecules with existing 3D conformers and return
new molecule values, so the input coordinates are left unchanged.

```rust
use cosmolkit::{
    Molecule, mmff_has_all_molecule_params, mmff_optimize_molecule,
    uff_has_all_molecule_params, uff_optimize_molecule,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let molecule = Molecule::from_smiles("CCO")?.with_hydrogens()?.sanitize()?;

    let mut builder = molecule.to_builder();
    builder.add_3d_conformer(vec![
        [0.000, 0.000, 0.000],
        [1.540, 0.000, 0.000],
        [2.100, 1.200, 0.000],
        [-0.600, 0.900, 0.000],
        [-0.600, -0.900, 0.000],
        [0.000, 0.000, 1.000],
        [1.900, -0.900, 0.000],
        [1.700, 0.000, 1.000],
        [2.900, 1.200, 0.000],
    ])?;
    let molecule = builder.build()?;

    if uff_has_all_molecule_params(&molecule)? {
        let result = uff_optimize_molecule(&molecule, 200, 10.0, -1, true)?;
        println!("UFF energy: {:.6}", result.energy);
    }

    if mmff_has_all_molecule_params(&molecule)? {
        let result = mmff_optimize_molecule(&molecule, "MMFF94", 200, 100.0, -1, true)?;
        println!("MMFF94 needs_more: {}", result.needs_more);
    }

    Ok(())
}
```

## Molecular Alignment And RMSD

Ordinary molecular alignment uses the source-backed RDKit MolAlign boundary.
Transform queries, best RMSD, coordinate-frame RMSD, and all-conformer pair
measurements are read-only. Coordinate changes are exposed only through a
value-style method or an explicit trailing-underscore method.

```rust
use cosmolkit::{AlignmentAtomMap, AlignmentParameters, Molecule};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let reference = Molecule::from_smiles("CCC")?.with_only_3d_conformer(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
        true,
    )?;
    let probe = Molecule::from_smiles("CCC")?.with_only_3d_conformer(
        vec![[3.0, -2.0, 1.0], [4.0, -2.0, 1.0], [3.0, 0.0, 1.0]],
        true,
    )?;
    let params = AlignmentParameters {
        atom_map: Some(
            (0..3)
                .map(|index| AlignmentAtomMap {
                    probe_atom: index,
                    reference_atom: index,
                })
                .collect(),
        ),
        ..Default::default()
    };

    let measured = probe.alignment_transform_to(&reference, &params)?;
    let (aligned, applied) = probe.with_alignment_to(&reference, &params)?;
    assert_eq!(probe.conformers_3d()[0].coordinates()[0], [3.0, -2.0, 1.0]);
    assert!(measured.rmsd < 1.0e-8 && applied.rmsd < 1.0e-8);
    assert_eq!(aligned.conformers_3d()[0].coordinates()[0], [0.0, 0.0, 0.0]);
    Ok(())
}
```

Weighted and reflected alignment, automatic or explicit atom maps, conformer
IDs, iteration limits, best-map selection, and conformer-set alignment use
typed parameter objects. O3A and MMFF/Crippen scoring are separate capabilities
and are not implied by this ordinary MolAlign API.

## Examples

```bash
cargo run -p cosmolkit --example smiles_write_options
cargo run -p cosmolkit --example draw_svg
cargo run -p cosmolkit --example draw_png
cargo run -p cosmolkit --example sdf_to_smiles
cargo run -p cosmolkit --example protein_from_pdb
cargo run -p cosmolkit --example read_xyz
cargo run -p cosmolkit --example molalign_rmsd
cargo run -p cosmolkit --example conformer_generation
cargo run -p cosmolkit --example forcefield_optimization
```

## Contributor Validation

Core validation should use operation-contract checks:

```bash
cargo check -p cosmolkit-core --features op-contracts-strict
cargo test -p cosmolkit-core --release --features op-contracts-strict
cargo check -p cosmolkit-py
cargo fmt --all
```

Use debug-profile test filters for small local iterations. Use release mode with
`op-contracts-strict` for large local runs, parity suites, and CI; release-mode
testing keeps operation contracts and runtime invariants enabled through the
strict feature set.

Python binding validation:

```bash
uv sync --group dev
.venv/bin/maturin develop --manifest-path python/Cargo.toml
.venv/bin/pytest
```

The planned crate split is a staged internal architecture migration, not the
current workspace layout. Its target is for `cosmolkit` to own the public
`Molecule` API and operation runtime while shared model values and source-backed
algorithms move behind explicit crate boundaries. It is not intended to break
the existing supported external API: users should continue importing from
`cosmolkit`, without changing normal molecule workflows merely because an
implementation moves between crates. Internal implementation crates will not
become competing public `Molecule` entry points. Any unavoidable public change
will be handled separately through the normal versioning and deprecation
policy.
