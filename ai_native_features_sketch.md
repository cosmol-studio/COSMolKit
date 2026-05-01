# COSMolKit AI-Native Feature Sketch

## Goal

Add a set of AI-native molecular representations and utilities on top of the RDKit-compatible chemistry core.

The goal is not to replace the parity baseline. Instead, COSMolKit should use the verified core as a reliable foundation for model-ready graph, geometry, tokenization, and diffusion workflows that are not first-class concepts in traditional cheminformatics APIs.

These features should make COSMolKit useful for:

- molecular graph neural networks
- 3D molecular transformers
- diffusion and flow-matching conformer models
- autoregressive molecular generation
- molecular editing and interpolation
- compact conformer and geometry representation

## Design Principles

- Keep RDKit parity work isolated from AI-specific abstractions.
- Preserve COSMolKit's deterministic transformation model: normal molecule operations return new objects and do not mutate inputs.
- Prefer explicit, versioned feature schemas such as `cosmol-v1` over ad hoc arrays.
- Preserve chirality and stereochemistry semantics through every exported representation.
- Make tensor export deterministic and reproducible across Rust and Python.
- Expose Python APIs that are immediately usable with PyTorch, JAX, NumPy, and graph-learning libraries.
- Treat internal coordinates and torsions as first-class molecular data, not only derived helper values.

## 1. AI-Native Graph Export

COSMolKit should provide model-ready graph export instead of requiring users to manually assemble atom and bond arrays.

Example API:

```python
data = mol.to_graph(
    node_features="cosmol-v1",
    edge_features="geometry",
    include_pos=True,
    include_chirality=True,
    include_torsion=True,
)
```

Possible exported fields:

- node feature matrix with atomic number, degree, charge, valence, aromaticity, hybridization, ring membership, hydrogen counts, isotope, radical state, chirality, and CIP information
- edge index and edge feature matrix with bond order, aromaticity, conjugation, stereo, direction, ring membership, and optional geometric distance
- optional 2D or 3D coordinates
- optional torsion indices and torsion features
- optional fragment, ring, and rotatable-bond annotations

Potential Python output adapters:

- plain dictionaries of NumPy arrays
- PyTorch Geometric `Data`
- DGL graph objects
- JAX-friendly array dictionaries

Schema versioning should be explicit:

```python
graph = mol.to_graph(node_features="cosmol-v1")
graph.schema
```

This makes downstream models reproducible when feature definitions evolve.

## 2. Internal Coordinate System

COSMolKit should expose a unified internal-coordinate layer for chemistry-aware geometry manipulation.

Example API:

```python
internal = mol.to_internal_coordinates()
cart = internal.to_cartesian()
```

Supported representations should include:

- Z-matrix export and import
- bond-angle-torsion trees
- ring-aware internal coordinates
- torsion graphs
- rotatable bond decomposition
- fragment-local coordinate frames

This can support:

- conformer generation
- diffusion models
- molecular editing
- conformer interpolation
- conformer compression
- 3D tokenization

Internal coordinates should be ring-aware. Ring closures, constrained torsions, and fused-ring systems need explicit handling instead of being silently flattened into an acyclic tree.

## 3. Torsion and Chirality-Aware Diffusion Utilities

COSMolKit can provide reusable components for diffusion or flow-matching models on molecular geometry.

Example API:

```python
noise = cosmolkit.diffusion.sample_torsion_noise(torsions, t)
x_t = cosmolkit.diffusion.q_sample_internal(x_0, t)
loss = cosmolkit.diffusion.torsion_loss(pred, target)
```

Important utilities:

- periodic angle loss
- sin/cos torsion encoding
- chirality-preserving reconstruction
- ring torsion constraints
- rotatable-bond masks
- torsion noise schedules
- internal-coordinate perturbation helpers

Torsion loss should account for periodicity. A naive squared angle difference is incorrect around the `-pi` / `pi` boundary:

```python
loss = (pred_angle - true_angle) ** 2
```

Instead, periodic loss can use:

```python
loss = 1 - cos(pred_angle - true_angle)
```

or an equivalent stable vector formulation over `sin` and `cos` encodings.

Reconstruction utilities should explicitly validate whether tetrahedral chirality and constrained ring geometry are preserved after denoising or coordinate conversion.

## 4. Molecular Tokenization for AI

COSMolKit should support tokenization schemes beyond SMILES strings.

Example API:

```python
tokens = mol.to_tokens(
    scheme="cosmol-3d",
    include_torsion=True,
    include_chirality=True,
)
```

Potential token families:

- SMILES tokens
- atom tokens
- bond tokens
- fragment tokens
- torsion tokens
- 3D geometry tokens
- pharmacophore tokens
- ring and scaffold tokens

These tokenization schemes can support:

- molecular language models
- 3D molecular transformers
- diffusion transformers
- autoregressive conformer generation
- edit-based molecular generation

Tokenizers should be versioned and reversible where possible:

```python
tokens = mol.to_tokens(scheme="cosmol-3d-v1")
mol2 = cosmolkit.Molecule.from_tokens(tokens, scheme="cosmol-3d-v1")
```

Reversibility may be exact for graph-level tokens and approximate or constrained for discretized geometry tokens. The API should make that distinction explicit.

## Suggested Implementation Stages

### Stage 1: Stable Graph Export

- Define `GraphExport` data structures in Rust.
- Add `cosmol-v1` node and edge feature schemas.
- Expose Python dictionary / NumPy export first.
- Add deterministic tests for feature ordering and shape stability.

### Stage 2: Torsion and Internal Coordinate Primitives

- Add rotatable-bond detection.
- Add torsion enumeration with atom-index quartets.
- Add bond-angle-torsion trees for acyclic fragments.
- Add explicit unsupported-path errors for complex ring cases until ring-aware support lands.

### Stage 3: Geometry-Aware ML Helpers

- Add periodic angle loss utilities.
- Add sin/cos torsion encoders.
- Add chirality-preservation checks.
- Add ring torsion masks and constraint metadata.

### Stage 4: Tokenization Schemes

- Add versioned token schema definitions.
- Start with graph and torsion tokens.
- Add optional geometry quantization.
- Provide roundtrip tests for reversible schemes.

## Interaction with RDKit Compatibility

RDKit compatibility remains the correctness floor for chemistry perception, file IO, and graph semantics.

AI-native features should build on that floor but do not need to imitate RDKit APIs. Where RDKit behavior determines chemical facts such as aromaticity, valence, chirality, or canonical traversal, COSMolKit should keep using source-level parity tests. Where COSMolKit introduces model-specific abstractions, the project should define its own stable schemas and tests.

The intended split is:

- RDKit-compatible core for chemical correctness
- COSMolKit-native AI layers for graph tensors, internal coordinates, torsions, diffusion helpers, and tokenization

This gives COSMolKit a distinct role: a strict, reliable chemistry engine that also exposes modern representations designed for machine learning workflows.
