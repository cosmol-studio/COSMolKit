Fingerprints
============

.. meta::
   :description: COSMolKit molecular fingerprint APIs for source-backed Morgan, Pattern, AtomPair, Topological Torsion, MACCS, RDKFingerprint, Avalon, and Layered bit and count vectors with explicit RDKit parity boundaries.

COSMolKit exposes fixed-length bit vectors plus the source-defined sparse bit
and count forms used by AtomPair. The exposed Morgan and MACCS branches are
covered by strict RDKit bit-identical parity tests. The source-backed
topological, Avalon, Pattern, Layered, AtomPair, and Topological Torsion
implementations follow their pinned upstream algorithms and exact
maintained-corpus validation.
Similarity-shape correlation or structurally similar hashing is not a
compatibility claim. The Python ``Fingerprint`` object is a sparse view over
the binary vector: ``on_bits()`` returns the bit indexes whose value is 1. It
is not a dense floating-point neural embedding.

Scaffold molecules
------------------

These value-returning transforms leave the original molecule unchanged:

.. code-block:: python

   import cosmolkit as ck

   mol = ck.mol_from_smiles("O=C1CCCCC1")
   mol.murcko_scaffold().to_smiles()   # "C1CCCCC1"
   mol.net_scaffold().to_smiles()      # "*=C1CCCCC1"
   mol.murcko_decompose().to_smiles()  # "O=C1CCCCC1"

``murcko_scaffold`` follows RDKit MolHash ``MurckoScaffold``: repeatedly remove
terminal atoms and transfer their bond-order contribution to hydrogen counts.
``net_scaffold`` follows MolHash ``ExtendedMurcko``: retain ring/linker atoms
and replace their directly attached substituents with dummy atoms.
``murcko_decompose`` follows ChemTransforms ``MurckoDecompose``: retain rings,
their shortest-path linkers, and atoms double-bonded to retained atoms.

The MolHash transforms require calculated hydrogen counts, and
``murcko_decompose`` requires initialized ring information for nonempty input.
Default SMILES parsing supplies both. Missing source state raises a typed
operation error; these methods do not silently sanitize unprepared input.

Topological Torsion fingerprints
--------------------------------

Topological Torsion is a separate fingerprint family from
``Molecule.fingerprint_topological()``. The latter is RDKit's path/subgraph
``RDKFingerprintMol`` algorithm; Topological Torsion enumerates ordered atom
paths of a configured length and encodes or hashes their atom invariants. The
two names, parameters, and outputs are intentionally not interchangeable.

The modern generator exposes all four RDKit vector forms. Sparse count output
retains unfolded 64-bit torsion ids, sparse bit output applies the generator's
bit-domain and count-simulation rules, count output folds counts into
``fp_size``, and bit output returns an explicit ``Fingerprint`` of that size:

.. code-block:: python

   import cosmolkit as ck

   mol = ck.mol_from_smiles("CCCCO")
   settings = ck.TopologicalTorsionParams(
       include_chirality=False,
       torsion_atom_count=4,
       count_simulation=True,
       count_bounds=[1, 2, 4, 8],
       fp_size=2048,
   )
   generator = ck.TopologicalTorsionFingerprintGenerator(params=settings)

   sparse_count = mol.fingerprint_topological_torsion_sparse_count_with_generator(generator)
   sparse_bit = mol.fingerprint_topological_torsion_sparse_with_generator(generator)
   count = mol.fingerprint_topological_torsion_count_with_generator(generator)
   bit = mol.fingerprint_topological_torsion_with_generator(generator)

   print(sparse_count.nonzero_elements())
   print(sparse_bit.on_bits())
   print(count.nonzero_elements())
   print(bit.on_bits())

The ``TopologicalTorsionParams`` passed to the constructor controls
``include_chirality``, ``torsion_atom_count``, ``count_simulation``,
``count_bounds``, ``fp_size``, ``bits_per_feature``, and
``only_shortest_paths``. ``generator.settings()`` returns a live settings view
whose writable fields update the generator. Scalar calls accept a
``TopologicalTorsionCallParams`` object with ``from_atoms``,
``ignore_atoms``, and ``custom_atom_invariants``. ``from_atoms`` selects paths
by endpoint, while ``ignore_atoms`` excludes paths containing any selected
atom. ``only_shortest_paths`` keeps only source-defined shortest atom paths.
Custom invariants replace the default atom-code generator for that call.

The default atom-code generator reads the molecule's explicit-valence cache,
as does RDKit. Molecules returned by the default ``ck.mol_from_smiles()``
path and by ``with_hydrogens()`` are ready for fingerprinting because both
operations commit the source-defined explicit-valence state. A molecule parsed
with ``sanitize=False`` or produced by another operation that invalidates
valence state must be sanitized first, for example with
``mol = mol.sanitize()``, unless the call supplies
``custom_atom_invariants``. Missing cache state is reported as a typed
exception rather than being recomputed implicitly inside the fingerprint call.

Topological Torsion provenance uses the shared ``FingerprintAdditionalOutput`` container.
Allocate only the outputs needed before the call. This family populates
``atom_to_bits``, ``atom_counts``, ``bit_paths``, and ``atoms_per_bit``;
``bit_info_map`` belongs to other generator families and remains empty here,
matching RDKit. Parity tests still request and compare that empty container so
the absence of entries is verified rather than assumed.

.. code-block:: python

   output = ck.FingerprintAdditionalOutput()
   output.allocate_atom_to_bits()
   output.allocate_atom_counts()
   output.allocate_bit_paths()
   output.allocate_atoms_per_bit()

   mol.fingerprint_topological_torsion_with_generator(generator, output=output)
   print(output.atom_to_bits())
   print(output.atom_counts())
   print(output.bit_paths())
   print(output.atoms_per_bit())

Bulk generator methods preserve input order and accept ``num_threads``:

.. code-block:: python

   molecules = [
       ck.mol_from_smiles("CCCC"),
       ck.mol_from_smiles("CCCCC"),
   ]
   bits = generator.fingerprints(molecules, num_threads=2)

   assert bits[0].on_bits() == molecules[0].fingerprint_topological_torsion_with_generator(generator).on_bits()
   assert bits[1].on_bits() == molecules[1].fingerprint_topological_torsion_with_generator(generator).on_bits()

Additional vector forms have their own parameter semantics:

- ``fingerprint_topological_torsion_sparse_count_legacy()`` returns the unfolded
  sparse-count vector and preserves RDKit's ``2**width - 1`` reported-size
  compatibility behavior.
- ``fingerprint_topological_torsion_count_legacy()`` returns folded sparse
  counts.
- ``fingerprint_topological_torsion_legacy()`` preserves the
  ``bits_per_entry`` block sizing and its distinct four-slot
  ``[1, 2, 4, 8]`` versus non-four-slot threshold rule.

All adapters terminate in the same Rust atom-code, path,
torsion-code/hash, and shared fingerprint-vector core. They do not call RDKit
at runtime. Invalid sizes, unsupported source-undefined shifts, malformed
JSON, bad selections, and incompatible custom invariants fail visibly as
typed Python exceptions instead of selecting another fingerprint algorithm.
Rust exposes the corresponding ``FingerprintError``; Rust ``MoleculeBatch``
conveniences additionally report computation failures as indexed
``BatchValidationError`` records while retaining existing invalid inputs as
``None`` in their original positions.

Validation covers vector forms, provenance, JSON, selections, chirality,
paths, counts, collisions, and errors against the pinned RDKit reference.
See ``VALIDATION.md`` in the source repository for the validation scope.

Source-backed topological and Avalon fingerprints
-------------------------------------------------

``Molecule.fingerprint_topological()`` and
``Molecule.fingerprint_avalon()`` execute source-backed Rust implementations.
The standard corpus workflow compares Avalon against pinned RDKit with one
fixed-seed parameter combination per molecule.

The two APIs are source-backed and return a fresh explicit bit vector:

- ``fingerprint_topological()`` requires the complete RDKit
  ``RDKFingerprintMol``/RDKitFP generator behavior, including branched-path
  enumeration, source random-bit generation, density folding, atom
  invariants, and the exposed path and atom-selection parameters.
- ``fingerprint_avalon()`` follows the Avalon/reaccs bit-vector path,
  including ``bitFlags``, ``isQuery``, hydrogen handling, aromaticity passes,
  and byte-rounded vector semantics. ``resetVect`` is an internal adapter
  detail and is not exposed on COSMolKit's value-returning API.

Comparison uses exact bit vectors, not similarity correlation or partial bit
agreement. Python defaults to ``0xF07FFF``, matching RDKit's Python wrapper.
MOL conversion leaves the original molecule unchanged. The fingerprint path
handles molecules without stored coordinates and does not require callers to
generate 2D coordinates first.

.. code-block:: python

   mol = ck.mol_from_smiles("c1ccccc1O")

   topological = mol.fingerprint_topological(
       min_path=1,
       max_path=7,
       fp_size=2048,
       num_bits_per_feature=2,
   )
   avalon = mol.fingerprint_avalon(
       n_bits=512,
       is_query=False,
       bit_flags=0xF07FFF,
   )

   print(topological.on_bits())
   print(avalon.on_bits())

Connected fragments
-------------------

.. code-block:: python

   mol = ck.mol_from_smiles("CC.O.[Na+]")
   fragments = mol.fragments()
   print([part.to_smiles() for part in fragments])  # ["CC", "O", "[Na+]"]
   print(mol.largest_fragment().to_smiles())        # "CC"

Fragments are sanitized, retain source component order and copy their
conformer rows. Both methods return new molecules without changing ``mol``.
``largest_fragment()`` chooses the most atoms and the last component on ties;
it is not RDKit MolStandardize's chemically ranked fragment chooser. Empty
input returns an empty fragment list; requesting its largest fragment raises
a typed operation error.

Layered fingerprints
--------------------

``Molecule.fingerprint_layered()`` implements RDKit's Layered fingerprint
algorithm version ``0.7.0``. RDKit labels this algorithm experimental;
COSMolKit preserves that upstream metadata even though the exposed parameter
and result types are ordinary project-native APIs. The six active source layer
bits are:

.. list-table::
   :header-rows: 1

   * - Flag
     - Layer
     - Encoded state
   * - ``0x01``
     - topology
     - path-local atom degree and bond-neighbor topology
   * - ``0x02``
     - bond order
     - normalized bond type plus topology
   * - ``0x04``
     - atom type
     - atomic numbers plus topology
   * - ``0x08``
     - ring presence
     - sparse ring-bond membership
   * - ``0x10``
     - ring size
     - minimum source SSSR ring size
   * - ``0x20``
     - aromaticity
     - query-aware endpoint aromaticity plus topology

``0x07`` is the source substructure prefix and ``0x3f`` selects all active
layers. The default ``0xffffffff`` retains all ten source flag slots; high bits
have no encoder in RDKit and therefore add no components. They are accepted,
not repurposed as COSMolKit extensions.

Paths are bond paths with inclusive ``min_path`` and ``max_path`` bounds. The
defaults are 1 through 7, 2,048 output bits, and branched enumeration.
``branched_paths=False`` selects linear bond paths. ``from_atoms=None`` means
no root restriction. An explicit empty list is different: it is a present
empty root selection and returns an empty fingerprint. Duplicate and reordered
roots retain the source aggregation order.

.. code-block:: python

   import cosmolkit as ck

   molecule = ck.mol_from_smiles("c1ccccc1O")
   substructure = molecule.fingerprint_layered(
       layers=0x07,
       min_path=1,
       max_path=7,
       fp_size=2048,
   )

   even_mask = ck.Fingerprint.from_on_bits(257, range(0, 257, 2))
   counted = molecule.fingerprint_layered_with_output(
       layers=0x3F,
       min_path=2,
       max_path=4,
       fp_size=257,
       atom_counts=[10] * molecule.num_atoms(),
       set_only_bits=even_mask,
       branched_paths=True,
       from_atoms=[0],
   )

   print(substructure.on_bits())
   print(counted.fingerprint().on_bits())
   print(counted.atom_counts())

``set_only_bits`` must have exactly ``fp_size`` bits. It gates projected bits
before insertion and before counting. A supplied ``atom_counts`` vector is a
seed, not a zeroed output buffer; its length must be at least the molecule atom
count. Every atom in a path is incremented once when any layer projection for
that path passes the mask, even when several layers set bits or bit collisions
occur. Omitting counts makes ``atom_counts()`` return ``None``.

Invalid zero/reversed path bounds, zero width, short count vectors,
width-mismatched masks, out-of-range roots, and ring-preparation failures are
reported as typed exceptions. The operation is read-only: repeated,
interleaved, batch, and concurrent calls share the same Rust scalar core and
do not mutate molecule topology, properties, coordinates, or derived caches.

Pinned RDKit contains one upstream defect in the unrooted linear call: it
requests atom-index paths and later consumes them as bond indices, which can
terminate the process for ordinary acyclic molecules. COSMolKit does not copy
that crash. It applies the source header's documented bond-path semantics, the
same semantics used by the valid rooted linear branch. This is an explicit
process-safety compatibility difference, not a chemistry fallback.

Pattern fingerprints
--------------------

``Molecule.fingerprint_pattern()`` reproduces RDKit's ordinary-molecule
Pattern fingerprint using the fixed, source-ordered table of 13 SMARTS
queries. The default result has 2,048 bits. ``tautomeric=True`` enables the
source's tautomer-aware structural hash for single, double, and aromatic bond
states:

.. code-block:: python

   mol = ck.mol_from_smiles("c1ccccc1O")

   pattern = mol.fingerprint_pattern()
   tautomeric_pattern = mol.fingerprint_pattern(
       n_bits=2048,
       tautomeric=True,
   )

   print(pattern.on_bits())
   print(tautomeric_pattern.on_bits())

The call is value-style and does not mutate the molecule. Query-bearing
molecules are supported and follow the source's Pattern-specific atom and bond
suppression rules. The compiled SMARTS table is shared across repeated,
concurrent, scalar, and batch calls; Pattern calls can be interleaved with the
other fingerprint families without shared option state.

``n_bits`` must be greater than zero. SMARTS-compilation and substructure
matching failures remain typed exceptions rather than selecting another
fingerprint algorithm. RDKit labels Pattern fingerprint version ``1.0.0``
experimental; COSMolKit preserves that upstream metadata while its documented
ordinary-molecule boundary is tested exactly against the pinned source.

The pinned RDKit ordinary overload accepts ``atomCounts`` and ``setOnlyBits``
but only validates their sizes; their values are otherwise unused and
unchanged. COSMolKit intentionally omits those inert arguments instead of
presenting them as functional options. RDKit's separate ``MolBundle`` overload
intersects member fingerprints. COSMolKit does not expose it because an
ordered ``MoleculeBatch`` returns one result per record and is not a bundle.

The committed focused, small, and 5,000-row matrices cover query inputs,
tautomer mode, boundary and collision-prone widths, inert-argument behavior,
and exact validation errors. The complete ChEMBL 37 audit compares ten full
profiles across all 2,897,804 mutually parseable records: 28,978,040 exact
vectors with zero mismatch.

AtomPair fingerprints
---------------------

``Molecule`` exposes the four source generator result forms through one
project-native implementation:

- ``fingerprint_atom_pair_sparse_count()`` returns the source-width count map;
- ``fingerprint_atom_pair_count()`` returns the folded count map;
- ``fingerprint_atom_pair_sparse()`` returns the source-width bit set;
- ``fingerprint_atom_pair()`` returns the fixed-width explicit bit vector.

The shared parameter surface covers 2D topological or 3D conformer distances,
minimum and maximum distance, chirality, count simulation and bounds,
``num_bits_per_feature``, ``from_atoms``, ``ignore_atoms``, conformer selection,
and custom atom invariants. A 3D call requires a conformer and sets
``use_2d=False`` explicitly.

.. code-block:: python

   mol = ck.mol_from_smiles("CCCO")

   explicit = mol.fingerprint_atom_pair()
   sparse_count = mol.fingerprint_atom_pair_sparse_count()
   folded_count = mol.fingerprint_atom_pair_count()
   sparse_bits = mol.fingerprint_atom_pair_sparse()

   print(explicit.on_bits())
   print(sparse_count.nonzero_elements())
   print(folded_count.nonzero_elements())
   print(sparse_bits.on_bits())

Pass a ``FingerprintAdditionalOutput`` collector to
``fingerprint_atom_pair_with_params()`` for exact atom counts, atom-to-bit
rows, endpoint-pair bit information, and atoms-per-bit provenance:

.. code-block:: python

   output = ck.FingerprintAdditionalOutput()
   output.allocate_atom_counts()
   output.allocate_atom_to_bits()
   output.allocate_bit_info_map()
   output.allocate_atoms_per_bit()
   fp = mol.fingerprint_atom_pair_with_params(ck.AtomPairFingerprintParams(), output)

   print(fp.on_bits())
   print(output.atom_counts())
   print(output.atom_to_bits())
   print(output.bit_info_map())
   print(output.atoms_per_bit())

The complete ChEMBL 37 validation compares 40 vectors and one full provenance
output for every one of the 2,897,804 mutually parseable records: 118,809,964
exact comparisons with zero mismatch. The maintained 5,000-row, ten-profile
matrix remains the committed continuous regression gate.

Topological provenance
-----------------------

``fingerprint_topological_with_output()`` returns a
``TopologicalFingerprintResult``. Request ``atom_bits`` and/or ``bit_info`` to
receive the matching source provenance outputs:

.. code-block:: python

   result = mol.fingerprint_topological_with_output(
       fp_size=2048,
       atom_bits=True,
       bit_info=True,
   )

   print(result.fingerprint().on_bits())
   print(result.atom_bits())
   print(result.bit_info())

Single Molecules
----------------

.. code-block:: python

   import cosmolkit as ck

   mol = ck.mol_from_smiles("c1ccccc1O")
   fp = mol.fingerprint_morgan(radius=2, fp_size=2048)

   print(fp.n_bits())
   print(fp.on_bits())

Tanimoto similarity is computed directly on ``Fingerprint`` values:

.. code-block:: python

   phenol = ck.mol_from_smiles("c1ccccc1O").fingerprint_morgan()
   benzene = ck.mol_from_smiles("c1ccccc1").fingerprint_morgan()

   print(phenol.tanimoto(benzene))

Additional Output
-----------------

Pass a ``FingerprintAdditionalOutput`` collector to ``fingerprint_morgan()``
to obtain provenance alongside the returned fingerprint:

.. code-block:: python

   output = ck.FingerprintAdditionalOutput()
   output.allocate_atom_counts()
   output.allocate_atom_to_bits()
   output.allocate_bit_info_map()
   output.allocate_atoms_per_bit()
   fp = mol.fingerprint_morgan(
       radius=2, fp_size=2048, additional_output=output,
   )

   print(fp.on_bits())
   print(output.atom_counts())
   print(output.atom_to_bits())
   print(output.bit_info_map())
   print(output.atoms_per_bit())

Supported Parameters
--------------------

The short Morgan methods accept a parameter record or configuration keywords:

.. code-block:: python

   fp = mol.fingerprint_morgan(radius=2, fp_size=256, include_chirality=True)

   params = ck.MorganFingerprintParams(from_atoms=[0])
   params.generator.radius = 2
   params.generator.fp_size = 256
   params.from_atoms.append(1)
   fp = mol.fingerprint_morgan(params)

Nested parameters and configuration lists are live views of their parent
record. Edits go through native conversion and validation; rejected edits leave
the record unchanged. Constructor inputs are copied, so editing a parameter
record does not modify a list or generator previously passed into it. Do not
combine a parameter record with its configuration keywords, or ``generator=``
with generator keywords such as ``radius=``.

The Python binding exposes the source-backed Morgan generator branches covered
by exact RDKit bit parity:

- ``radius`` and ``fp_size``
- ``include_chirality`` and ``use_bond_types``
- ``count_simulation`` and ``count_bounds``
- ``only_nonzero_invariants``
- ``include_redundant_environments``
- ``from_atoms`` and ``ignore_atoms``
- ``custom_atom_invariants`` and ``custom_bond_invariants``
- ``atom_invariants_generator="connectivity" | "morgan" | "feature" | "fcfp"``
- ``atom_invariants_include_ring_membership``
- ``bond_invariants_generator="morgan" | "default" | "bond"``
- ``bond_invariants_use_bond_types``
- ``bond_invariants_use_chirality``
- ``num_bits_per_feature``

The list is the tested boundary, not a claim that every RDKit fingerprint
generator or every input-state preparation branch is implemented. Unsupported
branches propagate an error; they do not silently fall back to another
fingerprint algorithm.

Batch Fingerprints
------------------

``MoleculeBatch`` exposes matching batch APIs. Invalid records kept with
``errors="keep"`` produce ``None`` in the corresponding output position.

.. code-block:: python

   import cosmolkit as ck

   batch = ck.mols_from_smiles_list(
       ["CCO", "not-smiles", "CCCO"],
       errors="keep",
   ).with_parallel_jobs(8)
   fps = batch.fingerprint_morgan_list(fp_size=2048)
   atom_pairs = batch.fingerprint_atom_pair_list(
       fp_size=2048,
       n_jobs=8,
       progress_bar=False,
   )
   layered = batch.fingerprint_layered_list(
       layers=0x3F,
       min_path=1,
       max_path=7,
       fp_size=2048,
       n_jobs=8,
       progress_bar=False,
   )

   print([fp.on_bits() if fp is not None else None for fp in fps])
   print([fp.on_bits() if fp is not None else None for fp in atom_pairs])
   print([fp.on_bits() if fp is not None else None for fp in layered])
