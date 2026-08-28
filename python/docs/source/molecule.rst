Molecule Values
===============

.. meta::
   :description: COSMolKit molecular graph APIs for value semantics, SMILES and SMARTS, stereochemistry, coordinates, conformers, substructure search, editing, InChI, and binary serialization.

``Molecule`` objects behave as value-style molecule values. Transformation
methods return new molecule objects and leave the original object unchanged.
Internally COSMolKit uses copy-on-write (COW) storage to share unchanged data
efficiently. This is intentionally different from common RDKit Python
workflows, where code often mutates an existing molecule or ``RWMol`` directly.

.. code-block:: python

   from cosmolkit import Molecule

   mol = Molecule.from_smiles("CCO")
   mol_h = mol.with_hydrogens()

   assert mol is not mol_h
   print(mol.to_smiles())
   print(mol_h.to_smiles())

Do not write code that assumes ``mol.with_hydrogens()`` changes ``mol``. Keep
the returned value and pass that value to later operations.

Common transformations include:

- ``with_hydrogens()``
- ``without_hydrogens()``
- ``with_kekulized_bonds()``
- ``with_2d_coordinates()``
- ``with_chiral_tags_from_structure()``

Elements and Periodic-Table Metadata
------------------------------------

``Element`` is a stable integer enum whose values are atomic numbers. Zero is
the dummy atom, and 1 through 118 represent H through Og. Symbol lookup is
case-sensitive and metadata comes from the same source-aligned periodic table
used by the chemistry core:

.. code-block:: python

   from cosmolkit import Element, element_from_symbol, get_element_info

   assert Element.C == 6
   assert element_from_symbol("Cl") == Element.CL

   chlorine = get_element_info(Element.CL)
   assert chlorine.symbol() == "Cl"
   assert chlorine.atomic_number() == 17

``ELEMENT_MAP`` maps canonical symbols, the dummy symbol ``"*"``, and retained
source aliases to ``Element`` values. Molecule atoms expose the same identity
through ``atom.element()``.

Read Finalization
-----------------

Molfile and SDF molecule readers use the RDKit-source-backed finalization path
for modeled parser behavior. With the default ``sanitize=True`` and
``remove_hs=True``, readers parse the CTAB, process molfile/SDF properties,
assign modeled stereochemistry, remove hydrogens through the RDKit-aligned
hydrogen-removal path, sanitize, and assign final stereochemistry.

Passing ``sanitize=False`` preserves the parsed molecule state for a later
value-style sanitize operation:

.. code-block:: python

   raw = Molecule.read_mol("input.mol", sanitize=False)
   sanitized = raw.sanitize()

   assert raw is not sanitized

Passing ``remove_hs=False`` preserves explicit hydrogens for later value-style
hydrogen removal:

.. code-block:: python

   with_h = Molecule.read_sdf("input.sdf", remove_hs=False)
   heavy = with_h.without_hydrogens()

   assert with_h is not heavy

The same delayed-operation pattern applies to ``read_mol_from_str()`` and
``read_sdf_from_str()``.

In-Place Operations
-------------------

Performance-sensitive code can opt into explicit in-place mutation. Every
public ``Molecule`` in-place method ends with ``_``; the trailing underscore has
no other ``Molecule`` API meaning.

.. code-block:: python

   mol = Molecule.from_smiles("CCO")
   mol.add_hydrogens_()
   mol.compute_2d_coordinates_()

Common in-place operations include:

- ``add_hydrogens_()``
- ``remove_hydrogens_()``
- ``kekulize_()``
- ``sanitize_()``
- ``compute_2d_coordinates_()``
- ``assign_chiral_tags_from_structure_()``

If an in-place method returns an error, the molecule is not guaranteed to equal
its pre-call value. Use the value-style method when failure-preserving behavior
is required.

SMILES Output
-------------

``to_smiles()`` returns a SMILES string:

.. code-block:: python

   mol = Molecule.from_smiles("F[C@H](Cl)Br")

   print(mol.to_smiles())
   print(mol.to_smiles(isomeric_smiles=False))

SMILES writer options are available on both single molecules and batches:

.. code-block:: python

   benzene = Molecule.from_smiles("c1ccccc1")
   ethanol = Molecule.from_smiles("CCO")

   print(benzene.to_smiles(kekule=True))
   print(ethanol.to_smiles(all_bonds_explicit=True))
   print(ethanol.to_smiles(canonical=False, rooted_at_atom=2))

Explicit Editing
----------------

Use ``Molecule.edit()`` when you want to stage changes and commit them as one
new molecule:

.. code-block:: python

   editor = mol.edit()
   cl = editor.add_atom("Cl")
   editor.add_bond(0, cl, order="single")

   edited = editor.commit()

Serialization
-------------

``Molecule`` supports Python ``pickle`` for in-process persistence and
inter-process transfer. The pickle state carries a COSMolKit pickle schema and
the versioned core molecule archive payload, so future incompatible payload
changes can be rejected explicitly instead of being decoded as the wrong
structure.

.. code-block:: python

   import pickle

   mol = Molecule.from_smiles("F[C@H](Cl)[13CH3:7]").with_2d_coordinates()
   restored = pickle.loads(pickle.dumps(mol, protocol=pickle.HIGHEST_PROTOCOL))

   print(restored.to_smiles(canonical=False))

Advanced callers can use ``mol_to_binary()`` and ``mol_from_binary()`` to
inspect or persist the COSMolKit molecule archive directly. Python
applications should prefer ``pickle`` unless they specifically need the raw
archive payload. The versioned archive preserves graph, coordinate, property,
and materialized derived chemistry state so supported operations retain the
same behavior after restoration; legacy archive readers remain available for
older payloads:

.. code-block:: python

   payload = mol.mol_to_binary()
   restored = Molecule.mol_from_binary(payload)

   assert restored.to_smiles(canonical=False) == mol.to_smiles(canonical=False)

Depictions
----------

Molecules with 2D coordinates can be exported as SVG or PNG:

.. code-block:: python

   mol = Molecule.from_smiles("c1ccccc1O").with_2d_coordinates()

   svg = mol.to_svg(width=400, height=300)
   mol.write_svg("python/examples/output/phenol.svg", width=400, height=300)
   mol.write_png("python/examples/output/phenol.png", width=400, height=300)

Stereo
------

COSMolKit keeps the atom-level CW/CCW chiral tag path available. This is the
closest representation to the explicit chiral information carried by SMILES or
RDKit atoms:

.. code-block:: python

   from cosmolkit import ChiralTag, Molecule

   mol = Molecule.from_smiles("F[C@H](Cl)Br")

   for atom in mol.atoms():
       if atom.chiral_tag() != ChiralTag.CHI_UNSPECIFIED:
           print(atom.idx(), atom.chiral_tag().name)

   print(mol.find_chiral_centers(include_unassigned=False))

Atom and bond enum-valued fields return Python ``IntEnum`` members, so callers
can compare or match against ``ChiralTag``, ``BondOrder``, ``BondDirection``,
and ``BondStereo`` instead of spelling chemistry states as strings. Read-only
maps such as ``BOND_ORDER_MAP`` and ``CHIRAL_TAG_MAP`` are available when a
string name from an external source needs to be converted to the enum member.

When code needs COSMolKit's ordered-ligand tetrahedral representation, use
``tetrahedral_stereo()``. The returned ligand order is the stereochemical
value, not just the atom adjacency order. This makes it useful both for
finding the four ligands around a center and for comparing whether two records
represent the same tetrahedral configuration. Equivalent even permutations are
canonicalized to one numeric representative. The precise contract is in
`dev/tetrahedral_stereo.md <https://github.com/cosmol-studio/COSMolKit/blob/main/dev/tetrahedral_stereo.md>`__.

.. code-block:: python

   mol = Molecule.from_smiles("F[C@H](Cl)Br")

   print(mol.tetrahedral_stereo())
   print(Molecule.from_smiles("F[C@@H](Cl)Br").tetrahedral_stereo())
   print(mol.with_hydrogens().tetrahedral_stereo())
   print(Molecule.from_smiles("F[C@](Cl)(Br)I").tetrahedral_stereo())
   print(Molecule.from_smiles("F[C@@](Cl)(Br)I").tetrahedral_stereo())

``None`` in the ligand list represents an implicit hydrogen ligand. It does
not mean the ligand slot is empty. If hydrogens are materialized with
``with_hydrogens()``, that hydrogen ligand is returned as an atom index.

Potential Stereo And Stereoisomer Enumeration
----------------------------------------------

``analyze_potential_stereo()`` returns an isolated molecule state and ordered,
typed potential-stereo records without mutating the source. Each record reports
its stereo type, specified state, atom or bond center, descriptor, permutation,
and ordered controlling atoms.

``stereoisomers()`` returns a lazy iterator. The default options enumerate only
unassigned atom and double-bond stereo, include enhanced stereo-group flippers,
deduplicate by canonical isomeric SMILES, and yield at most 1,024 outputs. A
molecule with no selected center yields one isolated molecule value.

.. code-block:: python

   from cosmolkit import Molecule, StereoisomerOptions

   source = Molecule.from_smiles("CC(F)C(Cl)Br")
   analysis = source.analyze_potential_stereo()

   print([(item.center_kind, item.center_index) for item in analysis.stereo_info])
   print(source.stereoisomer_count())

   options = StereoisomerOptions(max_isomers=4, rand=0xF00D)
   for isomer in source.stereoisomers(options):
       print(isomer.to_smiles())

   assert source.to_smiles() == "CC(F)C(Cl)Br"

``max_isomers`` bounds successful outputs. When it is smaller than the
configuration space, ``rand`` follows Python ``random.Random`` semantics; a
``random.Random`` instance or subclass can supply ``getrandbits()`` lazily.
``try_embedding=True`` applies the source-defined one-conformer geometry
filter, and embedding failures do not consume the successful-output limit.
``only_unassigned=False`` also enumerates assigned centers,
``only_stereo_groups=True`` restricts selection to enhanced stereo groups, and
``unique=False`` retains source configurations that canonicalize to the same
isomeric SMILES.

Iterator errors are raised when the corresponding output is requested. The
source molecule remains unchanged, including when candidate discovery,
custom random generation, finalization, or embedding fails. This API matches
the pinned RDKit Python ``EnumerateStereoisomers`` boundary; the separate newer
C++ enumerator and atropisomer enumeration are not part of this surface.

Atom Chiral Tags From 3D Coordinates
------------------------------------

``with_chiral_tags_from_structure()`` assigns atom chiral tags from one stored
3D conformer and returns a new molecule. ``conf_id=-1`` selects the first
conformer; pass a nonnegative conformer id to select a specific conformer.
Existing tags are replaced by default and can be preserved with
``replace_existing_tags=False``. Coordinates, atom and bond ordering, and
unrelated properties are preserved.

.. code-block:: python

   import numpy as np
   from cosmolkit import Molecule

   mol = Molecule.from_smiles("C(F)(Cl)Br").with_only_3d_conformer(
       np.array(
           [
               [0.0, 0.0, 0.0],
               [1.0, 0.0, 0.0],
               [0.0, 1.0, 0.0],
               [0.0, 0.0, 1.0],
           ]
       )
   )
   assigned = mol.with_chiral_tags_from_structure()

   print(assigned.atoms()[0].chiral_tag())

The in-place form is ``assign_chiral_tags_from_structure_()``. Missing
conformers and invalid source state raise structured Python exceptions without
committing partial molecule changes. A selected non-3D conformer is a
source-defined no-op. This stable ``supported_with_rdkit_parity`` capability
matches all 77 fixed full-state oracle records exactly against RDKit 2026.03.1
``assignChiralTypesFrom3D``. Its boundary includes tetrahedral C/S/Se centers,
environment-enabled square-planar, trigonal-bipyramidal, and octahedral
centers, property updates, no-op paths, and defined errors. It does not include the
broader ``assignStereochemistryFrom3D`` workflow, 3D double-bond direction or
E/Z assignment, CIP orchestration, or distinct-substituent validation.

Modern CIP Labels
-----------------

``with_cip_labels()`` assigns molecular-context CIP descriptors and returns a
new molecule. ``assign_cip_labels_()`` is the explicit in-place form. After
assignment, atom and bond records expose ``cip_descriptor()`` and
``cip_neighbor_order()``; ``cip_computed()`` reports the molecule-level
completion state.

.. code-block:: python

   from cosmolkit import Molecule

   mol = Molecule.from_smiles("C[C@H](F)Cl")
   labeled = mol.with_cip_labels()
   print(labeled.atoms()[1].cip_descriptor())

   alkene = Molecule.from_smiles("F/C=C/F")
   alkene.assign_cip_labels_()
   print(alkene.bonds()[1].cip_descriptor())

Pass ``atoms=[...]`` or ``bonds=[...]`` for selected assignment. Priorities
remain molecular-context dependent; selection controls which configurations
are labeled, not how their priorities are calculated. The source-backed
boundary includes tetrahedral ``R/S/r/s``, double-bond ``E/Z``, and
atropisomeric ``M/P/m/p`` descriptors constructed by the pinned RDKit
``findConfigs`` dispatcher. The maintained ChEMBL 37 phase compares full,
selected-atom, selected-bond, and empty-selection state across 2,854,362
eligible molecules with exact agreement.

Conformer Generation And Force-Field Optimization
-------------------------------------------------

Conformer generation APIs create native 3D conformers through the source-ported
distance-geometry path. The default value-style operation uses ETKDGv3 and
returns a new molecule value.

.. code-block:: python

   from cosmolkit import EmbedParameters, Molecule

   mol = Molecule.from_smiles("CC(=O)NC").with_hydrogens()

   params = EmbedParameters.etkdg_v3()
   params.random_seed = 0xF00D
   params.num_threads = 1
   params.track_failures = True

   embedded = mol.with_3d_conformer(params)

   print(embedded.num_conformers())
   print(embedded.coordinates_3d())
   print(params.failures)

Manual Coordinate Assignment
----------------------------

Coordinate setters take a complete coordinate block and validate it before the
molecule value is updated. This avoids exposing partially edited conformer
state through Python.

.. code-block:: python

   import numpy as np
   from cosmolkit import Molecule

   mol = Molecule.from_smiles("CCO")

   coords_2d = np.array(
       [
           [0.0, 0.0],
           [1.5, 0.0],
           [2.1, 1.2],
       ]
   )
   drawn = mol.with_2d_coordinates(coords_2d)

   coords_3d = np.array(
       [
           [0.0, 0.0, 0.0],
           [1.5, 0.0, 0.0],
           [2.1, 1.2, 0.4],
       ]
   )
   placed = mol.with_added_3d_conformer(coords_3d)

   shifted = placed.with_3d_coordinates(coords_3d + [0.0, 0.0, 1.0])
   single = placed.with_only_3d_conformer(coords_3d + [0.0, 0.0, 2.0])
   cleared = placed.with_cleared_3d_conformers()

   print(drawn.coordinates_2d())
   print(placed.num_conformers())
   print(shifted.coordinates_3d())
   print(single.num_conformers())
   print(cleared.num_conformers())

The in-place forms follow COSMolKit's trailing-underscore convention:

.. code-block:: python

   mol.set_2d_coordinates_(coords_2d)
   conf_id = mol.add_3d_conformer_(coords_3d)
   mol.set_3d_coordinates_(coords_3d + [0.0, 0.0, 1.0], conformer_index=conf_id)
   mol.clear_3d_conformers_()
   conf_id = mol.set_only_3d_conformer_(coords_3d)

2D assignment accepts shape ``(num_atoms, 2)`` or ``(num_atoms, 3)``. For
three-column input, ``z_policy`` controls whether the z column is ignored,
required to be zero, or rejected:

.. code-block:: python

   coords_2d_with_zero_z = np.column_stack([coords_2d, np.zeros(mol.num_atoms())])
   mol.with_2d_coordinates(coords_2d_with_zero_z, z_policy="require_zero")

3D assignment accepts only shape ``(num_atoms, 3)``. All coordinate values must
be finite, and row counts must match ``mol.num_atoms()``.
``with_only_3d_conformer(coords)`` is the direct value-style equivalent of
RDKit ``RemoveAllConformers(); AddConformer(conf, assignId=True)`` for manual
coordinate assignment; the in-place form is ``set_only_3d_conformer_(coords)``
and returns conformer id ``0``.

For multi-conformer generation, explicit seeds are deterministic. RMS pruning,
sequential seed expansion, and terminal-group symmetrization for pruning follow
the source-ported RDKit path.

.. code-block:: python

   params = EmbedParameters.etkdg()
   params.random_seed = 123
   params.num_threads = 1
   params.prune_rms_thresh = 0.5
   params.enable_sequential_random_seeds = True

   pruned = mol.with_3d_conformers(5, params)
   print(pruned.num_conformers())

Molecular Alignment And RMSD
----------------------------

MolAlign measurement methods do not change either molecule. Use
``alignment_transform_to()`` for an explicit or automatic first-match
transform, ``best_alignment_to()`` or ``best_rmsd_to()`` for symmetry-aware
best-map alignment, and ``coordinate_rmsd_to()`` to measure coordinates in
their existing frames without alignment.

.. code-block:: python

   import numpy as np
   from cosmolkit import AlignmentAtomMap, AlignmentParameters, Molecule

   reference_coords = np.array(
       [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]]
   )
   probe_coords = reference_coords + np.array([3.0, -2.0, 1.0])
   reference = Molecule.from_smiles("CCC").with_only_3d_conformer(reference_coords)
   probe = Molecule.from_smiles("CCC").with_only_3d_conformer(probe_coords)
   params = AlignmentParameters(
       atom_map=[AlignmentAtomMap(index, index) for index in range(3)]
   )

   measured = probe.alignment_transform_to(reference, params)
   aligned, applied = probe.with_alignment_to(reference, params)

   assert np.array_equal(probe.coordinates_3d(), probe_coords)
   assert measured.rmsd() < 1.0e-8
   assert applied.rmsd() < 1.0e-8
   assert np.allclose(aligned.coordinates_3d(), reference_coords)

``with_alignment_to()`` returns a new molecule and result. ``align_to_()`` is
the explicit in-place form. ``all_conformer_best_rmsds()`` is read-only and
returns RDKit's triangular pair order; ``with_aligned_conformers()`` and
``align_conformers_()`` expose value-style and in-place conformer-set
alignment. Its ``AllConformerRmsdParameters`` exposes exactly the source call's
map, weight, match-limit, terminal-symmetry, hydrogen, and threading controls;
reflection and iteration count are fixed by that RDKit path. Weighted and
reflected pair alignment, stored conformer IDs, match limits, terminal-group
symmetry, and iteration limits are available through the other typed parameter
classes. O3A and MMFF/Crippen alignment scoring are separate capabilities
outside this ordinary MolAlign boundary.

UFF and MMFF optimization APIs operate on existing or generated 3D conformers
and return new molecule values through result objects. They do not mutate the
source molecule.

.. code-block:: python

   from cosmolkit import Molecule

   mol = Molecule.from_smiles("CCO").with_hydrogens().with_3d_conformer()

   if mol.has_uff_params():
       result = mol.with_uff_optimized(max_iters=200)
       optimized = result.molecule()

       print(not result.needs_more())
       print(result.status_code())
       print(result.energy())
       print(optimized.coordinates_3d())

   if mol.has_mmff_params():
       result = mol.with_mmff_optimized(mmff_variant="MMFF94", max_iters=200)
       optimized = result.molecule()

       print(not result.needs_more())
       print(result.status_code())

Substructure And SMARTS
-----------------------

Substructure matching functions accept ordinary molecules or query-bearing
``Molecule`` values compiled by the canonical SMARTS parser:

.. code-block:: python

   import cosmolkit

   mol = Molecule.from_smiles("CCO")
   query = Molecule.from_smiles("CO")

   print(cosmolkit.has_substruct_match(mol, query))
   print(cosmolkit.get_substruct_match(mol, query).atom_mapping())

``parse_smarts()`` returns an ordinary ``Molecule`` carrying the compiled query
graph. The same value can be passed directly to the substructure functions.

.. code-block:: python

   smarts = cosmolkit.parse_smarts("[#6]-O")

   print(smarts.num_atoms())
   print(smarts.num_bonds())
   print(smarts.to_smarts())
   print(cosmolkit.has_substruct_match(mol, smarts))

Query-bearing molecules can be written as ordinary SMARTS or CXSMARTS without
converting them to a separate query type:

.. code-block:: python

   labeled = cosmolkit.parse_smarts("[#6] |$site$|")

   print(labeled.to_smarts())
   print(labeled.to_cx_smarts())
