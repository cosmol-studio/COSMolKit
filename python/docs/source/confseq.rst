ConfSeq
=======

.. meta::
   :description: Decode ConfSeq molecular conformations with COSMolKit using source-aligned distance geometry or the optional fast geometry coordinate backend.

ConfSeq decoding uses an input SMILES string and a ConfSeq string to recover a
3D molecule. COSMolKit exposes two template backends:

``distance_geometry``
   Reference backend. It uses the RDKit-derived distance-geometry path and is
   the default.

``fast_geometry``
   Lightweight backend. It only replaces the initial coordinate source; after
   that, it runs the same ConfSeq token decode as ``distance_geometry``. It is
   opt-in and may raise ``ValueError`` for unsupported FastGeometry cases.

Decode one tokenized corpus record with the reference backend. The spaces are
part of the ConfSeq token stream, not cosmetic formatting: tokens such as
``<112> |`` and bracketed atoms must remain separated. The ConfSeq string
contains both the molecular topology and recovered numeric angle/dihedral
tokens, so normal callers pass only this one string.

.. code-block:: python

   import cosmolkit

   confseq = (
       "N C <112> | ( = O ) <84> C <111> | <-45> n 1 c c ( <21> N <123> | "
       "<173> C <116> | ( = O ) <120> c 2 c c c c c 2 <112> C <112> | <-66> "
       "N 2 <-172> C ( = O ) <-2> C <0> N <3> C <174> 2 = O ) c n 1"
   )

   mol = cosmolkit.confseq.decode(
       confseq,
       optimize_with_uff=False,
       template_backend="distance_geometry",
   )

   print(mol.num_atoms())
   print(mol.coordinates_3d().shape)

Decode the same record with FastGeometry:

.. code-block:: python

   fast = cosmolkit.confseq.decode(
       confseq,
       optimize_with_uff=False,
       template_backend="fast_geometry",
   )

   print(fast.num_conformers())

Batch decoding accepts the same backend option:

.. code-block:: python

   confseq_1 = (
       "N C <113> | ( = O ) <-54> C <113> | <44> n 1 c c ( <-132> N <128> | "
       "<-9> C <115> | ( = O ) <-118> c 2 c c c c c 2 <-95> C <116> | <117> "
       "N 2 <176> C ( = O ) <6> C <-8> N <11> C <-176> 2 = O ) c n 1"
   )

   mols = cosmolkit.confseq.decode_batch(
       [confseq, confseq_1],
       errors="keep",
       n_jobs=8,
       optimize_with_uff=False,
       template_backend="distance_geometry",
   )

   print([mol is not None for mol in mols])

Use ``template_backend="distance_geometry"`` when measuring ConfSeq reference
behavior or comparing FastGeometry output. Use
``template_backend="fast_geometry"`` only when testing the fast
initialization path explicitly.
