COSMolKit Python
================

COSMolKit is a Python package for molecule graph workflows, SMILES/SDF/MOL2/XYZ
IO, coordinate access, native 3D conformer generation, UFF/MMFF optimization,
Morgan fingerprints, molecule depiction, Python pickle round-tripping, SMARTS
parse metadata, protein structure access, and high-throughput batch processing.

Important API model: COSMolKit presents value-style ``Molecule`` objects.
Transform methods return new ``Molecule`` objects and do not mutate the
original object in place. Internally the library uses copy-on-write (COW)
storage to share unchanged data efficiently. This differs from common RDKit
Python workflows where code often updates an existing molecule or ``RWMol``
directly.

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   quickstart
   confseq
   molecule
   batch
   fingerprints
   protein
   io
   api
