COSMolKit — Rust-native cheminformatics for Python
==================================================

.. meta::
   :description: COSMolKit is a Rust-native cheminformatics and structural biology toolkit with Python bindings for SMILES, SMARTS, fingerprints, conformers, UFF/MMFF, InChI, molecular IO, and batch workflows.

COSMolKit provides first-class Python bindings to a Rust-native cheminformatics and structural biology toolkit. The Python API covers molecular graph operations, SMILES/SMARTS and molecular file workflows, coordinate access, native 3D conformer generation, UFF/MMFF optimization, InChI, fingerprints, molecular descriptors, molecule depiction, protein structures, and high-throughput batch processing.

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
   descriptors
   protein
   io
   api
