COSMolKit — Rust-native cheminformatics for Python
==================================================

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

Project and Rust crates
-----------------------

COSMolKit is developed in the
`COSMolKit source repository <https://github.com/cosmol-studio/COSMolKit>`_.

- `cosmolkit <https://crates.io/crates/cosmolkit>`_ —
  Rust facade crate and
  `Rust API documentation <https://docs.rs/cosmolkit/latest/cosmolkit/>`_.
- `cosmolkit-core <https://crates.io/crates/cosmolkit-core>`_ —
  core molecular graph, chemistry, IO, and operations.
- `cosmolkit-inchi <https://crates.io/crates/cosmolkit-inchi>`_ —
  pure-Rust InChI implementation.
- `cosmolkit-ringdecomposer <https://crates.io/crates/cosmolkit-ringdecomposer>`_ —
  ring decomposition support.
- `COSMolKit Web Tools <https://github.com/cosmol-studio/cosmolkit-tools-web>`_ —
  browser-native cheminformatics tools powered by WebAssembly.
