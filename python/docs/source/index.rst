COSMolKit Python
================

COSMolKit is a Python package for molecule graph workflows, SMILES/SDF IO,
coordinate access, molecule depiction, and high-throughput batch processing.

Important API model: COSMolKit uses copy-on-write (COW) molecule values.
Transform methods return new ``Molecule`` objects and do not mutate the
original object in place. This differs from common RDKit Python workflows where
code often updates an existing molecule or ``RWMol`` directly.

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   quickstart
   molecule
   batch
   io
   api
