Installation
============

.. meta::
   :description: Install the COSMolKit Rust-native cheminformatics toolkit for Python from PyPI and verify the package version and NumPy coordinate support.

Install COSMolKit from PyPI:

.. code-block:: bash

   pip install cosmolkit

Verify the installation:

.. code-block:: python

   import cosmolkit

   print(cosmolkit.__version__)

COSMolKit depends on NumPy for array outputs such as coordinates and distance
bounds matrices.
