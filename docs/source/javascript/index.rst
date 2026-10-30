JavaScript / TypeScript
=======================

.. meta::
   :description: Choose a COSMolKit WebAssembly feature package, then parse, search, draw, fingerprint, optimize, and react molecules from JavaScript or TypeScript.

Choose a package and version
----------------------------

Start with the **full package** if you want all available WebAssembly APIs:

.. code-block:: console

   npm install @cosmol-studio/cosmolkit

For a smaller application, choose one of these prebuilt combinations. All names
below are under the ``@cosmol-studio/`` scope; install and import the same name.

.. list-table:: Prebuilt feature combinations
   :header-rows: 1
   :widths: 29 32 39

   * - Package
     - Public features
     - Use it for
   * - ``cosmolkit`` (default)
     - ``full``
     - All operations in this guide
   * - ``cosmolkit-core``
     - ``core``
     - SMILES, molecular IO, hydrogen handling, fragments, batch processing
   * - ``cosmolkit-core-search``
     - ``core + search``
     - SMARTS and substructure matching
   * - ``cosmolkit-core-fingerprints``
     - ``core + fingerprints``
     - Fingerprints and similarity
   * - ``cosmolkit-core-analysis``
     - ``core + search + fingerprints + descriptors``
     - Search, similarity and molecular descriptors together
   * - ``cosmolkit-core-reaction``
     - ``core + reaction + search``
     - SMIRKS reactions and their SMARTS matching
   * - ``cosmolkit-core-depict``
     - ``core + depict``
     - SVG drawing and 2D coordinate generation
   * - ``cosmolkit-core-3d``
     - ``core + conformer``
     - 3D conformers, forcefields and alignment
   * - ``cosmolkit-core-bio``
     - ``core + bio``
     - PDB/mmCIF structures and selections
   * - ``cosmolkit-core-inchi``
     - ``core + inchi``
     - InChI and InChIKey conversion

``reaction`` includes ``search``; ``conformer`` includes forcefields and
alignment. Other domains are not implicitly enabled: for example,
``core-analysis`` does not include drawing. Use ``full`` when your application
needs a combination not listed above. Installing multiple smaller packages
does not merge their compiled APIs into one module.

Each feature combination has its own package name. Installation selects the
stable release; an exact version or a committed lockfile makes the choice
reproducible.

.. code-block:: console

   npm install @cosmol-studio/cosmolkit-core-analysis

For an exact release, append ``@<published-version>`` and use
``npm install --save-exact``. Feature names are **not** version tags: use
``@cosmol-studio/cosmolkit-core-bio``, not ``@cosmol-studio/cosmolkit@core-bio``.
An upgrade keeps the selected combination because its package name stays the
same. The documentation's version menu selects documentation, not your installed
npm version; check that version's generated TypeScript declarations for its API.

Initialize once
---------------

The package is an ES module with TypeScript declarations included. In a browser
application with a bundler that serves WebAssembly assets, initialize before
creating objects. All following examples use the full package:

.. code-block:: javascript

   import init, * as ck from "@cosmol-studio/cosmolkit";

   await init();

For a smaller package, change that import's module name, for example to
``"@cosmol-studio/cosmolkit-core-analysis"``. Only its compiled capabilities
appear in its exports and declarations. Keep the generated JavaScript and
``.wasm`` asset from the same package version together.

Top-level constructors such as ``ck.molFromSmiles(text)`` and
``ck.bioFromMmcif(text)`` call the corresponding class factories
(``ck.Molecule.fromSmiles`` and ``ck.BioStructure.fromMmcif``). They accept the
same options, return the same types and preserve the same errors. Their
availability follows the selected package's features.

In Node.js, supply the WebAssembly bytes instead of relying on browser fetching:

.. code-block:: javascript

   import { readFileSync } from "node:fs";
   import init, * as ck from "@cosmol-studio/cosmolkit";

   const wasmUrl = import.meta.resolve(
     "@cosmol-studio/cosmolkit/cosmolkit_wasm_bg.wasm",
   );
   await init({ module_or_path: readFileSync(new URL(wasmUrl)) });

Parse and inspect a molecule
----------------------------

Requires ``core``. SMILES parsing sanitizes by default; invalid input throws
rather than returning a partially valid molecule. Atom indices are zero-based.

.. code-block:: javascript

   const mol = ck.Molecule.fromSmiles("OCC");
   try {
     console.log(mol.toSmiles()); // CCO
     console.log(mol.numAtoms, mol.numBonds); // 3, 2
   } finally {
     mol.free();
   }

.. code-block:: javascript

   try {
     const invalid = ck.Molecule.fromSmiles("C==C");
     invalid.free();
   } catch (error) {
     console.error(error); // Parsing/chemistry error, not a successful result
   }

New values, hydrogens and atom properties
-----------------------------------------

Requires ``core``. ``withHydrogens()`` returns a new molecule and leaves its
receiver unchanged. ``addHydrogens_()`` modifies the receiver. Both use the same
Rust chemistry; the value-returning form uses block-level copy-on-write.

.. code-block:: javascript

   const mol = ck.Molecule.fromSmiles("CCO");
   const hydrogenated = mol.withHydrogens();
   const tagged = mol.withAtomPropertyInt(0, "tracking_id", 42);
   try {
     console.log(mol.numAtoms, hydrogenated.numAtoms); // 3, 9
     console.log(tagged.atomPropertyInt(0, "tracking_id")); // 42, a number
     console.log(mol.atomPropertyInt(0, "tracking_id")); // undefined
     mol.addHydrogens_();
     console.log(mol.numAtoms); // 9
   } finally {
     tagged.free();
     hydrogenated.free();
     mol.free();
   }

Fingerprints and arrays
-----------------------

.. code-block:: javascript

   const mol = ck.molFromSmiles("CCO");
   const fingerprint = mol.fingerprintMorgan({ radius: 2, fpSize: 2048 });
   try {
     console.log(fingerprint.nBits, fingerprint.onBits);
   } finally {
     fingerprint.free();
     mol.free();
   }

Release owned WASM objects with ``free()`` when finished. Ordinary JavaScript
strings, numbers and arrays do not need it. Keep objects alive until their last
use; see :doc:`workflows` for complete IO, search, drawing, 3D, batch and BIO
examples with resource handling.

.. toctree::
   :maxdepth: 1

   workflows
   api
