JavaScript / TypeScript
=======================

.. meta::
   :description: Choose a COSMolKit WebAssembly feature package and release channel, then parse, search, draw, fingerprint, optimize, and react molecules from JavaScript or TypeScript.

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

**Feature combination and release version are separate choices.** Each
combination has its own package name. ``latest`` selects its stable release;
``rc`` selects its release candidate. An exact version or a committed lockfile
makes the choice reproducible.

.. code-block:: console

   # Stable analysis package
   npm install @cosmol-studio/cosmolkit-core-analysis

   # Release-candidate analysis package
   npm install @cosmol-studio/cosmolkit-core-analysis@rc

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

In Node.js, supply the WebAssembly bytes instead of relying on browser fetching:

.. code-block:: javascript

   import { readFileSync } from "node:fs";
   import init, * as ck from "@cosmol-studio/cosmolkit";

   const wasmUrl = import.meta.resolve(
     "@cosmol-studio/cosmolkit/wasm_wasm_bg.wasm",
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
     console.log(mol.numAtoms(), mol.numBonds()); // 3, 2
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
receiver unchanged. ``addHydrogens()`` modifies the receiver. Both use the same
Rust chemistry; the value-returning form uses block-level copy-on-write.

.. code-block:: javascript

   const mol = ck.Molecule.fromSmiles("CCO");
   const hydrogenated = mol.withHydrogens();
   const tagged = mol.withAtomProperty(0, "tracking_id", 42);
   try {
     console.log(mol.numAtoms(), hydrogenated.numAtoms()); // 3, 9
     console.log(tagged.atomProperty(0, "tracking_id")); // 42, a number
     console.log(mol.atomProperty(0, "tracking_id")); // null
     mol.addHydrogens();
     console.log(mol.numAtoms()); // 9
   } finally {
     tagged.free();
     hydrogenated.free();
     mol.free();
   }

Read and write molecular text
-----------------------------

MOL/SDF parsing and writing belong to ``core``. Automatic generation of missing
2D coordinates additionally requires ``depict``; this example uses ``full``.
For an existing coordinate-bearing MOL/SDF, parse the string with
``Molecule.fromMol(text)`` or ``Molecule.fromSdf(text)``.

.. code-block:: javascript

   const mol = ck.Molecule.fromSmiles("CCO");
   const positioned = mol.with2dCoordinates();
   try {
     const sdf = positioned.toSdf(); // String, not a filesystem write
     const restored = ck.Molecule.fromSdf(sdf);
     try {
       console.log(restored.toSmiles()); // CCO
     } finally {
       restored.free();
     }
   } finally {
     positioned.free();
     mol.free();
   }

In a browser, obtain input with ``await file.text()`` or
``await (await fetch(url)).text()``, and pass the resulting string to a parser. Download output
using a ``Blob`` and your application's download UI. Native filesystem methods
are not a browser file picker. Binary archives (``toBinary``/``fromBinary``) are
not compiled into WebAssembly, including ``full``.

SMARTS and substructure matching
--------------------------------

Requires ``search``. Pass a plain options object, or use an explicit
``SubstructMatchParams`` instance. Results map query atoms to molecule atoms.

.. code-block:: javascript

   const mol = ck.Molecule.fromSmiles("CCO");
   const query = ck.QueryGraph.fromSmarts("[#6]-[#8]");
   try {
     console.log(mol.hasSubstructMatch(query)); // true
     const matches = mol.substructMatches(query, {
       useChirality: true,
       uniquify: true,
       maxMatches: 100,
     });
     try {
       console.log(matches.map((match) => match.atomMapping())); // [[1, 2]]
     } finally {
       for (const match of matches) match.free();
     }
   } finally {
     query.free();
     mol.free();
   }

Fingerprints and similarity
---------------------------

Requires ``fingerprints``. A ``Fingerprint`` is an owned object, not a dense
array. ``onBits()`` returns the set-bit indices; ``nBits()`` is its full length.

.. code-block:: javascript

   const ethanol = ck.Molecule.fromSmiles("CCO");
   const propanol = ck.Molecule.fromSmiles("CCCO");
   const left = ethanol.fingerprintMorgan();
   const right = propanol.fingerprintMorgan();
   try {
     console.log(left.tanimoto(right));
     const dense = new Uint8Array(left.nBits());
     for (const bit of left.onBits()) dense[bit] = 1;
     console.log(dense); // One 0/1 entry per fingerprint bit
   } finally {
     right.free();
     left.free();
     propanol.free();
     ethanol.free();
   }

Descriptors
-----------

Requires ``descriptors``; choose ``core-analysis`` for descriptors, search and
fingerprints together.

.. code-block:: javascript

   const mol = ck.Molecule.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
   const crippen = mol.crippenDescriptors();
   try {
     console.table({
       formula: mol.molecularFormula(),
       molecularWeight: mol.molecularWeight(),
       exactMass: mol.exactMolecularWeight(),
       logP: crippen.logp,
       tpsa: mol.tpsa(),
       donors: mol.numHbd(),
       acceptors: mol.numHba(),
     });
   } finally {
     crippen.free();
     mol.free();
   }

Draw SVG and retain 2D coordinates
----------------------------------

Requires ``depict``. ``toSvg()`` returns markup without storing drawing
coordinates back into the original molecule. Use ``with2dCoordinates()`` if
you also need a molecule containing those coordinates.

.. code-block:: javascript

   const mol = ck.Molecule.fromSmiles("c1ccccc1O");
   const positioned = mol.with2dCoordinates();
   try {
     const svg = positioned.toSvg(320, 240);
     console.log(svg); // Insert into your application's molecule container
     console.log(positioned.coordinates2d()); // Flattened x, y pairs
     console.log(mol.coordinates2d().length); // 0: original unchanged
   } finally {
     positioned.free();
     mol.free();
   }

3D conformers, forcefields and RMSD
-----------------------------------

Requires ``conformer``; install ``core-3d`` or ``full``. Add explicit hydrogens,
generate a conformer, then optimize it. MMFF requires its own parameterization;
do not silently substitute UFF when MMFF parameters are unavailable.

.. code-block:: javascript

   const mol = ck.Molecule.fromSmiles("CCO");
   const hydrogenated = mol.withHydrogens();
   const defaults = ck.EmbedParams.etkdgV3();
   const params = defaults.withJson('{"randomSeed":42}');
   const conformer = hydrogenated.with3dConformerWithParams(params);
   try {
     console.log(conformer.numConformers3d()); // 1
     console.log(conformer.coordinates3d(0)); // Flattened x, y, z triples
     const uff = conformer.withUffOptimized();
     const optimized = uff.molecule();
     try {
       console.log(uff.energy(), uff.needsMore());
       console.log(optimized.bestRmsdTo(conformer));
     } finally {
       optimized.free();
       uff.free();
     }
     if (conformer.mmffHasAllMoleculeParams()) {
       const mmff = conformer.withMmffOptimized();
       const mmffMol = mmff.molecule();
       try {
         console.log(mmff.statusCode(), mmff.needsMore());
       } finally {
         mmffMol.free();
         mmff.free();
       }
     }
   } finally {
     conformer.free();
     params.free();
     defaults.free();
     hydrogenated.free();
     mol.free();
   }

``needsMore()`` reports that the optimizer has not converged within its limit;
it is not a reason to label the returned coordinates converged. For multiple
conformers use ``with3dConformersWithParams(count, params)`` and the corresponding
``withUffOptimizedConformers()`` or ``withMmffOptimizedConformers()`` operation.

Run a reaction
--------------

Requires ``reaction``, which also includes ``search``. Supply reactants in the
template order. The outer array contains alternative product sets; each inner
array contains the products from one execution. An empty outer array means no
products, not a parse error.

.. code-block:: javascript

   const rxn = ck.Reaction.fromSmirks("[C:1][O:2]>>[C:1]=[O:2]");
   const reactant = ck.Molecule.fromSmiles("CCO");
   const params = new ck.ReactionRunParams();
   params.maxProducts = 100;
   params.copyAtomProperties = true;
   try {
     const products = rxn.run([reactant], params);
     try {
       console.log(products.map((set) => set.map((mol) => mol.toSmiles())));
       console.log(reactant.toSmiles()); // CCO: input unchanged
     } finally {
       for (const set of products) for (const product of set) product.free();
     }
   } finally {
     params.free();
     reactant.free();
     rxn.free();
   }

Process a batch
---------------

Batch construction and SMILES output require ``core``; the fingerprint step
also requires ``fingerprints``. Results retain input order. APIs returning
nullable entries use ``null`` for failed records; inspect ``errors()`` and the
valid/invalid counts rather than silently dropping those positions.

.. code-block:: javascript

   const batch = ck.MoleculeBatch.fromSmilesList(["CCO", "CCCO", "c1ccccc1"]);
   try {
     console.log(batch.toSmilesList());
     console.log(batch.validCount(), batch.invalidCount()); // 3, 0
     const fingerprints = batch.fingerprintMorganList();
     try {
       console.log(fingerprints.map((fp) => fp === null ? null : fp.onBits()));
     } finally {
       for (const fp of fingerprints) fp?.free();
     }
   } finally {
     batch.free();
   }

Fragments and scaffolds
-----------------------

Requires ``core``. Fragment extraction returns separate owned molecules;
``largestFragment()`` selects one, while ``murckoScaffold()`` extracts a scaffold.

.. code-block:: javascript

   const mol = ck.Molecule.fromSmiles("CCc1ccccc1.[Na+]");
   const fragments = mol.fragments();
   const largest = mol.largestFragment();
   const scaffold = largest.murckoScaffold();
   try {
     console.log(fragments.map((fragment) => fragment.toSmiles()));
     console.log(largest.toSmiles(), scaffold.toSmiles());
   } finally {
     scaffold.free();
     largest.free();
     for (const fragment of fragments) fragment.free();
     mol.free();
   }

InChI conversion
----------------

Requires ``inchi``; install ``core-inchi`` or ``full``.

.. code-block:: javascript

   const mol = ck.Molecule.fromSmiles("CCO");
   try {
     const inchi = mol.toInchi();
     console.log(inchi, mol.toInchiKey());
     const restored = ck.Molecule.fromInchi(inchi);
     try {
       console.log(restored.toSmiles()); // CCO
     } finally {
       restored.free();
     }
   } finally {
     mol.free();
   }

Biological structures and selections
------------------------------------

Requires ``bio``; install ``core-bio`` or ``full``. Read PDB/mmCIF text and select
using a Gemmi-style CID (model / chain / residue / atom). This small PDB example
contains one residue in chain A:

.. code-block:: javascript

   const pdb = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \nEND\n";
   const structure = ck.BioStructure.fromPdb(pdb);
   const selection = ck.BioSelection.fromCid("/1/A/1");
   const selected = structure.withSelection(selection);
   try {
     console.log(structure.numModels(), structure.numChains(), structure.numResidues());
     console.log(Array.from(structure.selectedAtomIds(selection))); // [0]
     console.log(selected.toPdb()); // New structure; original unchanged
   } finally {
     selected.free();
     selection.free();
     structure.free();
   }

For mmCIF use ``BioStructure.fromMmcif(text)`` and ``structure.toMmcif()``.

Ownership and API reference
---------------------------

Release owned WebAssembly objects with ``free()`` when finished, including
objects returned inside arrays and result wrappers. Do not use a handle after
freeing it. Use ``try/finally`` around long-lived application work; ordinary
strings, numbers and detached JavaScript arrays do not need ``free()``. The
short examples above explicitly release their successfully created handles.

JavaScript and TypeScript use the same runtime API; TypeScript additionally
checks the package's generated declarations. Methods use camelCase. Consult the
reference for exact overloads, configuration types and which operations mutate
their receiver; do not infer that a Python signature is a JavaScript signature.

The reference below is generated from the full distribution's actual ``.d.ts``
entry point using TypeDoc and sphinx-js. Smaller distributions expose only the
APIs compiled into their combination.

.. toctree::
   :maxdepth: 1

   javascript-api
