JavaScript Workflows
====================

The examples below use the full package. Initialize it once before running
any workflow; see :doc:`index` for installation and feature combinations.

.. code-block:: javascript

   import init, * as ck from "@cosmol-studio/cosmolkit";
   await init();

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
       console.log(matches.map((match) => Array.from(match.atomMapping))); // [[1, 2]]
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
array. ``onBits`` returns the set-bit indices; ``nBits`` is its full length.

.. code-block:: javascript

   const ethanol = ck.Molecule.fromSmiles("CCO");
   const propanol = ck.Molecule.fromSmiles("CCCO");
   const left = ethanol.fingerprintMorgan();
   const right = propanol.fingerprintMorgan();
   try {
     console.log(left.tanimoto(right));
     const dense = new Uint8Array(left.nBits);
     for (const bit of left.onBits) dense[bit] = 1;
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
     console.log(positioned.coordinates2d); // Readonly [x, y] pairs
     console.log(mol.coordinates2d); // undefined: original unchanged
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
   const params = ck.EmbedParams.etkdgV3();
   params.randomSeed = 42;
   const conformer = hydrogenated.with3dConformerWithParams(params);
   try {
     console.log(conformer.num3dConformers); // 1
     console.log(conformer.coordinates3d(0)); // Readonly [x, y, z] triples
     const uff = conformer.withUffOptimized();
     const optimized = uff.molecule;
     try {
       console.log(uff.energy, uff.needsMore);
       console.log(optimized.bestRmsdTo(conformer));
     } finally {
       optimized.free();
       uff.free();
     }
     if (conformer.mmffHasAllMoleculeParams()) {
       const mmff = conformer.withMmffOptimized();
       const mmffMol = mmff.molecule;
       try {
         console.log(mmff.statusCode, mmff.needsMore);
       } finally {
         mmffMol.free();
         mmff.free();
       }
     }
   } finally {
     conformer.free();
     params.free();
     hydrogenated.free();
     mol.free();
   }

``needsMore`` reports that the optimizer has not converged within its limit;
it is not a reason to label the returned coordinates converged. For multiple
conformers use ``with3dConformersWithParams(count, params)`` and the corresponding
``withUffOptimizedConformers()`` or ``withMmffOptimizedConformers()`` operation.

To collect embedding diagnostics, enable tracking and request a result:

.. code-block:: javascript

   const input = ck.molFromSmiles("CCO");
   try {
     const mol = input.withHydrogens();
     try {
       const params = ck.EmbedParams.etkdgV3();
       try {
         params.trackFailures = true;
         params.randomSeed = 42;
         const result = mol.with3dConformersResultWithParams(8, params);
         try {
           console.log(result.failureCounts());
           console.log(result.failureCount(ck.EmbedFailureCause.InitialCoords));
         } finally {
           result.free();
         }
       } finally {
         params.free();
       }
     } finally {
       mol.free();
     }
   } finally {
     input.free();
   }

``failureCounts()`` returns an array, or ``undefined`` when tracking is disabled.
Counters describe intermediate failure events, not missing conformers. The
input parameters remain unchanged. ``EndOfEnum`` is not a valid query cause.

Persistent force fields and constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Requires ``conformer``; use ``core-3d`` or ``full``. A persistent
``MolecularForceField`` owns its coordinates independently of the source
molecule. Create it with ``uffForceField()`` or
``mmffForceField({ mmffVariant: ck.MmffVariant.Mmff94 })``;
``ck.MmffVariant.Mmff94s`` selects the other MMFF model.
UFF, MMFF94 and MMFF94s support the same four soft constraints.

Methods accept either a parameter instance or a plain options object. All
constraint fields are required; the numerical values below are explicit example
choices, not implicit defaults. Parameter fields are writable.

.. code-block:: javascript

   const input = ck.molFromSmiles("CCCC");
   const hydrogenated = input.withHydrogens();
   const mol = hydrogenated.with3dConformer({ randomSeed: 42 });
   const field = mol.uffForceField();
   const anchor = new ck.ForceFieldPositionConstraintParams({maxDisplacement: 0.2, forceConstant: 100});
   let outcome;
   try {
     field.addDistanceConstraint_(0, 1, {
       relative: false, minDistance: 1.3, maxDistance: 1.6, forceConstant: 100,
     });
     field.addAngleConstraint_(0, 1, 2, {
       relative: true, minAngleDeg: -5, maxAngleDeg: 5, forceConstant: 1,
     });
     field.addTorsionConstraint_(0, 1, 2, 3, {
       relative: true, minDihedralDeg: -10, maxDihedralDeg: 10, forceConstant: 1,
     });
     anchor.forceConstant = 200;
     field.addPositionConstraint_(0, anchor);
     outcome = field.minimize_({ maxIterations: 20 });
     console.log(outcome.energy, outcome.converged);
     console.log(field.positions); // Owned coordinate snapshot, not a mutable view
   } finally {
     outcome?.free();
     anchor.free();
     field.free();
     mol.free();
     hydrogenated.free();
     input.free();
   }

The four methods are ``addDistanceConstraint_``, ``addAngleConstraint_``,
``addTorsionConstraint_`` and ``addPositionConstraint_``. Distances and
displacements use angstroms; angle and torsion bounds use degrees. Force
constants use the corresponding RDKit contribution's native convention,
without unit conversion. Constraints add energy penalties rather than clamping
coordinates; use ``setFixedAtoms_([atomIndex])`` to prevent an atom from moving.

``relative: true`` captures the current distance, angle or signed dihedral when
the constraint is added. Position anchors are captured at addition as well;
later ``setPosition_()``, ``setPositions_()`` and ``minimize_()`` calls do not
recapture either baseline. Torsion bounds preserve wrapped interval semantics.
Invalid indices, bounds or angle ranges throw typed ``ForceFieldError`` values
without appending a contribution. The source molecule remains unchanged;
explicitly use ``mol.with3dCoordinates(field.positions)`` to create a new
molecule with the resulting coordinates.

The :download:`complete JavaScript example <../../../wasm/examples/vanilla/forcefield_optimization.js>`
demonstrates all three models, fixed atoms and explicit coordinate installation.

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
optional entries use ``undefined`` for failed records; inspect ``errors`` and the
valid/invalid counts rather than silently dropping those positions.

.. code-block:: javascript

   const batch = ck.MoleculeBatch.fromSmilesList(["CCO", "CCCO", "c1ccccc1"]);
   try {
     console.log(batch.toSmiles());
     console.log(batch.validCount, batch.invalidCount); // 3, 0
     const fingerprints = batch.fingerprintMorgan();
     try {
       console.log(fingerprints.len, fingerprints.nBits);
       console.log(fingerprints.toUint8Array()); // Dense row-major bit matrix
     } finally {
       fingerprints.free();
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
     console.log(structure.numModels, structure.numChains, structure.numResidues);
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

See :doc:`api` for complete generated signatures.
