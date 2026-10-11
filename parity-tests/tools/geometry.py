"""Direct native geometry calls; no expected-coordinate fitting or rounding."""
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor, rdDistGeom, rdMolAlign, rdForceFieldHelpers as FF
from rdkit.Chem.Draw import rdMolDraw2D
from chemistry import state

TASKS = {"coordinates_2d", "svg", "bounds_matrix", "embedding", "alignment",
         "uff_available", "mmff_available", "uff_optimize", "mmff_optimize",
         "uff_optimize_conformers", "mmff_optimize_conformers", "uff_force_field", "mmff_force_field"}

def xyz(count, seed):
    rows = []
    for _ in range(count):
        row = []
        for _ in range(3):
            seed = (seed * 6364136223846793005 + 1442695040888963407) & ((1 << 64) - 1)
            row.append((((seed >> 32) % 257) - 128) / 16.0)
        rows.append(row)
    return rows

def conformers(mol):
    return [[c.GetId(), c.Is3D(), [[float(v) for v in c.GetAtomPosition(i)]
                                 for i in range(mol.GetNumAtoms())]] for c in mol.GetConformers()]

def with_coordinates(mol, seed, count):
    result = Chem.Mol(mol)
    result.RemoveAllConformers()
    for index in range(count):
        c = Chem.Conformer(result.GetNumAtoms())
        c.Set3D(True)
        for i, row in enumerate(xyz(result.GetNumAtoms(), seed + index)):
            c.SetAtomPosition(i, row)
        result.AddConformer(c, assignId=True)
    return result

def namespace(svg):
    # Same approved product identity projection as Rust. No numerical rewrite.
    return (svg.replace("xmlns:rdkit='http://www.rdkit.org/xml'", "xmlns:tool='__tool_namespace__'")
            .replace("xmlns:ck='https://kit.cosmol.org/'", "xmlns:tool='__tool_namespace__'")
            .replace("rdkit:", "tool:").replace("ck:", "tool:"))

def alignment(mol, p):
    probe, reference = with_coordinates(mol, p["seed"], 3), with_coordinates(mol, p["seed"] + 1, 1)
    original, original_reference = conformers(probe), conformers(reference)
    atom_map = [(i, i) for i in range(mol.GetNumAtoms())] if p["explicit_map"] else []
    weights = [1.0 + (i % 3) / 2.0 for i in range(mol.GetNumAtoms())] if p["weighted"] else []
    options = dict(weights=weights, reflect=p["reflect"], maxIters=p["max_iterations"])
    best = dict(map=[atom_map] if atom_map else [], maxMatches=p["max_matches"],
                symmetrizeConjugatedTerminalGroups=p["symmetrize"], weights=weights,
                numThreads=1)
    operation = p["operation"]
    after = original
    if operation in ("transform", "align"):
        rms, transform = rdMolAlign.GetAlignmentTransform(probe, reference, atomMap=atom_map, **options)
        # The source wrapper exposes no chosen map for this API. When omitted,
        # obtain the actual source match, not a guessed input-position map.
        mapping = atom_map or list(enumerate(probe.GetSubstructMatch(reference)))
        result = [rms, transform.tolist(), [list(pair) for pair in mapping]]
        if operation == "align":
            target = Chem.Mol(probe)
            rdMolAlign.AlignMol(target, reference, atomMap=atom_map, **options)
            after = conformers(target)
    elif operation == "best":
        rms, transform, mapping = rdMolAlign.GetBestAlignmentTransform(probe, reference, reflect=p["reflect"], maxIters=p["max_iterations"], **best)
        # CK's RMS-only query accepts the same reflection/iteration parameters;
        # the legacy GetBestRMS wrapper fixes those to false/50. Use the exact
        # native parameterized alignment result, without a source mutation.
        minimum = rdMolAlign.GetBestAlignmentTransform(probe, reference, reflect=p["reflect"], maxIters=p["max_iterations"], **best)[0]
        result = [[rms, transform.tolist(), [list(pair) for pair in mapping]], minimum]
    elif operation == "coordinate":
        result = rdMolAlign.CalcRMS(probe, reference, map=best["map"], maxMatches=p["max_matches"],
                                   symmetrizeConjugatedTerminalGroups=p["symmetrize"], weights=weights)
    elif operation == "all":
        values = rdMolAlign.GetAllConformerBestRMS(probe, **best)
        pairs = [(i, j) for i in range(3) for j in range(i)]
        result = [[i, j, value] for (i, j), value in zip(pairs, values)]
    elif operation == "conformers":
        target, result = Chem.Mol(probe), []
        rdMolAlign.AlignMolConformers(target, weights=weights, reflect=p["reflect"],
                                     maxIters=p["max_iterations"], RMSlist=result)
        after = conformers(target)
    else:
        raise ValueError("unknown alignment operation")
    return [result, after, original, original_reference]

def forcefield(mol, task, p):
    mol = Chem.AddHs(mol) if p["add_hs"] else Chem.Mol(mol)
    mmff = task.startswith("mmff_")
    available = FF.MMFFHasAllMoleculeParams(mol) if mmff else FF.UFFHasAllMoleculeParams(mol)
    if task.endswith("available"):
        return available
    if not available:
        return [False]
    mol = with_coordinates(mol, p["seed"], 2 if task.endswith("conformers") else 1)
    original = conformers(mol)
    def field(target):
        # AtomTyper.cpp::MMFFMolProperties mutates target (_MMFFSanitized,
        # Kekulize and setMMFFAromaticity). Prepare only for the requested
        # factory, never as an extra precondition of an optimizer call.
        props = FF.MMFFGetMoleculeProperties(target, mmffVariant=p["variant"]) if mmff else None
        f = (FF.MMFFGetMoleculeForceField(target, props, nonBondedThresh=p["threshold"], ignoreInterfragInteractions=p["ignore_interfragment"])
             if mmff else FF.UFFGetMoleculeForceField(target, vdwThresh=p["threshold"], ignoreInterfragInteractions=p["ignore_interfragment"]))
        if f is None:
            raise RuntimeError("available force-field factory returned None")
        f.Initialize()
        return f
    outcomes = []
    if task.endswith("force_field"):
        f = field(mol)
        fixed = [0] if p["fixed_first"] and mol.GetNumAtoms() else []
        for i in fixed:
            f.AddFixedPoint(i)
        def snapshot():
            gradient = f.CalcGrad()
            return [f.CalcEnergy(), [list(gradient[i:i+3]) for i in range(0, len(gradient), 3)], conformers(mol)[0][2]]
        outcomes.append(snapshot())
        for i, row in enumerate(xyz(mol.GetNumAtoms(), p["seed"] + 9)):
            mol.GetConformer().SetAtomPosition(i, row)
        # ForceField.cpp::calcEnergy() does not reset the cached distance
        # matrix; calcEnergy(double *pos) does. Refresh the native cache after
        # writing its coordinate pointers, as CK set_positions_ does, without
        # rebuilding contributions or nonbonded pair selection.
        if p["coordinate_update"] != "initialize_distances":
            raise ValueError("invalid coordinate-update policy")
        f.Initialize()
        outcomes.append(snapshot())
        status = f.Minimize(maxIts=p["iterations"], forceTol=p["force_tolerance"], energyTol=p["energy_tolerance"])
        outcomes.extend([status == 0, snapshot(), fixed])
        # Owned CK handle does not write the source molecule's coordinate block.
        return [True, outcomes, original, original]
    if mmff and p["optimizer_observation"] != "isolated_queries_v1":
        raise ValueError("invalid MMFF optimizer observation policy")
    for _ in range(2):
        options = dict(maxIters=p["iterations"], ignoreInterfragInteractions=p["ignore_interfragment"])
        if mmff:
            options.update(mmffVariant=p["variant"], nonBondedThresh=p["threshold"])
        else:
            options.update(vdwThresh=p["threshold"])
        if task.endswith("conformers"):
            result = (FF.MMFFOptimizeMoleculeConfs(mol, numThreads=1, **options) if mmff
                      else FF.UFFOptimizeMoleculeConfs(mol, numThreads=1, **options))
            result = [list(row) for row in result]
        else:
            status = FF.MMFFOptimizeMolecule(mol, **options) if mmff else FF.UFFOptimizeMolecule(mol, **options)
            # CK's energy query prepares detached state, not the molecule to
            # be optimized next. Keep this native copy alive through CalcEnergy:
            # ForceField holds pointers into its conformer.
            observed = Chem.Mol(mol) if mmff else mol
            result = [status, field(observed).CalcEnergy()]
        outcomes.append([result, conformers(mol)])
    return [True, outcomes, original, conformers(mol)]

def geometry(mol, task, p):
    if task in ("alignment",):
        return alignment(mol, p)
    if task.startswith(("uff_", "mmff_")):
        return forcefield(mol, task, p)
    mol = Chem.Mol(mol)
    if task == "coordinates_2d":
        rdDepictor.Compute2DCoords(mol, canonOrient=p["canonical"], nFlipsPerSample=p["flips"],
                                  nSample=p["samples"], sampleSeed=p["seed"], permuteDeg4Nodes=p["permute"],
                                  useRingTemplates=p["ring_templates"])
        return [state(mol), [[float(mol.GetConformer().GetAtomPosition(i)[j]) for j in (0, 1)] for i in range(mol.GetNumAtoms())]]
    if task == "svg":
        draw = rdMolDraw2D.MolDraw2DSVG(p["width"], p["height"], -1, -1, p["no_freetype"])
        rdMolDraw2D.PrepareAndDrawMolecule(draw, mol)
        draw.FinishDrawing()
        return namespace(draw.GetDrawingText()).encode()
    if task == "bounds_matrix":
        if p["add_hs"]:
            mol = Chem.AddHs(mol)
        return rdDistGeom.GetMoleculeBoundsMatrix(mol).tolist()
    if task == "embedding":
        if p["add_hs"]:
            mol = Chem.AddHs(mol)
        method = p["method"]
        params = {"dg": rdDistGeom.EmbedParameters, "etdg": rdDistGeom.ETDG,
                  "etkdg": rdDistGeom.ETKDG, "etkdg_v2": rdDistGeom.ETKDGv2,
                  "etkdg_v3": rdDistGeom.ETKDGv3}[method]()
        params.randomSeed, params.numThreads = p["seed"], 1
        params.maxIterations, params.useRandomCoords = p["max_iterations"], p["random_coords"]
        params.pruneRmsThresh, params.timeout = p["prune_rms"], p["timeout"]
        params.trackFailures = True
        ids = ([rdDistGeom.EmbedMolecule(mol, params)] if p["count"] == 1 else
               list(rdDistGeom.EmbedMultipleConfs(mol, p["count"], params)))
        return [ids, conformers(mol), list(params.GetFailureCounts())]
    raise ValueError(f"unknown geometry task {task}")
