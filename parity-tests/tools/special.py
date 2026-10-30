"""Fixed input adapters call native public APIs directly; no old runner imports."""
import json
import tempfile
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
import bio
import search_tautomer

TASKS = {"tetrahedral_geometry", "structure_tags", "tautomer_focused", "tautomer_long_conjugated", "mmff_builtin",
         "forcefield_properties", "mcs_upstream", "mcs_jnk1", "bio_mmcif_switches", "conformer_fixed19", "conformer_library", "molalign_focused", "forcefield_optimizers"}


def structure_tags(row):
    import struct
    from geometry import conformers

    def set_property(obj, key, value):
        if isinstance(value, bool): obj.SetBoolProp(key, value)
        elif isinstance(value, str): obj.SetProp(key, value)
        elif key == "_chiralPermutation": obj.SetUnsignedProp(key, value)
        elif isinstance(value, int): obj.SetIntProp(key, value)
        else: obj.SetDoubleProp(key, value)

    def properties(obj):
        return [[k.encode(), v] for k, v in sorted(obj.GetPropsAsDict(includePrivate=True, includeComputed=False, autoConvertStrings=False).items())]

    def snapshot(m):
        atoms = [[a.GetAtomicNum(), a.GetFormalCharge(), a.GetNumExplicitHs(), a.GetNumImplicitHs(),
                  int(a.GetChiralTag()), properties(a)] for a in m.GetAtoms()]
        bonds = [[b.GetBeginAtomIdx(), b.GetEndAtomIdx(), int(b.GetBondType()), int(b.GetBondDir()),
                  properties(b)] for b in m.GetBonds()]
        return [atoms, bonds, conformers(m), properties(m)]

    m = Chem.RWMol()
    for r in row["atoms"]:
        a = Chem.Atom(r["atomic_number"])
        a.SetFormalCharge(r.get("formal_charge", 0))
        a.SetNumExplicitHs(r.get("explicit_hydrogens", 0))
        a.SetChiralTag(getattr(Chem.ChiralType, r.get("chiral_tag", "CHI_UNSPECIFIED")))
        if r.get("chiral_permutation") is not None: a.SetUnsignedProp("_chiralPermutation", r["chiral_permutation"])
        for k, v in r.get("props", {}).items(): set_property(a, k, v)
        m.AddAtom(a)
    for r in row["bonds"]:
        m.AddBond(r["begin"], r["end"], getattr(Chem.BondType, r.get("type", "SINGLE")))
        b = m.GetBondWithIdx(m.GetNumBonds()-1)
        b.SetBondDir(getattr(Chem.BondDir, r.get("direction", "NONE")))
        if r.get("unknown_stereo") is not None: b.SetIntProp("_UnknownStereo", r["unknown_stereo"])
        for k, v in r.get("props", {}).items(): set_property(b, k, v)
    for r in row["conformers"]:
        c = Chem.Conformer(m.GetNumAtoms())
        c.SetId(r["id"])
        c.Set3D(r["is_3d"])
        for i, bits in enumerate(r["coordinate_bits"]):
            c.SetAtomPosition(i, [struct.unpack("<d", struct.pack("<Q", v))[0] for v in bits])
        m.AddConformer(c, assignId=False)
    for k, v in row["molecule_props"].items(): set_property(m, k, v)
    m.UpdatePropertyCache(strict=False)
    before = snapshot(m)
    previous = Chem.GetAllowNontetrahedralChirality()
    environment = row["environment"]
    Chem.SetAllowNontetrahedralChirality(environment["mode"] == "unset" or environment["value"] != "0")
    try:
        try:
            Chem.AssignAtomChiralTagsFromStructure(m, confId=row["conf_id"], replaceExistingTags=row["replace_existing_tags"])
            status = True
        except ValueError as e:
            if str(e) != "Bad Conformer Id": raise
            status = ["conformer_not_found", row["conf_id"]]
            return [status, before, None]
        except RuntimeError as e:
            if str(e) != "Cannot normalize a zero length vector": raise
            return [["zero_length_vector"], before, None]
        return [status, before, snapshot(m)]
    finally:
        Chem.SetAllowNontetrahedralChirality(previous)


def optimizer(row):
    from geometry import conformers
    m = Chem.MolFromSmiles(row["single_cxsmiles"])
    mmff = row["case"]["forcefield"] == "mmff"
    def field(mol):
        if mmff:
            return AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94"), nonBondedThresh=100.0, confId=0, ignoreInterfragInteractions=True)
        return AllChem.UFFGetMoleculeForceField(mol, vdwThresh=100.0, confId=0, ignoreInterfragInteractions=True)
    initial_m = Chem.Mol(m)
    f = field(initial_m)
    f.Initialize()
    f.Minimize(maxIts=0)
    initial = [float(f.CalcEnergy()), list(f.CalcGrad()), conformers(m)]
    single_m = Chem.Mol(m)
    f = field(single_m)
    f.Initialize()
    status = f.Minimize(maxIts=2)
    single = [status, float(f.CalcEnergy()), conformers(single_m)]
    multi = Chem.Mol(m)
    multi.RemoveAllConformers()
    for id, is3d, xyz in row["multi_conformers"]:
        c = Chem.Conformer(m.GetNumAtoms())
        c.SetId(id)
        c.Set3D(is3d)
        for i,p in enumerate(xyz): c.SetAtomPosition(i,p)
        multi.AddConformer(c,assignId=False)
    before = conformers(multi)
    if mmff:
        result = AllChem.MMFFOptimizeMoleculeConfs(multi,numThreads=1,maxIters=2,mmffVariant="MMFF94",nonBondedThresh=100.0,ignoreInterfragInteractions=True)
        available = AllChem.MMFFHasAllMoleculeParams(Chem.Mol(m))
    else:
        result = AllChem.UFFOptimizeMoleculeConfs(multi,numThreads=1,maxIters=2,vdwThresh=100.0,ignoreInterfragInteractions=True)
        available = AllChem.UFFHasAllMoleculeParams(Chem.Mol(m))
    return [available,initial,single,before,[[list(r) for r in result],conformers(multi)],conformers(m)]


def alignment_molecule(row):
    m = Chem.MolFromSmiles(row["smiles"])
    m.RemoveAllConformers()
    for r in row["conformers"]:
        c = Chem.Conformer(m.GetNumAtoms())
        c.SetId(r["id"])
        c.Set3D(True)
        for i, xyz in enumerate(r["coordinates"]): c.SetAtomPosition(i, xyz)
        m.AddConformer(c, assignId=False)
    return m


def alignment(case):
    from rdkit.Chem import rdMolAlign as A
    from geometry import conformers
    source = alignment_molecule(case.get("probe", case.get("molecule")))
    reference = alignment_molecule(case["reference"]) if "reference" in case else None
    before, before_ref = conformers(source), None if reference is None else conformers(reference)
    out = []
    for p in case["calls"]:
        probe, after = Chem.Mol(source), before
        ids = dict(prbCid=p.get("probe_conformer_id", -1), refCid=p.get("reference_conformer_id", -1))
        weights = p.get("weights", [])
        basic = dict(weights=weights, reflect=p.get("reflect", False), maxIters=p.get("max_iterations", 50))
        best = dict(map=p.get("atom_maps", []), maxMatches=p.get("max_matches", 1_000_000),
                    symmetrizeConjugatedTerminalGroups=p.get("symmetrize_conjugated_terminal_groups", True), weights=weights)
        operation = p["operation"]
        try:
            if operation == "alignment_transform":
                rms, matrix = A.GetAlignmentTransform(probe, reference, atomMap=p.get("atom_map", []), **ids, **basic)
                mapping = p.get("atom_map") or list(enumerate(probe.GetSubstructMatch(reference)))
                result = [rms, matrix.tolist(), [list(pair) for pair in mapping]]
            elif operation == "best_alignment":
                rms, matrix, mapping = A.GetBestAlignmentTransform(probe, reference, reflect=basic["reflect"], maxIters=basic["maxIters"], numThreads=p.get("num_threads", 1), **ids, **best)
                result = [rms, matrix.tolist(), [list(pair) for pair in mapping]]
            elif operation == "coordinate_rmsd":
                result = A.CalcRMS(probe, reference, prbId=ids["prbCid"], refId=ids["refCid"], **best)
            elif operation == "all_conformer_best_rms":
                rms = A.GetAllConformerBestRMS(probe, numThreads=p.get("num_threads", 1), **best)
                confs = list(probe.GetConformers())
                pairs = [(confs[i].GetId(), confs[j].GetId()) for i in range(len(confs)) for j in range(i)]
                result = [[i,j,v] for (i,j),v in zip(pairs,rms)]
            elif operation == "align_conformers":
                result = []
                args = dict(atomIds=p.get("atom_indices", []), confIds=p.get("conformer_ids", []), RMSlist=result, **basic)
                A.AlignMolConformers(probe, **args)
                after = conformers(probe)
            else:
                raise ValueError("unknown alignment operation")
            result = [True, result]
        except Exception as e:
            message = str(e)
            if isinstance(e, ValueError) and message == "Bad Conformer Id":
                present = {c.GetId() for c in probe.GetConformers()}
                requested = ids["prbCid"] if ids["prbCid"] >= 0 and ids["prbCid"] not in present else ids["refCid"]
                result = [False, ["conformer_not_found", requested]]
            elif "Incorrect number of weights specified" in message:
                result = [False, ["weight_count_mismatch", len(p["atom_map"]), len(weights)]]
            elif message == "No sub-structure match found between the reference and probe mol":
                result = [False, ["no_substructure_match"]]
            else:
                raise
        out.append([result, after, before, before_ref])
    return out


def embedding(task, row):
    from rdkit.Chem import rdDistGeom as D
    from rdkit.Geometry import Point3D
    from geometry import conformers
    c = row["case"] if task == "conformer_fixed19" else {"source":row["case"], "source_kind":"smiles_with_hydrogens", "preset_name":"ETKDGv3", "attrs":{"maxIterations":3,"randomSeed":61453,"numThreads":1,"timeout":0}, "mode":"single"}
    if c["source_kind"] == "fixture_mol":
        m = Chem.MolFromMolBlock(c["mol_block"], sanitize=True, removeHs=False)
    else:
        m = Chem.MolFromSmiles(c["source"])
        if m is not None and c["source_kind"] == "smiles_with_hydrogens":
            m = Chem.AddHs(m)
    if m is None:
        return None
    p = D.EmbedParameters() if c["preset_name"] == "DG" else getattr(D, c["preset_name"])()
    for k,v in c["attrs"].items(): setattr(p,k,v)
    if "coord_map" in c:
        p.SetCoordMap({int(k):Point3D(*v) for k,v in c["coord_map"].items()})
    if "cpci" in c:
        p.SetCPCI({tuple(map(int,k.split("-"))):v for k,v in c["cpci"].items()})
    status = int(D.EmbedMolecule(m,p)) if c["mode"] == "single" else list(D.EmbedMultipleConfs(m,c["num_confs"],p))
    return [status, conformers(m), list(p.GetFailureCounts()) if p.trackFailures else None]


def mmff(mol, variant):
    available = AllChem.MMFFHasAllMoleculeParams(Chem.Mol(mol))
    props = AllChem.MMFFGetMoleculeProperties(Chem.Mol(mol), mmffVariant=variant)
    atoms = None if props is None else [[props.GetMMFFAtomType(i), props.GetMMFFFormalCharge(i),
                                       props.GetMMFFPartialCharge(i)] for i in range(mol.GetNumAtoms())]
    return [mol.GetNumAtoms(), available, props is not None, atoms]


def mcs(row):
    mols = []
    for recipe in row["molecules"]:
        if recipe["format"] == "smiles":
            p = Chem.SmilesParserParams()
            p.sanitize, p.removeHs = recipe["sanitize"], recipe["remove_hs"]
            m = Chem.MolFromSmiles(recipe["text"], p)
        else:
            m = Chem.MolFromMolBlock(recipe["text"], sanitize=recipe["sanitize"], removeHs=recipe["remove_hs"])
        if m is None:
            raise ValueError("native rejected fixed MCS input")
        mols.append(m)
    params = rdFMCS.MCSParameters()
    fields = {"store_all":"StoreAll", "maximize_bonds":"MaximizeBonds", "threshold":"Threshold",
              "timeout":"Timeout", "verbose":"Verbose", "initial_seed":"InitialSeed"}
    atom = {"match_valences":"MatchValences", "match_chiral_tag":"MatchChiralTag", "match_formal_charge":"MatchFormalCharge",
            "ring_matches_ring_only":"RingMatchesRingOnly", "complete_rings_only":"CompleteRingsOnly", "match_isotope":"MatchIsotope", "max_distance":"MaxDistance"}
    bond = {"ring_matches_ring_only":"RingMatchesRingOnly", "complete_rings_only":"CompleteRingsOnly", "match_fused_rings":"MatchFusedRings",
            "match_fused_rings_strict":"MatchFusedRingsStrict", "match_stereo":"MatchStereo"}
    for k, v in row["case"]["parameters"].items():
        if k in fields:
            setattr(params, fields[k], v)
        elif k == "atom_comparator":
            params.AtomTyper = {"any":rdFMCS.AtomCompare.CompareAny,"elements":rdFMCS.AtomCompare.CompareElements,
                               "isotopes":rdFMCS.AtomCompare.CompareIsotopes,"any_heavy_atom":rdFMCS.AtomCompare.CompareAnyHeavyAtom}[v]
        elif k == "bond_comparator":
            params.BondTyper = {"any":rdFMCS.BondCompare.CompareAny,"order":rdFMCS.BondCompare.CompareOrder,"order_exact":rdFMCS.BondCompare.CompareOrderExact}[v]
        else:
            target, names = (params.AtomCompareParameters, atom) if k == "atom_compare_parameters" else (params.BondCompareParameters, bond)
            for name, value in v.items():
                setattr(target, names[name], value)
    r = rdFMCS.FindMCS(mols, params)
    if r.canceled:
        return [False]
    q = r.queryMol
    query = None if q is None else [q.GetNumAtoms(), q.GetNumBonds(), Chem.MolToSmarts(q).encode(), [m.HasSubstructMatch(q) for m in mols]]
    return [True, r.numAtoms, r.numBonds, r.smartsString.encode(), [s.encode() for s in sorted(r.degenerateSmartsQueryMolDict)], query]


def cif(row):
    import gemmi
    if gemmi.__version__ != bio.GEMMI:
        raise RuntimeError("unexpected Gemmi version")
    flags = ("atoms block_name entry database_status author cell symmetry entity entity_poly struct_ref chem_comp exptl diffrn reflns refine title_keywords ncs struct_asym origx struct_conf struct_sheet struct_biol assembly conn cis scale atom_type entity_poly_seq tls software group_pdb auth_all").split()
    with tempfile.TemporaryDirectory() as folder:
        path = Path(folder) / Path(row["case"]["input"]).name
        path.write_text(row["text"])
        structure = gemmi.read_structure(str(path), merge_chain_parts=False)
    values = []
    for all_groups in (False, True):
        for flag in flags:
            groups = gemmi.MmcifOutputGroups(all_groups)
            setattr(groups, flag, True if flag == "auth_all" else not all_groups)
            values.append(structure.make_mmcif_document(groups).as_string().encode())
    return values


def observe(task, text):
    row = json.loads(text)
    if task == "tetrahedral_geometry":
        mol = Chem.MolFromSmiles(row["smiles"])
        if mol is None:
            raise ValueError("geometry input rejected")
        centers = [a.GetIdx() for a in mol.GetAtoms() if a.GetChiralTag() in
                   (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)]
        params = AllChem.ETKDGv3()
        params.randomSeed = row["seed"]
        embedding_status = AllChem.EmbedMolecule(mol, params)
        if embedding_status != 0:
            return ["embedding_failed", embedding_status]
        positions = mol.GetConformer().GetPositions().tolist()
        # Coordinates are an independent witness, not a CK-generated expectation.
        return [positions, [[center, True] for center in centers],
                [[center, True] for center in centers]]
    if task == "structure_tags":
        return structure_tags(row)
    if task == "forcefield_optimizers":
        return optimizer(row)
    if task == "molalign_focused":
        return alignment(row["case"])
    if task in ("conformer_fixed19", "conformer_library"):
        return embedding(task,row)
    c = row["case"]
    if task in ("mcs_upstream", "mcs_jnk1"):
        return mcs(row)
    if task == "bio_mmcif_switches":
        return cif(row)
    smiles = c if isinstance(c, str) else c["smiles"]
    p = Chem.SmilesParserParams()
    if task.startswith("tautomer_"):
        p.sanitize, p.removeHs = c["sanitize"], c["remove_hs"]
    m = Chem.MolFromSmiles(smiles, p)
    if m is None:
        return None
    if task.startswith("tautomer_"):
        return [search_tautomer.tautomer(m, {k:v for k,v in branch.items() if k != "name"}) for branch in row["branches"]]
    if task == "mmff_builtin":
        return mmff(m, c["variant"])
    if task == "forcefield_properties":
        return [[AllChem.UFFHasAllMoleculeParams(Chem.Mol(x)), mmff(x, "MMFF94")] for x in (m, Chem.AddHs(m))]
    raise ValueError(f"unknown special task {task}")
