"""Direct pinned RDKit descriptor calls; all floating bits remain intact."""
from rdkit.Chem import rdMolDescriptors as D, QED

COUNTS = {
    "num_heavy_atoms": D.CalcNumHeavyAtoms, "total_atom_count": D.CalcNumAtoms,
    "lipinski_hba": D.CalcNumLipinskiHBA, "lipinski_hbd": D.CalcNumLipinskiHBD,
    "num_heteroatoms": D.CalcNumHeteroatoms, "num_hba": D.CalcNumHBA,
    "num_hbd": D.CalcNumHBD, "num_amide_bonds": D.CalcNumAmideBonds,
    "num_atom_stereo_centers": D.CalcNumAtomStereoCenters,
    "num_unspecified_atom_stereo_centers": D.CalcNumUnspecifiedAtomStereoCenters,
    "num_rings": D.CalcNumRings, "num_heterocycles": D.CalcNumHeterocycles,
    "num_aromatic_rings": D.CalcNumAromaticRings, "num_saturated_rings": D.CalcNumSaturatedRings,
    "num_aliphatic_rings": D.CalcNumAliphaticRings,
    "num_aromatic_heterocycles": D.CalcNumAromaticHeterocycles,
    "num_aromatic_carbocycles": D.CalcNumAromaticCarbocycles,
    "num_aliphatic_heterocycles": D.CalcNumAliphaticHeterocycles,
    "num_aliphatic_carbocycles": D.CalcNumAliphaticCarbocycles,
    "num_saturated_heterocycles": D.CalcNumSaturatedHeterocycles,
    "num_saturated_carbocycles": D.CalcNumSaturatedCarbocycles,
    "num_spiro_atoms": D.CalcNumSpiroAtoms, "num_bridgehead_atoms": D.CalcNumBridgeheadAtoms,
}
SCALARS = {"chi_0":D.CalcChi0n, "chi_1":D.CalcChi1n,
           "fraction_csp3":D.CalcFractionCSP3, "kappa_1":D.CalcKappa1,
           "kappa_2":D.CalcKappa2, "kappa_3":D.CalcKappa3, "phi":D.CalcPhi,
           "qed":QED.qed}

def descriptor(mol, task, p):
    if task == "labute_asa_sequence":
        return [D.CalcLabuteASA(mol,include,force) for include,force in p["calls"]]
    if task in ("slogp_vsa_sequence","smr_vsa_sequence"):
        cold = D.SlogP_VSA_(mol)
        warm = D.SlogP_VSA_(mol)
        forced = D.SlogP_VSA_(mol,p["custom_bins"],True)
        return [cold,warm,forced] if task=="slogp_vsa_sequence" else [D.SMR_VSA_(mol),D.SMR_VSA_(mol,p["custom_bins"],True)]
    if task in ("chi_n_v_sequence","chi_n_n_sequence"):
        fn = D.CalcChiNv if task=="chi_n_v_sequence" else D.CalcChiNn
        return [[fn(mol,order,force) for order in p["orders"]] for force in p["forces"]]
    if task in ("slogp_vsa_bins","smr_vsa_bins"):
        from rdkit.Chem import MolSurf
        prefix,count = ("SlogP_VSA",12) if task=="slogp_vsa_bins" else ("SMR_VSA",10)
        return [getattr(MolSurf,prefix+str(i))(mol) for i in range(1,count+1)]
    if task in COUNTS: return COUNTS[task](mol)
    if task in SCALARS:
        if task == "chi_0":
            from rdkit.Chem.GraphDescriptors import Chi0
            return float(Chi0(mol))
        if task == "chi_1":
            from rdkit.Chem.GraphDescriptors import Chi1
            return float(Chi1(mol))
        return SCALARS[task](mol)
    if task == "molecular_weight": return D._CalcMolWt(mol,p["only_heavy"])
    if task == "exact_molecular_weight": return D.CalcExactMolWt(mol,p["only_heavy"])
    if task == "molecular_formula": return D.CalcMolFormula(mol,p["separate_isotopes"],p["abbreviate_h_isotopes"])
    if task.startswith("chi_"):
        _,order,kind = task.split("_")
        if order == "n": return getattr(D,"CalcChiN"+kind)(mol,p["order"],p["force"])
        return getattr(D,"CalcChi"+order+kind)(mol,p["force"])
    if task == "hall_kier_alpha":
        rows=[0.0]*mol.GetNumAtoms()
        return [D.CalcHallKierAlpha(mol,rows),rows]
    if task == "mqns": return list(D.MQNs_(mol,p["force"]))
    if task == "crippen": return D.CalcCrippenDescriptors(mol,p["include_hydrogens"],p["force"])
    if task == "labute_asa":
        total=D.CalcLabuteASA(mol,p["include_hydrogens"],p["force"])
        rows,h=D._CalcLabuteASAContribs(mol,p["include_hydrogens"],p["force"])
        return [total,total,list(rows),h]
    if task == "tpsa": return D.CalcTPSA(mol,p["force"],p["include_sulfur_phosphorus"])
    if task in ("slogp_vsa","smr_vsa"):
        fn=D.SlogP_VSA_ if task=="slogp_vsa" else D.SMR_VSA_
        return fn(mol,[] if p["bins"] is None else p["bins"],p["force"])
    if task == "num_rotatable_bonds":
        name={"default":"Default","non_strict":"NonStrict","strict":"Strict","strict_linkages":"StrictLinkages"}[p["definition"]]
        return D.CalcNumRotatableBonds(mol,getattr(D.NumRotatableBondsOptions,name))
    raise ValueError(f"unknown descriptor {task}")
