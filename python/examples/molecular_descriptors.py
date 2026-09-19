"""Compose the source-backed molecular descriptor families."""

import cosmolkit as ck


molecule = ck.Molecule.from_smiles("CC(O)c1ccncc1C(=O)NCCCl")

connectivity = {
    "chi_0": ck.calc_chi_0(molecule),
    "chi_3v": ck.calc_chi_3v(molecule),
    "kappa_2": ck.calc_kappa_2(molecule),
    "phi": ck.calc_phi(molecule),
}

counts = {
    "heteroatoms": ck.calc_num_heteroatoms(molecule),
    "rings": ck.calc_num_rings(molecule),
    "heterocycles": ck.calc_num_heterocycles(molecule),
    "stereocenters": ck.calc_num_atom_stereo_centers(molecule),
}

mqns = ck.calc_mqns(molecule)
asa, atom_asa, hydrogen_asa = ck.calc_labute_asa_contributions(molecule)
slogp_vsa = ck.calc_slogp_vsa(molecule)
smr_vsa = ck.calc_smr_vsa(molecule)

assert len(mqns) == 42
assert len(atom_asa) == molecule.num_atoms()
assert len(slogp_vsa) == 12
assert len(smr_vsa) == 10

print(connectivity)
print(counts)
print({"asa": asa, "hydrogen_asa": hydrogen_asa})
