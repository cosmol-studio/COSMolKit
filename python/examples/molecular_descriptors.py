"""Compose the source-backed molecular descriptor families."""

import cosmolkit


molecule = cosmolkit.Molecule.from_smiles("CC(O)c1ccncc1C(=O)NCCCl")

connectivity = {
    "chi_0": cosmolkit.calc_chi_0(molecule),
    "chi_3v": cosmolkit.calc_chi_3v(molecule),
    "kappa_2": cosmolkit.calc_kappa_2(molecule),
    "phi": cosmolkit.calc_phi(molecule),
}

counts = {
    "heteroatoms": cosmolkit.calc_num_heteroatoms(molecule),
    "rings": cosmolkit.calc_num_rings(molecule),
    "heterocycles": cosmolkit.calc_num_heterocycles(molecule),
    "stereocenters": cosmolkit.calc_num_atom_stereo_centers(molecule),
}

mqns = cosmolkit.calc_mqns(molecule)
asa, atom_asa, hydrogen_asa = cosmolkit.calc_labute_asa_contributions(molecule)
slogp_vsa = cosmolkit.calc_slogp_vsa(molecule)
smr_vsa = cosmolkit.calc_smr_vsa(molecule)

assert len(mqns) == 42
assert len(atom_asa) == molecule.num_atoms()
assert len(slogp_vsa) == 12
assert len(smr_vsa) == 10

print(connectivity)
print(counts)
print({"asa": asa, "hydrogen_asa": hydrogen_asa})
