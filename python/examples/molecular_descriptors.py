"""Compose the source-backed molecular descriptor families."""

import cosmolkit as ck


molecule = ck.mol_from_smiles("CC(O)c1ccncc1C(=O)NCCCl")

connectivity = {
    "chi_0": molecule.chi_0(),
    "chi_3v": molecule.chi_3_v(),
    "kappa_2": molecule.kappa_2(),
    "phi": molecule.phi(),
}

counts = {
    "heteroatoms": molecule.num_heteroatoms(),
    "rings": molecule.num_rings(),
    "heterocycles": molecule.num_heterocycles(),
    "stereocenters": molecule.num_atom_stereo_centers(),
}

mqns = molecule.mqns()
contributions = molecule.labute_asa_contributions()
asa = contributions.asa
atom_asa = contributions.atom_contributions
hydrogen_asa = contributions.hydrogen_contribution
slogp_vsa = molecule.slogp_vsa()
smr_vsa = molecule.smr_vsa()

assert len(mqns) == 42
assert len(atom_asa) == molecule.num_atoms()
assert len(slogp_vsa) == 12
assert len(smr_vsa) == 10

print(connectivity)
print(counts)
print({"asa": asa, "hydrogen_asa": hydrogen_asa})
