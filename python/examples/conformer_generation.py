"""Public Python API example: 3D conformer generation and force-field cleanup."""

from __future__ import annotations

import numpy as np

from cosmolkit import EmbedParameters, Molecule


mol = Molecule.from_smiles("CC(=O)NC").with_hydrogens()

single_params = EmbedParameters.etkdg_v3()
single_params.random_seed = 0xF00D
single_params.num_threads = 1
single_params.track_failures = True

single_result = mol.with_3d_conformer_result(single_params)
embedded = single_result.molecule()
print("single conformers:", embedded.num_conformers())
print("single conf id:", single_result.conf_id())
print("failure counters:", single_result.params().failures)

multi_params = EmbedParameters.etkdg()
multi_params.random_seed = 123
multi_params.num_threads = 1
multi_params.prune_rms_thresh = 0.5
multi_params.enable_sequential_random_seeds = True

multi_result = mol.with_3d_conformers_result(5, multi_params)
multi = multi_result.molecule()
print("pruned conformers:", multi.num_conformers())
print("kept conformer ids:", multi_result.conf_ids())

if embedded.has_uff_params():
    start = embedded.coordinates_3d().copy()
    uff = embedded.with_uff_optimized(max_iters=200)
    print("UFF converged:", not uff.needs_more())
    print("UFF status code:", uff.status_code())
    print("UFF energy:", uff.energy())
    print("UFF moved coordinates:", not np.allclose(start, uff.molecule().coordinates_3d()))

if embedded.has_mmff_params():
    mmff = embedded.with_mmff_optimized(max_iters=200)
    print("MMFF94 converged:", not mmff.needs_more())
    print("MMFF94 status code:", mmff.status_code())
