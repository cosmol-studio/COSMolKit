"""Focused adapter checks against native public entrypoints."""
import unittest
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize as T
from descriptors import descriptor
from search_tautomer import tautomer
from geometry import geometry,xyz,with_coordinates,conformers
from rdkit.Chem import rdForceFieldHelpers as FF
from prepare import encode

class NativeTasks(unittest.TestCase):
    def test_complete_descriptor_sequences_and_individual_bins(self):
        profiles = {
            "labute_asa_sequence":dict(calls=[[False,False],[True,False],[True,True],[False,False]]),
            "slogp_vsa_sequence":dict(custom_bins=[-.2,0.,.25,.25,.8],calls="cold_warm_forced_v1"),
            "smr_vsa_sequence":dict(custom_bins=[-.2,0.,.25,.25,.8],calls="cold_warm_forced_v1"),
            "chi_n_v_sequence":dict(forces=[False,False,True],orders=list(range(7))),
            "chi_n_n_sequence":dict(forces=[False,False,True],orders=list(range(7))),
            "slogp_vsa_bins":{}, "smr_vsa_bins":{},
        }
        for task,p in profiles.items():
            result=descriptor(Chem.MolFromSmiles("CC(=O)O"),task,p)
            self.assertEqual(len(result),{"labute_asa_sequence":4,"slogp_vsa_sequence":3,
                "smr_vsa_sequence":2,"slogp_vsa_bins":12,"smr_vsa_bins":10}.get(task,3))
    def test_geometry_recipes_and_native_entrypoints(self):
        self.assertEqual(xyz(2,17),xyz(2,17))
        self.assertNotEqual(xyz(2,17),xyz(2,18))
        mol=Chem.MolFromSmiles("CCO")
        for operation in ("transform","best","coordinate","align","all","conformers"):
            profile=dict(operation=operation,weighted=False,reflect=False,explicit_map=True,
                         symmetrize=True,max_iterations=50,max_matches=1000,seed=17)
            value=geometry(mol,"alignment",profile)
            self.assertEqual(len(value),4)
        for method in ("dg","etdg","etkdg","etkdg_v2","etkdg_v3"):
            value=geometry(mol,"embedding",dict(method=method,count=1,seed=20260306,add_hs=True,
                random_coords=False,max_iterations=100,prune_rms=-1.0,timeout=10))
            self.assertEqual(value[0],[0])
            self.assertEqual(len(value[1][0][2]),9)
        for task in ("uff_force_field","mmff_force_field","uff_optimize","mmff_optimize",
                     "uff_optimize_conformers","mmff_optimize_conformers"):
            value=geometry(mol,task,dict(variant="MMFF94" if task.startswith("mmff") else "UFF",
                iterations=2,threshold=100.0,ignore_interfragment=True,add_hs=True,seed=17,
                fixed_first=False,force_tolerance=1e-4,energy_tolerance=1e-6,
                coordinate_update="initialize_distances",optimizer_observation="isolated_queries_v1"))
            self.assertTrue(value[0])

    def test_mmff_optimizer_observations_do_not_prepare_the_live_input(self):
        # This fixed macrocycle exposes non-idempotent setMMFFAromaticity:
        # an extra MMFF property constructor changes the next optimization.
        text = ("CC1=C(/C=C2\\C(C)=C3/C(C)=C(/C=C4\\C(C)=C5/C(CCC(=O)O)="
                "C(C(C)=C5N4)C=C6\\C(CCC(=O)O)=C(C)C(=C6N2)C=C1)N3)C=C(\\C=C)C")
        for variant in ("MMFF94", "MMFF94s"):
            for task in ("mmff_optimize", "mmff_optimize_conformers"):
                with self.subTest(variant=variant, task=task):
                    mol = Chem.MolFromSmiles(text)
                    p = dict(variant=variant, iterations=2, threshold=100.,
                             ignore_interfragment=True, add_hs=True, seed=17,
                             optimizer_observation="isolated_queries_v1")
                    direct = with_coordinates(Chem.AddHs(mol),17,2 if task.endswith("conformers") else 1)
                    original = conformers(direct)
                    outcomes = []
                    for _ in range(2):
                        options = dict(maxIters=2, mmffVariant=variant, nonBondedThresh=100.,
                                       ignoreInterfragInteractions=True)
                        if task.endswith("conformers"):
                            result = [list(row) for row in FF.MMFFOptimizeMoleculeConfs(direct,numThreads=1,**options)]
                        else:
                            status = FF.MMFFOptimizeMolecule(direct,**options)
                            observed = Chem.Mol(direct)
                            props = FF.MMFFGetMoleculeProperties(observed,mmffVariant=variant)
                            field = FF.MMFFGetMoleculeForceField(observed,props)
                            field.Initialize()
                            result = [status,field.CalcEnergy()]
                        outcomes.append([result,conformers(direct)])
                    expected = [True,outcomes,original,conformers(direct)]
                    self.assertEqual(encode(geometry(mol,task,p)),encode(expected))

    def test_connectivity_descriptors_have_the_declared_float_result(self):
        for task in ("chi_0","chi_1"):
            value=descriptor(Chem.MolFromSmiles("[Na+].[Cl-]"),task,{"result_type":"f64"})
            self.assertIs(type(value),float)
            self.assertEqual(value,0.0)

    def test_source_loaded_v1_is_the_native_v1_catalog(self):
        p=dict(catalog="v1",max_tautomers=1000,max_transforms=1000,
               remove_sp3_stereo=True,remove_bond_stereo=True,
               remove_isotopic_hydrogens=True,reassign_stereo=True)
        for text in ("CC(C)=O","O=c1cccc[nH]1","CC(O)=N","C/C=C(O)/[2H]"):
            m=Chem.MolFromSmiles(text)
            expected=T.GetV1TautomerEnumerator().Enumerate(m)
            actual=tautomer(m,p)
            self.assertEqual(actual[0],int(expected.status))
            self.assertEqual([row[0] for row in actual[3]],[text.encode() for text in expected.smiles])
            self.assertEqual(actual[1],list(expected.modifiedAtoms))
            self.assertEqual(actual[2],list(expected.modifiedBonds))

if __name__=="__main__": unittest.main()
