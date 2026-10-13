import * as ck from 'cosmolkit-generated';
const mol=ck.Molecule.fromSmiles('C[C@H](O)Cl');
const structure=new ck.StructureTagParams();structure.conformerId=-1;structure.replaceExistingTags=false;
const potential=new ck.PotentialStereoParams();potential.clean=true;potential.allowNontetrahedral=false;
const options=new ck.CipLabelOptions();options.atoms=[1];options.bonds=null;options.maxRecursiveIterations=1000;
const result:ck.PotentialStereoResult=mol.potentialStereo(potential);
const rows:ck.PotentialStereoInfo[]=result.stereo;
const cleaned:ck.Molecule|null=result.cleanedMolecule;
const center:ck.PotentialStereoCenter=rows[0].centeredOn;
const controlling:Array<number|null>=rows[0].controllingAtoms;
const tetrahedral:ck.TetrahedralStereo[]=mol.tetrahedralStereo();
const ligand:ck.LigandRef=tetrahedral[0].ligands[0];
mol.withCipLabels({atoms:[1],maxRecursiveIterations:200});
mol.assignChiralTagsFromStructure_({conformerId:0});
mol.assignCipLabelsWithOptions_(options);
// @ts-expect-error atom IDs must be numbers
options.atoms=['1'];
// @ts-expect-error result snapshots are read-only
rows[0].permutation=2;
