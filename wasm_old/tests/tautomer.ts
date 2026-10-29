import * as ck from 'cosmolkit-generated';
const mol=ck.Molecule.fromSmiles('CC(=O)C');
const selected:ck.Molecule=ck.canonicalTautomerFromMolecules([mol,mol],{scorer:v=>v.tautomerScore().total});
ck.canonicalTautomerFromMoleculesWithParams([mol],new ck.TautomerParams());
declare const view:ck.TautomerMoleculeView;
declare const progress:ck.TautomerProgress;
const ownedView:ck.TautomerMoleculeView=view.toOwned();
const ownedProgress:ck.TautomerProgress=progress.toOwned();
const params=new ck.TautomerParams();params.maxTautomers=10;params.removeSp3Stereo=false;
params.callback=(source,progress)=>{const rows:Array<[string,ck.TautomerMoleculeView]>=progress.entries();const score:number=source.tautomerScore().total;return rows.length>0;};
params.scorer=view=>view.tautomerScore().total;
params.scoreParams={terms:[new ck.TautomerScoreTerm('carbonyl','[C]=[O]',7)]};
const result:ck.TautomerEnumeration=mol.enumerateTautomers(params);
const status:ck.TautomerEnumerationStatus=result.status;
const values:ck.Molecule[]=result.iter();const child:ck.Molecule|null=result.get(0);
mol.canonicalTautomer({maxTransforms:100,scorer:view=>view.numAtoms()});
mol.tautomerScore({terms:ck.defaultTautomerScoreTerms()});
// @ts-expect-error scorers must return a number
params.scorer=()=>true;
// @ts-expect-error limits must be numbers
params.maxTransforms='100';
