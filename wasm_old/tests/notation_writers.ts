import {Molecule,CxSmilesFields,CxCoordinateSelection,CxSmilesWriteParams,FragmentSmilesWriteParams,FragmentCxSmilesWriteParams,RandomSmilesWriteParams,SmartsWriteParams} from 'cosmolkit-generated';
const mol=Molecule.fromSmiles('CCO');
const params=new CxSmilesWriteParams();
params.fields=CxSmilesFields.NONE;
params.smiles={allBondsExplicit:true};
params.coordinateSelection=CxCoordinateSelection.twoD(0);
const cx:string=mol.toCxSmiles(params);
mol.toCxSmiles({fields:CxSmilesFields.ALL});
const fragment=new FragmentSmilesWriteParams([0,1]);
fragment.bonds=null;fragment.atomSymbols=['C','C','O'];
mol.toFragmentSmilesWithParams(fragment);
mol.toFragmentSmiles([0,1],{smiles:{canonical:false}});
mol.toFragmentCxSmiles([0,1],new FragmentCxSmilesWriteParams([0,1]));
const random:string[]=mol.toRandomSmiles(3,42,new RandomSmilesWriteParams());
mol.toRandomSmiles(3,42,{kekule:true});
mol.toSmarts({rootedAtAtom:0});mol.toCxSmarts(new SmartsWriteParams());
// @ts-expect-error index sequence must contain numbers
fragment.atoms=['0'];
// @ts-expect-error annotation mask is a typed value
params.fields=1;
// @ts-expect-error random seed is a number
mol.toRandomSmiles(3,'42');
void [cx,random];
