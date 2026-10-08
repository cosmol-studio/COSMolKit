import {Molecule,MmffPropertiesParams,MmffProperties,MmffAtomProperties,MmffVariant,MmffMolPropertiesError,UffParameterQueryError,UffParameterError,UffParameterErrorKind} from 'cosmolkit-generated';
const m=Molecule.fromSmiles('CCO'),p=new MmffPropertiesParams('MMFF94s'),defaultP=new MmffPropertiesParams();const uff:boolean=m.uffHasAllMoleculeParams(),mmff:boolean=m.mmffHasAllMoleculeParams(),a:MmffProperties=m.mmffProperties(),b:MmffProperties=m.mmffPropertiesWithParams(p),valid:boolean=a.isValid(),variant:MmffVariant=a.variant(),atoms:MmffAtomProperties[]=a.atoms();const type:number=a.atomType(0),formal:number=a.formalCharge(0),partial:number=a.partialCharge(0),at:number=atoms[0].atomType(),af:number=atoms[0].formalCharge(),ap:number=atoms[0].partialCharge();declare const e:UffParameterError;const kind:UffParameterErrorKind=e.kind();declare const me:MmffMolPropertiesError,ue:UffParameterQueryError;const cause:Error|null=ue.cause;
// @ts-expect-error no numeric coercion
new MmffPropertiesParams(94);
// @ts-expect-error immutable parameter
p.mmffVariant='MMFF94';
// @ts-expect-error atom index is a number
a.atomType('0');
