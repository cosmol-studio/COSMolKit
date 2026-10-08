import {BioStructure,Protein,BioMoleculeParams,BioMoleculeError,BioMoleculeConversionError,Molecule} from 'cosmolkit-generated';
const p=new BioMoleculeParams(false,false,4294967295,false);const sanitize:boolean=p.sanitize,remove:boolean=p.removeHs,flavor:number=p.flavor,proximity:boolean=p.proximityBonding;
for(const s of [BioStructure.fromPdb('END\n'),Protein.fromPdb('END\n')]){const m:Molecule=s.toMolecule(),q:Molecule=s.toMoleculeWithParams(p);}
declare const e:BioMoleculeError,c:BioMoleculeConversionError;const domain:string=e.domain,kind:string=c.kind,cause:Error|null=c.cause;
// @ts-expect-error flavor remains a checked number
new BioMoleculeParams(true,true,1n);
// @ts-expect-error full params require the typed projection
BioStructure.fromPdb('END\n').toMoleculeWithParams({sanitize:false});
// @ts-expect-error result is not a BIO hierarchy
const s:BioStructure=Protein.fromPdb('END\n').toMolecule();
