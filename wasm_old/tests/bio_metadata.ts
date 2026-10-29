import {BioStructure,BioMetadata,BioStructureSourceState,BioCrystalInfo,BioAssembly,BioNcsOperator,BioConnectionKind,BioSoftwareClassification,AtomAddress,ResidueAddress,BioTransform,BioEntityDbRef,BioSiftsUnpResidue,BioRefinementInfo,BioTlsGroup} from 'cosmolkit-generated';
const s=BioStructure.fromMmcif('data_x'),m:BioMetadata=s.metadata(),state:BioStructureSourceState=s.sourceState(),c:BioCrystalInfo|null=s.crystal(),a:BioAssembly[]=s.assemblies(),n:BioNcsOperator[]=s.ncsOperators(),info:Map<string,string>=state.info,conect:Map<number,number[]>=state.conectMap,raw:string[]=state.rawRemarks;declare const address:AtomAddress;const res:ResidueAddress=address.residue,num:number|null=res.sequenceNumber;declare const t:BioTransform;const matrix:number[][]=t.matrix,translation:number[]=t.translation;declare const ref:BioRefinementInfo;const tls:BioTlsGroup[]=ref.tlsGroups;declare const db:BioEntityDbRef;const begin:number|null=db.labelSeqBegin;declare const unp:BioSiftsUnpResidue;const residue:number|null=unp.residue();
// @ts-expect-error metadata snapshots are readonly
m.authors=[];
// @ts-expect-error explicit nullable crystal value
const bad:BioCrystalInfo=s.crystal();
void [m,c,a,n,info,conect,raw,res,num,matrix,translation,tls,begin,residue,bad,BioConnectionKind.Covale,BioSoftwareClassification.Refinement];
