import {Molecule,SdfRecord,QueryGraph,MoleculeProperties,MolBlockWriteParams,MolCoordinateSelection,SdfFormat,MolWriteError,CoordinateDimension} from 'cosmolkit-generated';
const s=MolCoordinateSelection.TwoD(0),auto=MolCoordinateSelection.Auto,three=MolCoordinateSelection.ThreeD(1),id:number|null=s.id,kind:string=s.kind,p=new MolBlockWriteParams(SdfFormat.V3000,false,true,true,6,s,true),format:SdfFormat=p.format,selection:MolCoordinateSelection=p.coordinateSelection;const force:boolean=p.force2d,stereo:boolean=p.includeStereo,kek:boolean=p.kekulize,precision:number=p.precision,coordinates:boolean=p.includeCoordinates;declare const props:MoleculeProperties;const r=SdfRecord.fromQueryGraph(QueryGraph.fromSmarts('*'),props),m=Molecule.fromSmiles('CC');for(const object of [m,r]){const mol:string=object.toMol(),sdf:string=object.toSdf(),customMol:string=object.toMolWithParams(p),customSdf:string=object.toSdfWithParams(p);}const two:string=m.toSdf2d(),threeText:string=m.toSdf3d(),twoCustom:string=m.toSdf2dWithParams(p),threeCustom:string=m.toSdf3dWithParams(p),writeMol:void=m.writeMol(''),writeSdf:void=m.writeSdf(''),writeMolCustom:void=m.writeMolWithParams('',p),writeSdfCustom:void=m.writeSdfWithParams('',p),path:string=m.writeSdfFiles('',null),customPath:string=m.writeSdfFilesWithParams('','x',p);
declare const e:MolWriteError;const domain:string=e.domain,errorKind:string=e.kind,message:string=e.message,cause:Error|null=e.cause,subset:string|null=e.subset,detail:string|null=e.valueDetail,twoCount:number|null=e.twoD,threeCount:number|null=e.threeD,dimension:CoordinateDimension|null=e.dimension,missingId:number|null=e.id;
// @ts-expect-error full typed policy required
m.toMolWithParams({});
// @ts-expect-error nullable id
const definite:number=auto.id;
// @ts-expect-error explicit nullable file name required
m.writeSdfFiles('');
// @ts-expect-error unsigned id is number
MolCoordinateSelection.TwoD(1n);
