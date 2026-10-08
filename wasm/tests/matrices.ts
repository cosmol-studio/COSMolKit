import {Molecule,DenseMatrix,DistanceMatrixParams,DistanceMatrix3dParams,MatrixError,BondOrder,KekulizeParams} from 'cosmolkit-generated';
const m=Molecule.fromSmiles('CO'),p=new DistanceMatrixParams(true,true),p3=new DistanceMatrix3dParams(null,true);
const a:DenseMatrix=m.distanceMatrix(),c:DenseMatrix=m.distanceMatrixWithParams(p),d:DenseMatrix=m.distanceMatrix3d(),f:DenseMatrix=m.distanceMatrix3dWithParams(p3);
const n:number=a.dimension(),values:Float64Array=a.values(),item:number|null=a.get(0,1),id:number|null=p3.conformerId,bo:boolean=p.useBondOrder,aw:boolean=p.useAtomWeights,aw3:boolean=p3.useAtomWeights;
declare const e:MatrixError;
const domain:string=e.domain,kind:string=e.kind,message:string=e.message,cause:Error|null=e.cause;
const position:number|null=e.position,atom:number|null=e.atom,atomCount:number|null=e.atomCount,first:number|null=e.firstPosition,second:number|null=e.secondPosition,bond:number|null=e.bond,bondCount:number|null=e.bondCount,endpoint:string|null=e.endpoint,order:BondOrder|null=e.order,conf:number|null=e.conformerId,dimension:number|null=e.dimension;
// @ts-expect-error typed topological params
m.distanceMatrixWithParams(p3);
// @ts-expect-error typed 3D params
m.distanceMatrix3dWithParams(new KekulizeParams());
// @ts-expect-error u32 index represented by number
a.get(1n,0);
// @ts-expect-error readonly parameter
p.useBondOrder=false;
