import assert from 'node:assert/strict';import {readFileSync} from 'node:fs';import test from 'node:test';import {pathToFileURL} from 'node:url';
const b=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);b.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
test('Topological distance matrices transport all weights, row-major values and independent result ownership',()=>{
 const defaults=new b.DistanceMatrixParams();assert.deepEqual([defaults.useBondOrder,defaults.useAtomWeights],[false,false]);
 const m=b.Molecule.fromSmiles('C=O'),before=m.toSmiles();assert.deepEqual(Array.from(m.distanceMatrix().values()),[0,1,1,0]);
 for(const bo of [false,true])for(const aw of [false,true]){
  const p=new b.DistanceMatrixParams(bo,aw),matrix=m.distanceMatrixWithParams(p),expected=[aw?1:0,bo?0.5:1,bo?0.5:1,aw?0.75:0];
  assert.ok(matrix instanceof b.DenseMatrix);assert.equal(matrix.dimension(),2);assert.deepEqual(Array.from(matrix.values()),expected);
  for(let i=0;i<2;i++)for(let j=0;j<2;j++)assert.equal(matrix.get(i,j),expected[i*2+j]);
  assert.equal(matrix.get(2,0),null);assert.equal(matrix.get(0,4294967295),null);
  const values=matrix.values();values[0]=99;assert.deepEqual(Array.from(matrix.values()),expected);
  assert.deepEqual(Array.from(m.distanceMatrixWithParams(p).values()),expected);assert.equal(m.toSmiles(),before);
 }
 const chain=b.Molecule.fromSmiles('CCO').distanceMatrix();assert.deepEqual(Array.from(chain.values()),[0,1,2,1,0,1,2,1,0]);
 assert.equal(b.Molecule.fromSmiles('C.O').distanceMatrix().get(0,1),100000000);
 const empty=b.Molecule.new().distanceMatrix();assert.equal(empty.dimension(),0);assert.deepEqual(Array.from(empty.values()),[]);assert.equal(empty.get(0,0),null);
 assert.equal(b.Molecule.fromSmiles('*').distanceMatrixWithParams(new b.DistanceMatrixParams(false,true)).get(0,0),Infinity);
 const detached=m.distanceMatrix();m.free();assert.equal(detached.get(0,1),1);
});
test('3D matrix parameters select exact conformers and preserve canonical geometry and typed errors',()=>{
 const p=new b.DistanceMatrix3dParams();assert.equal(p.conformerId,null);assert.equal(p.useAtomWeights,false);
 const m=b.Molecule.fromXyzBlock('3\nmatrix\nC 0 0 0\nN 3 0 0\nO 3 4 0\n'),coords=Array.from(m.coordinates3d(0));
 assert.deepEqual(Array.from(m.distanceMatrix3d().values()),[0,3,5,3,0,4,5,4,0]);
 for(const id of [null,0])for(const weights of [false,true]){
  const params=new b.DistanceMatrix3dParams(id,weights),matrix=m.distanceMatrix3dWithParams(params);
  assert.deepEqual(Array.from(matrix.values()),[weights?1:0,3,5,3,weights?6/7:0,4,5,4,weights?0.75:0]);assert.equal(params.conformerId,id);assert.deepEqual(Array.from(m.distanceMatrix3dWithParams(params).values()),Array.from(matrix.values()));
 }
 for(const [fn,kind,id] of [[()=>b.Molecule.fromSmiles('C').distanceMatrix3d(),'No3dConformer',null],[()=>m.distanceMatrix3dWithParams(new b.DistanceMatrix3dParams(99)),'ConformerNotFound',99]])assert.throws(fn,e=>{
  assert.equal(e.name,'MatrixError');assert.equal(e.domain,'matrices');assert.equal(e.kind,kind);assert.equal(e.conformerId,id);assert.equal(e.cause,null);assert.ok(e.detail instanceof b.MatrixError);assert.equal(e.detail.kind,kind);assert.equal(e.detail.conformerId,id);assert.equal(e.detail.message,e.message);assert.equal(e.detail.cause,null);
  for(const key of ['position','atom','atomCount','firstPosition','secondPosition','bond','bondCount','endpoint','order','dimension'])assert.equal(e.detail[key],null);return true;
 });
 assert.deepEqual(Array.from(m.coordinates3d(0)),coords);
});
test('Matrix host boundaries reject coercion and retain full canonical bond-order vocabulary',()=>{
 const m=b.Molecule.fromSmiles('C'),matrix=m.distanceMatrix();
 for(const n of [-1,1.5,4294967296,NaN,Infinity]){assert.throws(()=>matrix.get(n,0),RangeError);assert.throws(()=>matrix.get(0,n),RangeError);assert.throws(()=>new b.DistanceMatrix3dParams(n),RangeError);}
 for(const value of [{},null,undefined,new b.KekulizeParams()]){assert.throws(()=>m.distanceMatrixWithParams(value),TypeError);assert.throws(()=>m.distanceMatrix3dWithParams(value),TypeError);}
 assert.throws(()=>new b.DistanceMatrixParams('false'),TypeError);assert.throws(()=>new b.DistanceMatrixParams(false,1),TypeError);assert.throws(()=>new b.DistanceMatrix3dParams(null,'true'),TypeError);
 assert.throws(()=>{new b.DistanceMatrixParams().useBondOrder=true;},TypeError);
 const names=['Unspecified','Single','Double','Triple','Quadruple','Quintuple','Hextuple','OneAndHalf','TwoAndHalf','ThreeAndHalf','FourAndHalf','FiveAndHalf','Aromatic','Ionic','Hydrogen','ThreeCenter','DativeOne','Dative','DativeLeft','DativeRight','Other','Zero'];names.forEach((name,i)=>assert.equal(b.BondOrder[name],i));
});
