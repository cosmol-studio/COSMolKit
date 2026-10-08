import assert from "node:assert/strict";
import {readFileSync} from "node:fs";
import test from "node:test";
import {pathToFileURL} from "node:url";
const b=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);b.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
const batch=()=>b.MoleculeBatch.fromSmilesListWithParams(["CCO","[","O"],new b.SmilesParseParams(),new b.BatchParams(b.BatchErrorMode.KeepErrors));

test("full SMILES writer and query parameters retain frozen source defaults and checked values",()=>{
    const write=new b.SmilesWriteParams();
    const fields=["doIsomericSmiles","doKekule","canonical","cleanStereo","allBondsExplicit","allHydrogensExplicit","includeDativeBonds","ignoreAtomMapNumbers"],defaults=[true,false,true,true,false,false,true,false];
    fields.forEach((key,i)=>assert.equal(write[key],defaults[i]));assert.equal(write.rootedAtAtom,null);
    const inverse=new b.SmilesWriteParams(false,true,false,false,1,true,true,false,true);fields.forEach((key,i)=>assert.equal(inverse[key],!defaults[i]));assert.equal(inverse.rootedAtAtom,1);
    assert.equal(new b.SmilesWriteParams(undefined,undefined,undefined,undefined,null).rootedAtAtom,null);
    assert.throws(()=>new b.SmilesWriteParams("true"),TypeError);
    for(const n of [-1,1.5,4294967296,NaN,Infinity])assert.throws(()=>new b.SmilesWriteParams(undefined,undefined,undefined,undefined,n),RangeError);
    assert.throws(()=>{write.canonical=false;},TypeError);
    const query=new b.BatchQueryParams();assert.equal(query.nJobs,null);assert.equal(query.progressBar,null);assert.equal(query.progressCallback,null);
    const callback=()=>{};const configured=new b.BatchQueryParams(2,true,callback);assert.equal(configured.nJobs,2);assert.equal(configured.progressBar,true);assert.equal(configured.progressCallback,callback);
    assert.throws(()=>{configured.progressCallback=null;},TypeError);
    assert.throws(()=>new b.BatchQueryParams(null,null,{}),TypeError);
    assert.throws(()=>new b.BatchQueryParams(null,"false"),TypeError);
    assert.throws(()=>new b.BatchQueryParams(-1),RangeError);
});

test("all six batch queries retain ordered nulls and detached nested matrix rows",()=>{
    const value=batch(),query=new b.BatchQueryParams(),write=new b.SmilesWriteParams();
    assert.deepEqual(value.toSmilesList(),["CCO",null,"O"]);assert.deepEqual(value.toSmilesListWithParams(write,query),value.toSmilesList());
    const matrix=value.dgBoundsMatrixList(),explicit=value.dgBoundsMatrixListWithParams(query);assert.deepEqual(explicit,matrix);
    assert.equal(matrix.length,3);assert.equal(matrix[1],null);assert.equal(matrix[0].length,3);assert.equal(matrix[0][0].length,3);assert.equal(matrix[2].length,1);assert.equal(matrix[0][0][0],0);
    const first=matrix[0][0][1];assert.ok(first>0);matrix[0][0][1]=999;assert.equal(value.dgBoundsMatrixList()[0][0][1],first);
    const svg=value.toSvgList(120,100);assert.equal(svg[1],null);assert.ok(svg[0].includes("<svg"));assert.ok(svg[2].includes("<svg"));assert.deepEqual(value.toSvgListWithParams(120,100,query),svg);
    for(const n of [-1,1.5,4294967296,NaN,Infinity])assert.throws(()=>value.toSvgList(n,100),RangeError);
    assert.throws(()=>value.toSvgList(0,100),e=>{assert.equal(e.domain,"batch");assert.deepEqual(e.recordErrors.map(v=>v.index()),[0,2]);assert.equal(e.cause.cause.domain,"drawing");assert.equal(e.cause.cause.kind,"InvalidDimensions");assert.equal(e.cause.cause.width,0);assert.ok(e.cause.cause.detail instanceof b.DrawingError);return true;});
    assert.throws(()=>value.toSmilesListWithParams(new b.SmilesWriteParams(undefined,undefined,undefined,undefined,99),query),e=>e.errors===2&&e.cause.cause.domain==="smiles"&&e.cause.cause.kind==="Write");
    assert.throws(()=>value.dgBoundsMatrixListWithParams(new b.BatchQueryParams(0)),e=>e.recordErrors[0].operation()==="n_jobs");
    assert.deepEqual(b.MoleculeBatch.fromSmilesList([]).dgBoundsMatrixList(),[]);assert.deepEqual(value.validMask(),[true,false,true]);
});

test("synchronous progress callbacks cover every row, survive reentry and retain original exceptions",()=>{
    const value=batch(),write=new b.SmilesWriteParams();let ticks=0;const callback=()=>{ticks++;};const params=new b.BatchQueryParams(null,false,callback);
    value.toSmilesListWithParams(write,params);assert.equal(ticks,3);ticks=0;
    value.dgBoundsMatrixListWithParams(params);assert.equal(ticks,3);ticks=0;
    value.toSvgListWithParams(120,100,params);assert.equal(ticks,3);ticks=0;
    assert.throws(()=>value.toSvgListWithParams(0,100,params));assert.equal(ticks,3);
    ticks=0;b.MoleculeBatch.fromSmilesList([]).toSmilesListWithParams(write,params);assert.equal(ticks,0);
    let depth=0,nested=0;const reentrant=new b.BatchQueryParams(null,null,()=>{ticks++;if(depth===0){depth++;nested++;assert.deepEqual(value.toSmilesListWithParams(write,reentrant),["CCO",null,"O"]);depth--;}});
    assert.deepEqual(value.toSmilesListWithParams(write,reentrant),["CCO",null,"O"]);assert.equal(nested,3);assert.equal(ticks,12);
    const first=new Error("first callback failure"),second=new Error("later callback failure");let count=0;
    const failing=new b.BatchQueryParams(null,null,()=>{count++;throw count===1?first:second;});
    assert.throws(()=>value.toSmilesListWithParams(write,failing),error=>error===first);assert.equal(count,3);
    count=0;assert.throws(()=>value.toSvgListWithParams(0,100,failing),error=>error!==first&&error.domain==="batch");assert.equal(count,3);
    const symbol=Symbol("native thrown value");assert.throws(()=>value.toSmilesListWithParams(write,new b.BatchQueryParams(null,null,()=>{throw symbol;})),e=>e===symbol);
    assert.deepEqual(value.toSmilesList(),["CCO",null,"O"]);
});
