import assert from "node:assert/strict";
import {readFileSync} from "node:fs";
import test from "node:test";
import {pathToFileURL} from "node:url";
const b=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);b.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});

test("image parameters preserve all defaults, copied filenames, frozen nested execution and checked dimensions",()=>{
    const defaults=new b.BatchImageParams();assert.equal(defaults.format,"png");assert.equal(defaults.width,300);assert.equal(defaults.height,300);assert.equal(defaults.filenames,null);assert.equal(defaults.reportPath,null);assert.equal(defaults.execution.errors,null);
    const execution=new b.BatchParams(b.BatchErrorMode.KeepErrors,2,false),filenames=["ethanol",null,"water.svg"];
    const custom=new b.BatchImageParams("svg",120,100,execution,filenames,"counts.json");filenames[0]="changed";
    assert.equal(execution.nJobs,2);assert.equal(custom.execution.nJobs,2);assert.equal(custom.execution.progressBar,false);assert.equal(custom.execution.errors,b.BatchErrorMode.KeepErrors);
    assert.equal(custom.width,120);assert.equal(custom.height,100);assert.equal(custom.format,"svg");assert.equal(custom.reportPath,"counts.json");assert.deepEqual(custom.filenames,["ethanol",null,"water.svg"]);
    const copied=custom.filenames;copied[0]="changed";assert.deepEqual(custom.filenames,["ethanol",null,"water.svg"]);
    assert.deepEqual(new b.BatchImageParams(undefined,undefined,undefined,null,[]).filenames,[]);
    assert.throws(()=>{custom.format="png";},TypeError);assert.throws(()=>new b.BatchImageParams(7),TypeError);assert.throws(()=>new b.BatchImageParams("svg",120,100,{}),TypeError);
    assert.throws(()=>new b.BatchImageParams("svg",120,100,null,[7]),TypeError);assert.throws(()=>new b.BatchImageParams("svg",120,100,null,null,7),TypeError);
    for(const n of [-1,1.5,4294967296,NaN,Infinity])assert.throws(()=>new b.BatchImageParams("png",n),RangeError);
});

test("image export preserves the current WASM filesystem error and its complete typed cause",()=>{
    const batch=b.MoleculeBatch.fromSmilesList(["CCO"]);
    assert.throws(()=>batch.writeImages(""),e=>{
        assert.equal(e.domain,"batch");assert.equal(e.kind,"Validation");assert.equal(e.errors,1);assert.equal(e.cause.index,0);assert.equal(e.cause.operation,"batch.write_images");
        const image=e.cause.cause;assert.equal(image.name,"BatchImageError");assert.equal(image.kind,"Write");assert.ok(image.detail instanceof b.BatchImageError);
        const write=image.cause;assert.equal(write.name,"DrawingWriteError");assert.equal(write.kind,"Io");assert.equal(write.filename,"mol_0.png");assert.ok(write.detail instanceof b.DrawingWriteError);
        assert.equal(write.cause.domain,"io");assert.equal(write.cause.kind,"Unsupported");assert.equal(write.cause.errno,null);return true;
    });
});

test("KeepErrors retains failed write rows and report writing preserves the native filesystem failure",()=>{
    const keep=new b.BatchParams(b.BatchErrorMode.KeepErrors);
    const batch=b.MoleculeBatch.fromSmilesListWithParams(["CCO","[","O"],new b.SmilesParseParams(),keep);
    const report=batch.writeImagesWithParams("",new b.BatchImageParams("svg",120,100,keep,["ethanol",null,"water.svg"]));
    assert.ok(report instanceof b.BatchExportReport);assert.equal(report.total(),3);assert.equal(report.written,0);assert.equal(report.success(),0);assert.equal("skipped" in report,false);assert.equal(report.failed(),3);
    const errors=report.errors();assert.deepEqual(errors.map(row=>row.index()),[0,1,2]);assert.equal(report.failed(),errors.length);
    const inputError=errors[1];assert.equal(inputError.operation(),"batch.from_smiles_list");assert.equal(inputError.cause().domain,"smiles");assert.equal(inputError.cause().kind,"Parse");
    const writeErrors=[errors[0],errors[2]];
    for(const row of writeErrors){assert.equal(row.operation(),"batch.write_images");assert.equal(row.cause().name,"BatchImageError");assert.equal(row.cause().cause.cause.kind,"Unsupported");}
    assert.deepEqual(writeErrors.map(row=>row.cause().cause.filename),["ethanol.svg","water.svg"]);
    assert.throws(()=>report.writeReport("counts.json"),e=>e.domain==="batch"&&e.cause.operation==="write error report"&&e.cause.cause.domain==="io"&&e.cause.cause.kind==="Unsupported");
    const empty=b.MoleculeBatch.fromSmilesList([]).writeImages("");assert.equal(empty.total(),0);assert.equal(empty.success(),0);assert.equal(empty.failed(),0);assert.deepEqual(empty.errors(),[]);
});

test("invalid jobs and filename lengths fail before any filesystem operation",()=>{
    const batch=b.MoleculeBatch.fromSmilesList(["C"]);
    assert.throws(()=>batch.writeImagesWithParams("",new b.BatchImageParams("png",300,300,new b.BatchParams(b.BatchErrorMode.Strict,0))),e=>e.domain==="batch"&&e.recordErrors[0].operation()==="n_jobs");
    assert.throws(()=>batch.writeImagesWithParams("",new b.BatchImageParams("png",300,300,null,[])),e=>e.domain==="batch"&&e.recordErrors[0].operation()==="filenames");
});
