import assert from 'node:assert/strict';import {readFileSync} from 'node:fs';import test from 'node:test';import {pathToFileURL} from 'node:url';
const b=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);b.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
const raw=text=>b.Molecule.fromSmilesWithParams(text,new b.SmilesParseParams(false,undefined,undefined,undefined,false));
const outcome=fn=>{try{return {value:fn()};}catch(error){return {error};}};
test('All four sanitize operations retain stage selections, reusable parameters and failure atomicity',()=>{
 const flags=['NONE','ALL','CLEANUP','PROPERTIES','SYMM_RINGS','KEKULIZE','FIND_RADICALS','SET_AROMATICITY','SET_CONJUGATION','SET_HYBRIDIZATION','CLEANUP_CHIRALITY','ADJUST_HS','CLEANUP_ORGANOMETALLICS','CLEANUP_ATROPISOMERS'];
 for(const text of ['', 'CCO','c1ccccc1','[CH3]','C(F)(F)(F)(F)F','c','c1cccc1']){
  const original=raw(text),before=original.toSmiles();
  for(const flag of flags){const params=new b.SanitizeParams(b.SanitizeOperations[flag]),a=outcome(()=>original.sanitizeWithParams(params)),m=raw(text),i=outcome(()=>m.sanitizeWithParams_(params));
   assert.equal(Boolean(a.error),Boolean(i.error));
   if(a.error){assert.equal(a.error.name,'OperationError');assert.equal(a.error.kind,'Sanitize');assert.equal(i.error.message,a.error.message);assert.equal(i.error.cause.kind,a.error.cause.kind);assert.equal(m.toSmiles(),before);assert.ok(a.error.cause.detail instanceof b.SanitizeError);assert.equal(a.error.cause.detail.stage,a.error.cause.stage);}
   else{assert.equal(i.value,undefined);assert.equal(m.toSmiles(),a.value.toSmiles());assert.equal(original.sanitizeWithParams(params).toSmiles(),a.value.toSmiles());}
   assert.equal(original.toSmiles(),before);assert.equal(params.operations.bits(),b.SanitizeOperations[flag].bits());
  }
  const a=outcome(()=>original.sanitize()),m=raw(text),i=outcome(()=>m.sanitize_());assert.equal(Boolean(a.error),Boolean(i.error));if(a.error){assert.equal(a.error.message,i.error.message);assert.equal(m.toSmiles(),before);}else{assert.equal(m.toSmiles(),a.value.toSmiles());}
 }
 const m=raw('CCO'),out=m.sanitize();m.free();assert.equal(out.toSmiles(),'CCO');
});
test('Chemistry problems retain order, exact stages, typed valence and Kekulize causes and independent report lifetime',()=>{
 const input=raw('C(F)(F)(F)(F)F.C(F)(F)(F)(F)F'),before=input.toSmiles(),params=new b.SanitizeParams(b.SanitizeOperations.PROPERTIES),report=input.detectChemistryProblemsWithParams(params);
 assert.ok(report instanceof b.ChemistryProblemReport);assert.equal(report.problems.length,2);assert.deepEqual(report.problems.map(p=>p.error.cause.atom),[0,6]);
 for(const problem of report.problems){assert.ok(problem instanceof b.ChemistryProblem);assert.equal(problem.operation,b.SanitizeStage.Properties);const e=problem.error;assert.ok(e instanceof Error);assert.equal(e.name,'ChemistryProblemError');assert.equal(e.domain,'sanitize');assert.equal(e.kind,'Valence');assert.ok(e.detail instanceof b.ChemistryProblemError);assert.equal(e.detail.kind,'Valence');const cause=e.cause;assert.equal(cause.name,'ValenceError');assert.equal(cause.kind,'InvalidValence');assert.equal(cause.domain,'valence');assert.equal(cause.atomicNumber,6);assert.equal(cause.formalCharge,0);assert.equal(cause.phase,'Explicit');assert.equal(cause.calculated,5);assert.equal(typeof cause.reason,'string');assert.ok(cause.detail instanceof b.ValenceError);assert.equal(cause.detail.atom,cause.atom);assert.equal(cause.detail.calculated,5);assert.equal(cause.detail.order,null);assert.equal(cause.detail.cause,null);}
 const copy=report.problems;copy.length=0;assert.equal(report.problems.length,2);assert.equal(input.toSmiles(),before);assert.equal(input.detectChemistryProblemsWithParams(params).problems.length,2);input.free();assert.equal(report.problems[1].error.cause.atom,6);
 const kek=raw('c1cccc1'),kreport=kek.detectChemistryProblemsWithParams(new b.SanitizeParams(b.SanitizeOperations.KEKULIZE));assert.equal(kreport.problems.length,1);const kp=kreport.problems[0];assert.equal(kp.operation,b.SanitizeStage.Kekulize);assert.equal(kp.error.kind,'Kekulize');assert.equal(kp.error.cause.kind,'NotKekulizable');assert.ok(kp.error.cause.detail instanceof b.KekulizeError);assert.deepEqual(kp.error.cause.detail.problemAtoms,[0,1,2,3,4]);
 assert.equal(raw('CCO').detectChemistryProblems().problems.length,0);assert.equal(raw('C(F)(F)(F)(F)F').detectChemistryProblemsWithParams(new b.SanitizeParams(b.SanitizeOperations.NONE)).problems.length,0);
});
test('Sanitize parameter classes and invalid operation masks retain strict boundary errors',()=>{
 const m=raw('CCO'),before=m.toSmiles();for(const v of [{},undefined,null,new b.KekulizeParams()])for(const name of ['sanitizeWithParams','sanitizeWithParams_','detectChemistryProblemsWithParams'])assert.throws(()=>m[name](v),TypeError);
 assert.throws(()=>b.SanitizeOperations.fromBits(0x1000),e=>e.kind==='InvalidOperations'&&e.detail.bits===0x1000&&e.detail.unknownBits===0x1000&&e.detail.stage===null);assert.equal(m.toSmiles(),before);
});
