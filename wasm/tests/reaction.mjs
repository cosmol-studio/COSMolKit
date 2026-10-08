import assert from "node:assert/strict";
import {readFileSync} from "node:fs";
import test from "node:test";
import {pathToFileURL} from "node:url";
const b=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
b.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
const molecule=text=>b.Molecule.fromSmilesWithSanitize(text,true);
const reaction=()=>b.Reaction.fromSmirks("[C:1]>>[N:1]");

test("Reaction values preserve ordered templates, settings, initialization and text",()=>{
 const rxn=b.Reaction.fromSmirks("[C:1]>O>[N:1]");
 assert.deepEqual([rxn.numReactantTemplates(),rxn.numProductTemplates(),rxn.numAgentTemplates()],[1,1,1]);
 assert.equal(rxn.reactantTemplate(0).numAtoms(),1);
 assert.equal(rxn.reactantTemplates().length,1);
 assert.equal(rxn.productTemplates().length,1);
 assert.equal(rxn.agentTemplates().length,1);
 assert.equal(rxn.isInitialized(),false);
 const params=new b.ReactionValidationParams(true),report=rxn.validateWithParams(params);
 assert.equal(report.isValid,true);assert.equal(report.numErrors,0);assert.deepEqual(report.errors,[]);
 assert.equal(rxn.withInitialized().isInitialized(),true);
 assert.equal(rxn.withInitializedWithParams(params).isInitialized(),true);
 assert.equal(rxn.isInitialized(),false);
 for(const text of [rxn.toSmirks(),rxn.toSmirksWithParams(new b.ReactionWriteParams()),
 rxn.toCxSmirks(),rxn.toCxSmirksWithParams(new b.ReactionWriteParams())]){
   assert.equal(typeof text,"string");assert.equal(b.parseSmirks(text).numAgentTemplates(),1);
 }
 const query=b.QueryGraph.fromSmarts("[C:1]");
 const built=b.Reaction.fromTemplates([query],[query],[]);
 assert.equal(built.numAgentTemplates(),0);
 const extended=built.withReactantTemplate(query).withProductTemplate(query).withAgentTemplate(query);
 assert.deepEqual([extended.numReactantTemplates(),extended.numProductTemplates(),extended.numAgentTemplates()],[2,2,1]);
 assert.equal(extended.withoutAgents().removedTemplates.length,1);
 assert.equal(extended.withoutAgents().reaction.numAgentTemplates(),0);
 assert.equal(extended.numAgentTemplates(),1);
 // Reaction.h:376: manual ChemicalReaction defaults the flag to false.
 assert.equal(built.withImplicitProperties(true).implicitProperties(),true);
 assert.equal(built.implicitProperties(),false);
 assert.equal(built.withMatchParams(new b.SubstructMatchParams(7)).matchParams().maxMatches,7);
 for(const name of ["withoutUnmappedReactants","withoutUnmappedProducts"]){
   assert.ok(built[name]() instanceof b.ReactionTemplateRemoval);
   assert.ok(built[name+"WithParams"](new b.ReactionTemplateRemovalParams()) instanceof b.ReactionTemplateRemoval);
 }
 assert.throws(()=>b.Reaction.fromTemplates([{}],[query],[]),TypeError);
});

test("Reaction options retain typed coordinate IDs and refuse host coercion",()=>{
 const parse=new b.ReactionParseParams();
 assert.deepEqual([parse.useSmiles,parse.sanitize,parse.allowCxsmiles,parse.strictCxsmiles],[false,false,true,true]);
 assert.equal(parse.replacements.size,0);
 const params=new b.ReactionRunParams();
 assert.equal(params.maxProducts,1000);assert.deepEqual(params.coordinateSelections,[]);
 assert.equal(new b.ReactionApplyParams().removeUnmatchedAtoms,true);
 assert.equal(new b.ReactionTemplateRemovalParams().thresholdUnmappedAtoms,0.2);
 const auto=b.ReactionCoordinateSelection.auto(),two=b.ReactionCoordinateSelection.twoD(7),three=b.ReactionCoordinateSelection.threeD(42);
 assert.equal(auto.isAuto,true);assert.equal(auto.id,null);
 assert.equal(two.is2d,true);assert.equal(two.is3d,false);assert.equal(two.id,7);
 assert.equal(three.is3d,true);assert.equal(three.id,42);
 assert.equal(new b.ReactionSingleRunParams(three).coordinateSelection.id,42);
 assert.deepEqual(new b.ReactionRunParams(3,[two,three]).coordinateSelections.map(v=>v.id),[7,42]);
 const write=new b.ReactionWriteParams(undefined,undefined,0,undefined,undefined,b.CxSmilesFields.NONE,[two]);
 assert.equal(write.cxFields.bits(),0);assert.equal(write.rootedAtAtom,0);assert.equal(write.coordinateSelections[0].id,7);
 assert.equal(b.CxSmilesFields.COORDS.or(b.CxSmilesFields.ATOM_LABELS).bits(),5);
 assert.equal(b.CxSmilesFields.ALL.contains(b.CxSmilesFields.COORDS),true);
 assert.throws(()=>{params.maxProducts=4;},TypeError);
 assert.throws(()=>new b.ReactionSingleRunParams("Auto"),TypeError);
 assert.throws(()=>new b.ReactionParseParams("false"),TypeError);
 assert.throws(()=>new b.ReactionRunParams(-1),RangeError);
 assert.throws(()=>b.ReactionCoordinateSelection.threeD(1.5),RangeError);
 const replacements=new b.ReactionParseParams(undefined,undefined,new Map([["{C}","[C:1]"]]));
 assert.equal(b.parseSmirksWithParams("{C}>>[N:1]",replacements).numReactantTemplates(),1);
});

test("Product groups preserve sources and initialize the same Reaction handle",()=>{
 for(const explicit of [false,true]){
 const source=molecule("C"),rxn=reaction(),original=source.toSmiles();
 assert.equal(rxn.isInitialized(),false);
 const sets=explicit?source.reactionProductsWithParams(rxn,0,new b.ReactionSingleRunParams()):source.reactionProducts(rxn,0);
 assert.equal(rxn.isInitialized(),true);assert.equal(sets.length,1);assert.equal(sets[0].length,1);
 assert.equal(sets[0][0].toSmiles(),"N");assert.equal(source.toSmiles(),original);
 assert.deepEqual(molecule("O").reactionProducts(rxn,0),[]);
 const pair=b.Reaction.fromSmirks("[C:1].[C:2]>>[C:1].[C:2]");
 const pairs=source.reactionProductsFromInputs(pair,[source,source],new b.ReactionRunParams());
 assert.deepEqual(pairs.map(set=>set.map(m=>m.toSmiles())),[["C","C"]]);
 assert.equal(source.toSmiles(),original);
 assert.throws(()=>source.reactionProductsFromInputs(pair,[source],new b.ReactionRunParams()),e=>e.name==="OperationError"&&e.kind==="ReactionRun"&&e.cause.name==="ReactionRunError"&&e.cause.kind==="ReactantArity"&&e.cause.expected===2&&e.cause.actual===1&&e.cause.detail instanceof b.ReactionRunError);
 }
});

test("Apply value/inplace results retain bool semantics, COW and typed failure atomicity",()=>{
 for(const explicit of [false,true]){
 const source=molecule("C"),rxn=reaction();
 const result=explicit?source.applyReactionWithParams(rxn,new b.ReactionApplyParams()):source.applyReaction(rxn);
 assert.equal(result.changed,true);assert.equal(result.molecule.toSmiles(),"N");assert.equal(source.toSmiles(),"C");
 assert.equal(explicit?source.applyReactionWithParams_(rxn,new b.ReactionApplyParams()):source.applyReaction_(rxn),true);
 assert.equal(source.toSmiles(),"N");assert.equal(source.applyReaction_(rxn),false);
 assert.throws(()=>source.applyReaction_(b.Reaction.fromSmirks("[N:1]>>[N:1]C")),e=>e.name==="OperationError"&&e.kind==="ReactionApply"&&e.cause.name==="ReactionApplyError"&&e.cause.kind==="AddsProductAtom"&&e.cause.atom===1&&e.cause.detail instanceof b.ReactionApplyError);
 assert.equal(source.toSmiles(),"N");
 }
});

test("Reaction failures preserve error kinds, contextual fields and validation reports",()=>{
 assert.throws(()=>b.Reaction.fromSmirks("[C:1]>>["),e=>e.name==="ReactionParseError"&&e.kind==="Smarts"&&e.role===b.ReactionRole.Product&&e.template===0&&e.text==="["&&e.cause.name==="SmartsParseError"&&e.cause.kind==="UnclosedBracket");
 assert.throws(()=>b.Reaction.fromSmirks("CC"),e=>e.name==="ReactionParseError"&&e.kind==="Separators"&&e.count===0&&e.detail instanceof b.ReactionParseError&&e.detail.count===0);
 const rxn=reaction();
 assert.throws(()=>rxn.reactantTemplate(9),e=>e.name==="ReactionModelError"&&e.kind==="TemplateIndex"&&e.index===9&&e.count===1&&e.role===b.ReactionRole.Reactant);
 assert.throws(()=>molecule("C").reactionProducts(rxn,99),e=>e.name==="OperationError"&&e.kind==="ReactionRun"&&e.cause.kind==="ReactantTemplateIndex"&&e.cause.index===99);
 const empty=new b.Reaction(),report=empty.validateWithParams(new b.ReactionValidationParams(true));
 assert.equal(report.isValid,false);assert.equal(report.numErrors,2);
 assert.equal(report.errors[0].kind,b.ReactionValidationIssueKind.MissingReactants);
 assert.equal(report.errors[0].severity,b.ReactionValidationSeverity.Error);
 assert.equal(report.errors[0].role,b.ReactionRole.Reactant);assert.equal(report.errors[0].atom,null);
 assert.equal(typeof report.errors[0].detail,"string");
 assert.throws(()=>empty.withInitializedWithParams(new b.ReactionValidationParams(true)),e=>e.name==="ReactionInitializationError"&&e.kind==="Invalid"&&e.report.numErrors===2&&e.detail.report instanceof b.ReactionValidationReport);
});
