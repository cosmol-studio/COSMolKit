import {Molecule,QueryGraph,Reaction,ReactionApplyResult,ReactionRunParams,ReactionSingleRunParams,
ReactionApplyParams,ReactionWriteParams,ReactionParseParams,ReactionValidationParams,
ReactionCoordinateSelection,ReactionValidationReport,ReactionRole,ReactionParseError,parseSmirks,
parseSmirksWithParams,CxSmilesFields} from "cosmolkit-generated";
declare const molecule:Molecule;
const rxn:Reaction=parseSmirks("[C:1]>>[N:1]");
const parsed:Reaction=parseSmirksWithParams("[C:1]>>[N:1]",new ReactionParseParams());
const selected=new ReactionSingleRunParams(ReactionCoordinateSelection.threeD(42));
const id:number|null=selected.coordinateSelection.id;
const sets:Molecule[][]=molecule.reactionProductsWithParams(rxn,0,selected);
const multi:Molecule[][]=molecule.reactionProductsFromInputs(rxn,[molecule,molecule],new ReactionRunParams(3));
const runProducts:Molecule[][]=rxn.run([molecule],new ReactionRunParams());
// @ts-expect-error Reactants must be molecule objects.
rxn.run(["C"],new ReactionRunParams());
const result:ReactionApplyResult=molecule.applyReaction(rxn);
const value:Molecule=result.molecule;
const changed:boolean=molecule.applyReactionWithParams_(rxn,new ReactionApplyParams());
const report:ReactionValidationReport=rxn.validateWithParams(new ReactionValidationParams(true));
const role:ReactionRole|null=report.errors[0].role;
const templates:QueryGraph[]=rxn.reactantTemplates();
const text:string=rxn.toCxSmirksWithParams(new ReactionWriteParams(undefined,undefined,null,undefined,undefined,CxSmilesFields.NONE));
declare const error:ReactionParseError;
const count:number|null=error.count;
// @ts-expect-error A parameter object is required, not a string.
molecule.reactionProductsWithParams(rxn,0,"default");
// @ts-expect-error A typed coordinate selector is required.
new ReactionSingleRunParams("Auto");
// @ts-expect-error Value apply returns a report, not a bare molecule.
const wrong:Molecule=molecule.applyReaction(rxn);
// @ts-expect-error Immutable result properties cannot be assigned.
result.changed=false;
// @ts-expect-error Registered options are immutable.
new ReactionRunParams().maxProducts=1;
// @ts-expect-error Coordinate IDs require numbers, never strings.
ReactionCoordinateSelection.threeD("42");
// @ts-expect-error Template indices require numbers.
rxn.reactantTemplate("0");
void [parsed,id,sets,multi,runProducts,value,changed,report,role,templates,text,count,wrong];
