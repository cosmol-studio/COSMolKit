import { Molecule, ValenceModel, ValenceParams } from "cosmolkit-generated";

const params = new ValenceParams(ValenceModel.RdkitLike, false);
const model: ValenceModel = params.model;
const strict: boolean = params.strict;
declare const molecule: Molecule;
const prepared: Molecule = molecule.withAssignedValenceWithParams(params);
const assigned: void = molecule.assignValence_();
const assignedWithParams: void = molecule.assignValenceWithParams_(params);
// @ts-expect-error Parameter values are not strings.
new ValenceParams("RdkitLike");
// @ts-expect-error A preparation mutation does not return a molecule.
const invalid: Molecule = molecule.assignValence_();
// @ts-expect-error Explicit parameters are required.
molecule.withAssignedValenceWithParams();
params.strict = true;
void [model, strict, prepared, assigned, assignedWithParams];
