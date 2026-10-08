import { Molecule, ValenceModel, ValenceParams } from "cosmolkit-generated";

const params = new ValenceParams(ValenceModel.RdkitLike, false);
const model: ValenceModel = params.model;
const strict: boolean = params.strict;
declare const molecule: Molecule;
const prepared: Molecule = molecule.withAssignedValenceWithParams(params);
const assigned: void = molecule.assignValence();
const assignedWithParams: void = molecule.assignValenceWithParams(params);
// @ts-expect-error Parameter values are not strings.
new ValenceParams("RdkitLike");
// @ts-expect-error A preparation mutation does not return a molecule.
const invalid: Molecule = molecule.assignValence();
// @ts-expect-error Explicit parameters are required.
molecule.withAssignedValenceWithParams();
// @ts-expect-error Frozen parameter fields are read-only.
params.strict = true;
void [model, strict, prepared, assigned, assignedWithParams];
