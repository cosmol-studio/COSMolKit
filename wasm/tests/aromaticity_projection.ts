import { AromaticityModel, AromaticityParams, AromaticityError, Molecule } from "cosmolkit-generated";
const params = new AromaticityParams(AromaticityModel.Rdkit);
const defaultParams = new AromaticityParams();
const model: AromaticityModel = params.model;
// @ts-expect-error Parameters are frozen, matching the Python projection.
params.model = AromaticityModel.Simple;
const source = Molecule.fromSmiles("c1ccccc1");
const value: Molecule = source.withAssignedAromaticity();
const selected: Molecule = source.withAssignedAromaticityWithParams(params);
const inPlace: void = source.assignAromaticity();
const inPlaceSelected: void = source.assignAromaticityWithParams(defaultParams);
function inspect(error: AromaticityError) {
    const domain: string = error.domain, kind: string = error.kind, message: string = error.message;
    return [domain, kind, message, error.cause];
}
void [model, value, selected, inPlace, inPlaceSelected, inspect];
