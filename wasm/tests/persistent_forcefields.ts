import {Molecule, MolecularForceField, MmffForceFieldParams, UffForceFieldParams, ForceFieldMinimizeParams, ForceFieldMinimizeOutcome, ForceFieldEnergyGradient, MolecularForceFieldErrorKind, ForceFieldError} from 'cosmolkit-generated';
const molecule = Molecule.fromSmiles('CCO');
const mmff = new MmffForceFieldParams();
mmff.mmffVariant = 'MMFF94s';
mmff.conformerId = null;
const uff = new UffForceFieldParams();
uff.vdwThreshold = 10;
const handles: MolecularForceField[] = [molecule.mmffForceField(mmff), molecule.uffForceField(uff), molecule.mmffForceField({mmffVariant: 'MMFF94'}), molecule.uffForceFieldWithParams(uff)];
const params = new ForceFieldMinimizeParams();
params.maxIterations = 2;
params.energyTolerance = 1e-6;
const outcome: ForceFieldMinimizeOutcome = handles[0].minimize(params);
handles[0].minimize({maxIterations: 2});
const result: ForceFieldEnergyGradient = handles[0].energyGradient();
const gradient: number[][] = result.gradient();
handles[0].setPosition(0, new Float64Array([1, 2, 3]));
handles[0].setPositions(gradient);
handles[0].setFixedAtoms(new Uint32Array([0]));
declare const error: ForceFieldError;
const kind: MolecularForceFieldErrorKind = error.kind();
// @ts-expect-error tolerance must be a number
params.forceTolerance = 'small';
// @ts-expect-error pin indices are numeric
handles[0].setFixedAtoms(['0']);
void [outcome, kind];
