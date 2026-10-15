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
const outcome: ForceFieldMinimizeOutcome = handles[0].minimize_(params);
handles[0].minimize_({maxIterations: 2});
const result: ForceFieldEnergyGradient = handles[0].energyGradient();
const gradient: number[][] = result.gradient();
handles[0].setPosition_(0, new Float64Array([1, 2, 3]));
handles[0].setPositions_(gradient);
handles[0].setFixedAtoms_(new Uint32Array([0]));
declare const error: ForceFieldError;
const kind: MolecularForceFieldErrorKind = error.kind();
// @ts-expect-error tolerance must be a number
params.forceTolerance = 'small';
// @ts-expect-error pin indices are numeric
handles[0].setFixedAtoms_(['0']);
void [outcome, kind];
import {ForceFieldDistanceConstraintParams, ForceFieldAngleConstraintParams, ForceFieldTorsionConstraintParams, ForceFieldPositionConstraintParams} from 'cosmolkit-generated';
const distance = new ForceFieldDistanceConstraintParams(false, 1, 2, 10);
distance.forceConstant = 20;
handles[0].addDistanceConstraint_(0, 1, distance);
handles[0].addDistanceConstraint_(0, 1, {relative: false, minDistance: 1, maxDistance: 2, forceConstant: 10});
handles[0].addAngleConstraint_(0, 1, 2, new ForceFieldAngleConstraintParams(false, 80, 90, 10));
handles[0].addAngleConstraint_(0, 1, 2, {relative: false, minAngleDeg: 80, maxAngleDeg: 90, forceConstant: 10});
handles[0].addTorsionConstraint_(0, 1, 2, 3, new ForceFieldTorsionConstraintParams(false, -10, 10, 10));
handles[0].addTorsionConstraint_(0, 1, 2, 3, {relative: false, minDihedralDeg: -10, maxDihedralDeg: 10, forceConstant: 10});
handles[0].addPositionConstraint_(0, new ForceFieldPositionConstraintParams(0.1, 10));
handles[0].addPositionConstraint_(0, {maxDisplacement: 0.1, forceConstant: 10});
// @ts-expect-error required constructor fields cannot be omitted
new ForceFieldDistanceConstraintParams();
// @ts-expect-error forceConstant is required, not an invented default
handles[0].addDistanceConstraint_(0, 1, {relative: false, minDistance: 1, maxDistance: 2});
// @ts-expect-error unknown configuration fields are rejected
handles[0].addPositionConstraint_(0, {maxDisplacement: 0.1, forceConstant: 10, unknown: true});
