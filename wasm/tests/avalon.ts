import {Molecule, AvalonFingerprintParams, Fingerprint, AvalonFingerprintError} from 'cosmolkit-generated';
const mol = Molecule.fromSmiles('CCO'), params = new AvalonFingerprintParams();
params.nBits = 64; params.isQuery = false; params.bitFlags = 32767;
const fp: Fingerprint = mol.fingerprintAvalon(params);
const defaultFp: Fingerprint = mol.fingerprintAvalon();
const options: Fingerprint = mol.fingerprintAvalon({nBits:64, isQuery:false});
declare const error: AvalonFingerprintError;
const cause: Error|null = error.cause;
// @ts-expect-error boolean configuration must not accept numbers
params.isQuery = 1;
