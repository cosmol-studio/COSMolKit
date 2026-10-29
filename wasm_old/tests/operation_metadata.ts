import {featureSpecs, featureSpec, operationSpecs, operationSpec, supportMatrix, operationInvariantMatrix, parityMatrix, operationInvariant, operationParity, FeatureSpec, MoleculeOpSpec, FunctionStatus, ParityPolicy, SupportMatrixEntry, OperationInvariantEntry, ParityMatrixEntry} from 'cosmolkit-generated';
const feature: FeatureSpec|null=featureSpecs().next();
const found:FeatureSpec|null=featureSpec('core');
const specs:MoleculeOpSpec[]=operationSpecs();
const spec:MoleculeOpSpec|null=operationSpec('sanitize');
const status:FunctionStatus=specs[0].status;
const reference:string|null=status.reference;
const explanation:string|null=status.explanation;
const policy:ParityPolicy=specs[0].parity;
const access:{read:number;write:number}=specs[0].access;
const effects:{recompute:number;preserve:number;invalidate:number;operationDefined:number}=specs[0].derivedEffects;
const support:SupportMatrixEntry[]=supportMatrix();
const optional:MoleculeOpSpec|null=support[0].operation;
const invariant:OperationInvariantEntry[]=operationInvariantMatrix();
const parity:ParityMatrixEntry[]=parityMatrix();
const one:OperationInvariantEntry|null=operationInvariant('sanitize');
const pair:ParityMatrixEntry|null=operationParity('sanitize');
const version:string|null=parity[0].rdkitVersion;
// @ts-expect-error metadata is read-only
specs[0].mayMutate=255;
// @ts-expect-error feature names are read-only
feature!.name='changed';
