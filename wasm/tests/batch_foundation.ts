import { BatchErrorMode, BatchParams, BatchError, BatchValidationError, BatchRecord, MoleculeBatch, Molecule, SmilesParseParams } from "cosmolkit-generated";
const parse = new SmilesParseParams(true,true,true,true,false,false,false,new Map([["{R}","C"]]));
const map: Map<string,string> = parse.replacements;
const params = new BatchParams(BatchErrorMode.KeepErrors,null,false);
const jobs: number|null = params.nJobs, progress: boolean|null = params.progressBar;
// @ts-expect-error Frozen parameters match the Python value.
params.nJobs = 2;
const scalar: Molecule = Molecule.fromSmilesWithParams("{R}",parse);
const batch = MoleculeBatch.fromSmilesListWithParams(["C","["],parse,params);
const defaultBatch: MoleculeBatch = MoleculeBatch.fromSmilesList(["C"]);
const records: BatchRecord[] = [BatchRecord.molecule(scalar)];
const copied: MoleculeBatch = MoleculeBatch.fromRecords(records,BatchErrorMode.Strict);
const count: number = batch.len(), empty: boolean = batch.isEmpty(), mask: boolean[] = batch.validMask(), invalid: boolean[] = batch.invalidMask();
const valid: number = batch.validCount(), failures: number = batch.invalidCount(), errors: BatchError[] = batch.errors();
const values: (Molecule|null)[] = batch.toList();
const configured: MoleculeBatch = batch.withParallelJobs(2).withProgressBar(false).withValidRecords();
const storedJobs: number|null = batch.parallelJobs(), storedProgress: boolean|null = batch.progressBar();
const row: BatchError = errors[0], index: number = row.index(), operation: string = row.operation(), message: string = row.message();
const details: {index:number;operation:string;message:string} = row.asDict(), cause: Error|null = row.cause();
const failed: BatchRecord = BatchRecord.error(row), maybeMolecule: Molecule|null = failed.moleculeValue(), maybeError: BatchError|null = failed.errorValue();
function inspect(error: BatchValidationError) {
    const count: number = error.errors, reason: null = error.reason, records: BatchError[] = error.recordErrors, cause: Error|null = error.cause;
    return [count,reason,records,cause,error.domain,error.kind,error.message];
}
void [map,jobs,progress,defaultBatch,copied,count,empty,mask,invalid,valid,failures,values,configured,storedJobs,storedProgress,index,operation,message,details,cause,maybeMolecule,maybeError,inspect];
