import {AddHsParams, BatchErrorMode, BatchParams, Molecule, SmilesParseParams} from 'cosmolkit-generated';

const parser = new SmilesParseParams();
parser.sanitize = false;
parser.parseName = true;
// @ts-expect-error Configuration booleans do not accept strings.
parser.sanitize = 'false';
const hydrogen = new AddHsParams();
hydrogen.onlyOnAtoms = [2];
hydrogen.addCoords = false;
const batch = new BatchParams();
batch.nJobs = 1;
batch.errors = BatchErrorMode.KeepErrors;
batch.errors = null;
const molecule = Molecule.fromSmiles('CCO', {sanitize: true});
molecule.withHydrogens(hydrogen);
molecule.withHydrogens({onlyOnAtoms: [2]});
molecule.toSmiles({isomericSmiles: false});
