import * as ck from 'cosmolkit-generated';
declare const text: string;
const params = new ck.BioMoleculeParams(true, false, ck.PdbReadFlags.NONE, false);
const pdb: ck.Molecule = ck.molFromPdb(text, {removeHs: false});
const cif: ck.Molecule = ck.molFromMmcif(text, params);
const configured: ck.Molecule = ck.Molecule.fromMmcifWithParams(text, params);
const plain: ck.Molecule = ck.Molecule.fromPdb(text);
void [pdb, cif, configured, plain];
