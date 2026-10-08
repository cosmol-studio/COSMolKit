import {Molecule, InchiReadParams, InchiWriteParams, inchiToInchiKey} from "cosmolkit-generated";
const read = new InchiReadParams();
read.sanitize = false;
read.removeHydrogens = true;
const write = new InchiWriteParams("-AuxNone");
write.options = "";
const mol = Molecule.fromInchi("InChI=1S/CH4/h1H4", read);
const identifier: string = mol.toInchi(write);
const key: string = inchiToInchiKey(identifier);
mol.toInchi({options: "-AuxNone"});
Molecule.fromInchi(identifier, {sanitize: false, removeHydrogens: true});
mol.toInchiWithParams(write);
mol.toInchiKeyWithParams(write);
// @ts-expect-error options must be a string
write.options = 123;
// @ts-expect-error sanitize must be boolean
read.sanitize = "false";
// @ts-expect-error removeHydrogens must be boolean
read.removeHydrogens = 1;
void key;
