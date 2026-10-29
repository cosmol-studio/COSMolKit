import {Molecule} from 'cosmolkit-generated';
const mol = Molecule.fromSmiles('CC.O');
const parts: Molecule[] = mol.fragments();
const largest: Molecule = mol.largestFragment();
