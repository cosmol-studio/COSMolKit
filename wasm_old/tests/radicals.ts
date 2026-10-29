import {Molecule} from 'cosmolkit-generated';
const m=Molecule.fromSmiles('[C]'),value:Molecule=m.withAssignedRadicals(),mutation:void=m.assignRadicals_();
// @ts-expect-error mutation returns void
const incorrect:Molecule=m.assignRadicals_();
// @ts-expect-error no assignment options
m.withAssignedRadicals({});
