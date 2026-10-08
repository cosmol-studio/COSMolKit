import {BioStructure,BioStructureParts,BioCoordinateFormat} from 'cosmolkit-generated';
const s=BioStructure.fromPdb('END'),parts:BioStructureParts=s.intoParts(),format:BioCoordinateFormat=parts.inputFormat,restored:BioStructure=BioStructure.fromParts(parts);const ok:void=BioStructure.validateParts(parts),valid:void=restored.validate();
// @ts-expect-error readonly detached format
parts.inputFormat=BioCoordinateFormat.Unknown;
// @ts-expect-error no arbitrary unvalidated object accepted as parts
BioStructure.fromParts({inputFormat:BioCoordinateFormat.Pdb});
void [format,ok,valid];
