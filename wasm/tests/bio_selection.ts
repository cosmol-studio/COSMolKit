import {BioStructure,BioSelection,BioSelectionParseError,BioOperationError,BioSelectionMatchError,BioSelectionCopyError,BioSelectionCopyCause,ProteinProjectionError} from 'cosmolkit-generated';
const s=BioStructure.fromPdb('END'),selection:BioSelection=BioSelection.fromCid('/'),cid:string=selection.toCid(),ids:Uint32Array=s.selectedAtomIds(selection),value:BioStructure=s.withSelection(selection),moved:BioStructure=s.withTranslatedCoordinates(new Float64Array([1,2,3]));const retained:void=s.retainSelection_(selection),translated:void=s.translate_([1,2,3]);declare const e:BioSelectionParseError;const offset:number|null=e.pos,variant:string=e.variant,info:string|null=e.info;declare const copy:BioSelectionCopyError;const cause:BioSelectionCopyCause=copy.cause();declare const projection:ProteinProjectionError;const index:number|null=projection.index;declare const op:BioOperationError;declare const match:BioSelectionMatchError;const source:Error|null=op.cause;
// @ts-expect-error coordinate argument type
s.translate_('1,2,3');
// @ts-expect-error nullable parse offset
const bad:number=e.pos;
void [cid,ids,value,moved,retained,translated,offset,variant,info,cause,index,source,match,bad];
