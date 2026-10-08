import {QueryGraph,SmartsParseParams,SmartsParseError,parseSmarts,parseSmartsWithParams} from 'cosmolkit-generated';
const p=new SmartsParseParams(true,true,true,false,false,false,new Map([['{X}','C']]));const qs:QueryGraph[]=[QueryGraph.fromSmarts('C'),QueryGraph.fromSmartsWithParams('C',p),parseSmarts('C'),parseSmartsWithParams('C',p)];const name:string|null=qs[0].name(),size:number=qs[0].numAtoms(),bonds:number=qs[0].numBonds(),replacements:Map<string,string>=p.replacements;declare const error:SmartsParseError;const kind:string=error.kind;
// @ts-expect-error readonly parameters
p.mergeHs=true;
// @ts-expect-error checked replacement values
new SmartsParseParams(true,true,true,false,false,false,{X:3});
void [name,size,bonds,replacements,kind];
