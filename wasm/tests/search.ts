import {Molecule,QueryGraph,SmartsParseParams,SmartsParseError,parseSmarts,parseSmartsWithParams,CompiledQuery,MatchResult,SubstructMatchParams,SmartsWriteParams,SubstructMatchError,QueryCompileError,MatchError,SmartsWriteError,compileQuery,writeSmarts,writeCxSmarts} from 'cosmolkit-generated';
const parse=new SmartsParseParams(),q:QueryGraph=parseSmarts('CO'),q2:QueryGraph=parseSmartsWithParams('CO',parse),q3:QueryGraph=QueryGraph.fromSmarts('CO'),q4:QueryGraph=QueryGraph.fromSmartsWithParams('CO',parse),m=Molecule.fromSmiles('CCO'),params=new SubstructMatchParams(1000,true,false,false,false,false,true,1000,1,false,false,['atom'],['bond'],false,false,false),compiled:CompiledQuery=compileQuery(q),w=new SmartsWriteParams(true,true,true,null);
const first:MatchResult|null=m.substructMatch(q),all:MatchResult[]=m.substructMatches(q),has:boolean=m.hasSubstructMatch(q),explicit:MatchResult[]=m.substructMatchesWithParams(q,params),planned:MatchResult[]=m.substructMatchesCompiled(compiled),smarts:string=writeSmarts(q,w),cx:string=writeCxSmarts(q,w),mapping:number[]=all[0].atomMapping(),bonds:number[]=all[0].bondMapping(),pairs:[number,number][]=all[0].atomPairs(),order:number[]=compiled.atomOrder(),copy:QueryGraph=compiled.query(),n:number=compiled.numAtoms(),nb:number=compiled.numBonds();
const mm:number=params.maxMatches,u:boolean=params.uniquify,c:boolean=params.useChirality,e:boolean=params.useEnhancedStereo,s:boolean=params.specifiedStereoQueryMatchesUnspecified,qq:boolean=params.useQueryQueryMatches,r:boolean=params.recursionPossible,mr:number=params.maxRecursiveMatches,nt:number=params.numThreads,ac:boolean=params.aromaticMatchesConjugated,asd:boolean=params.aromaticMatchesSingleOrDouble,ap:string[]=params.atomProperties,bp:string[]=params.bondProperties,ea:boolean=params.extraAtomCheckOverridesDefaultCheck,eb:boolean=params.extraBondCheckOverridesDefaultCheck,g:boolean=params.useGenericMatchers;
const am:boolean=w.includeAtomMaps,iso:boolean=w.doIsomericSmiles,db:boolean=w.includeDativeBonds,root:number|null=w.rootedAtAtom;
declare const se:SubstructMatchError,ce:QueryCompileError,me:MatchError,we:SmartsWriteError,pe:SmartsParseError;const errors:string[]=[se.kind,ce.kind,me.kind,we.kind],cause:Error|null=se.cause,branch:string|null=se.branch,fn:string|null=se.rdkitFunction,position:number|null=pe.position,character:string|null=pe.character,context:string|null=pe.context,pd:string|null=pe.primitiveDetail,ring:number|null=pe.ring,atom:number|null=pe.atom,begin:number|null=pe.beginAtom,end:number|null=pe.endAtom,feature:string|null=pe.feature,carrier:number|null=pe.carrier;
// @ts-expect-error exact query class
m.substructMatches(compiled);
// @ts-expect-error nullable first match
const required:MatchResult=m.substructMatch(q);
// @ts-expect-error checked number thread count
new SubstructMatchParams(undefined,undefined,undefined,undefined,undefined,undefined,undefined,undefined,1n);
// @ts-expect-error immutable policies
params.numThreads=2;
// @ts-expect-error writer params required
writeSmarts(q);
