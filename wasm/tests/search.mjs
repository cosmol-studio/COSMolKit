import assert from 'node:assert/strict';import {readFileSync} from 'node:fs';import test from 'node:test';import {pathToFileURL} from 'node:url';
const b=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);b.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
const maps=rows=>rows.map(row=>row.atomMapping());
test('Short matching methods accept checked parameter instances and options without changing targets',()=>{
 const m=b.Molecule.fromSmiles('F[C@](Cl)(Br)I'),q=b.parseSmarts('F[C@@](Cl)(Br)I'),before=m.toSmiles(),params=new b.SubstructMatchParams(undefined,undefined,true);
 for(const configuration of [params,{useChirality:true}]){
  assert.equal(m.substructMatch(q,configuration),null);
  assert.equal(m.hasSubstructMatch(q,configuration),false);
  assert.deepEqual(m.substructMatches(q,configuration),[]);
 }
 assert.equal(m.hasSubstructMatch(q),true);
 assert.notEqual(m.substructMatch(q),null);
 const chain=b.Molecule.fromSmiles('CCC'),carbon=b.parseSmarts('C');
 for(const configuration of [new b.SubstructMatchParams(1),{maxMatches:1}])assert.equal(chain.substructMatches(carbon,configuration).length,1);
 for(const method of ['substructMatch','substructMatches','hasSubstructMatch']){
  assert.throws(()=>m[method](q,{unknownOption:true}),TypeError);
  assert.throws(()=>m[method](q,{useChirality:'true'}),TypeError);
  assert.throws(()=>m[method](q,{maxMatches:-1}),RangeError);
  for(const invalid of [null,[],true,'params',new Date()])assert.throws(()=>m[method](q,invalid),TypeError);
 }
 assert.equal(m.toSmiles(),before);
});
test('All five substructure methods preserve ordered mappings, compiled query lifetime and read-only targets',()=>{
 const m=b.Molecule.fromSmiles('CCOCCO'),before=m.toSmiles(),query=b.parseSmarts('CO'),compiled=b.compileQuery(query),params=new b.SubstructMatchParams();
 const expected=[[1,2],[3,2],[4,5]],rows=m.substructMatches(query);assert.deepEqual(maps(rows),expected);assert.deepEqual(m.substructMatch(query).atomMapping(),expected[0]);assert.equal(m.hasSubstructMatch(query),true);assert.deepEqual(maps(m.substructMatchesWithParams(query,params)),expected);assert.deepEqual(maps(m.substructMatchesCompiled(compiled)),expected);
 for(const row of rows){assert.ok(row instanceof b.MatchResult);assert.equal(row.bondMapping().length,1);assert.deepEqual(row.atomPairs(),row.atomMapping().map((v,i)=>[i,v]));const copy=row.atomMapping();copy[0]=99;assert.notEqual(row.atomMapping()[0],99);}
 const no=b.parseSmarts('N');assert.equal(m.substructMatch(no),null);assert.deepEqual(m.substructMatches(no),[]);assert.equal(m.hasSubstructMatch(no),false);
 assert.equal(compiled.numAtoms(),2);assert.equal(compiled.numBonds(),1);assert.deepEqual([...compiled.atomOrder()].sort(),[0,1]);const order=compiled.atomOrder();order[0]=99;assert.notEqual(compiled.atomOrder()[0],99);const copied=compiled.query();query.free();assert.deepEqual(maps(m.substructMatches(copied)),expected);assert.deepEqual(maps(m.substructMatchesCompiled(compiled)),expected);
 assert.equal(m.toSmiles(),before);m.free();assert.deepEqual(rows[0].atomMapping(),[1,2]);assert.equal(compiled.query().numAtoms(),2);
});
test('All sixteen matching options preserve defaults, copied collections, integer widths and actual policy behavior',()=>{
 const p=new b.SubstructMatchParams();const names=['maxMatches','uniquify','useChirality','useEnhancedStereo','specifiedStereoQueryMatchesUnspecified','useQueryQueryMatches','recursionPossible','maxRecursiveMatches','numThreads','aromaticMatchesConjugated','aromaticMatchesSingleOrDouble','atomProperties','bondProperties','extraAtomCheckOverridesDefaultCheck','extraBondCheckOverridesDefaultCheck','useGenericMatchers'];const defaults=[1000,true,false,false,false,false,true,1000,1,false,false,[],[],false,false,false];assert.deepEqual(names.map(n=>p[n]),defaults);
 const atoms=['a'],bonds=['b'],args=[2,false,true,true,true,true,false,3,-1,true,true,atoms,bonds,true,true,true],custom=new b.SubstructMatchParams(...args);atoms[0]='changed';bonds.push('changed');assert.deepEqual(custom.atomProperties,['a']);assert.deepEqual(custom.bondProperties,['b']);const returned=custom.atomProperties;returned.push('bad');assert.deepEqual(custom.atomProperties,['a']);assert.deepEqual(names.map(n=>custom[n]),[2,false,true,true,true,true,false,3,-1,true,true,['a'],['b'],true,true,true]);
 const m=b.Molecule.fromSmiles('CCOCCO'),q=b.parseSmarts('CO');assert.deepEqual(maps(m.substructMatchesWithParams(q,new b.SubstructMatchParams(1))),[[1,2]]);assert.deepEqual(maps(m.substructMatchesWithParams(q,new b.SubstructMatchParams(2))),[[1,2],[3,2]]);
 const ethane=b.Molecule.fromSmiles('CC'),eq=b.parseSmarts('CC');assert.equal(ethane.substructMatches(eq).length,1);assert.equal(ethane.substructMatchesWithParams(eq,new b.SubstructMatchParams(1000,false)).length,2);
 const recursive=b.parseSmarts('[$(CO)]'),r=b.Molecule.fromSmiles('CCO');assert.deepEqual(maps(r.substructMatches(recursive)),[[1]]);assert.deepEqual(r.substructMatchesWithParams(recursive,new b.SubstructMatchParams(undefined,undefined,undefined,undefined,undefined,undefined,false)),[]);
 for(const index of [0,7])for(const value of [-1,1.5,4294967296,NaN,Infinity]){const a=[];a[index]=value;assert.throws(()=>new b.SubstructMatchParams(...a),RangeError);}
 for(const value of [-2147483649,2147483648,1.5]){const a=[];a[8]=value;assert.throws(()=>new b.SubstructMatchParams(...a),RangeError);}
 const max=new b.SubstructMatchParams(4294967295,undefined,undefined,undefined,undefined,undefined,undefined,4294967295,-2147483648);assert.equal(max.numThreads,-2147483648);assert.equal(max.maxRecursiveMatches,4294967295);
 const bad=[];bad[11]=[1];assert.throws(()=>new b.SubstructMatchParams(...bad),TypeError);assert.throws(()=>new b.SubstructMatchParams(undefined,'false'),TypeError);
 for(const value of [{},null,undefined,new b.SmartsWriteParams()]){assert.throws(()=>m.substructMatchesWithParams(q,value),TypeError);assert.throws(()=>m.substructMatches(value),TypeError);assert.throws(()=>m.substructMatchesCompiled(value),TypeError);}
});
test('SMARTS writers preserve full options and source errors while shared parse errors retain typed fields',()=>{
 const defaults=new b.SmartsWriteParams();assert.deepEqual([defaults.includeAtomMaps,defaults.isomericSmiles,defaults.includeDativeBonds,defaults.rootedAtAtom],[true,true,true,null]);
 const query=b.QueryGraph.fromSmarts('[C:7]O');for(const maps of [false,true])for(const stereo of [false,true])for(const dative of [false,true])for(const root of [null,0,1]){
  const p=new b.SmartsWriteParams(maps,stereo,dative,root),text=b.writeSmarts(query,p);assert.equal(b.parseSmarts(text).numAtoms(),2);assert.equal(text.includes(':7'),maps);assert.equal(b.writeSmarts(query,p),text);assert.equal(b.parseSmarts(b.writeCxSmarts(query,p)).numBonds(),1);
 }
 assert.throws(()=>b.writeSmarts(query,new b.SmartsWriteParams(true,true,true,99)),e=>e.name==='SmartsWriteError'&&e.domain==='search'&&e.kind==='RootedAtomOutOfRange'&&e.detail instanceof b.SmartsWriteError);
 assert.throws(()=>b.parseSmarts('['),e=>e.name==='SmartsParseError'&&e.detail.position===0&&e.detail.character===null&&e.detail.feature===null);
 assert.throws(()=>b.writeSmarts(query,{}),TypeError);assert.throws(()=>b.compileQuery({}),TypeError);assert.throws(()=>new b.SmartsWriteParams('true'),TypeError);for(const v of [-1,1.5,4294967296])assert.throws(()=>new b.SmartsWriteParams(true,true,true,v),RangeError);
 for(const name of ['SubstructMatchError','MatchError','QueryCompileError','SmartsWriteError'])assert.equal(typeof Object.getOwnPropertyDescriptor(b[name].prototype,'kind').get,'function');
});
