import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import test from 'node:test';
import {pathToFileURL} from 'node:url';
const ck=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
ck.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
function collect(iterator){const values=[];for(let row=iterator.next();row!==null;row=iterator.next())values.push(row.toSmiles());return values;}

test('Lazy stereoisomers retain their source, exact bigint count and yielded state',()=>{
 const mol=ck.Molecule.fromSmiles('CC(O)Cl');assert.equal(mol.stereoisomerCount(),2n);
 const iterator=mol.enumerateStereoisomers({maxIsomers:0});assert.equal(iterator.yieldedCount,0);
 const before=mol.toSmiles();mol.free();
 assert.deepEqual(collect(iterator).sort(),['C[C@@H](O)Cl','C[C@H](O)Cl'].sort());
 assert.equal(iterator.yieldedCount,2);assert.equal(iterator.next(),null);assert.equal(iterator.next(),null);
 assert.equal(before,'CC(O)Cl');
 const assigned=ck.Molecule.fromSmiles('C[C@H](O)Cl');
 assert.equal(assigned.stereoisomerCount({onlyUnassigned:true}),1n);
 assert.equal(assigned.stereoisomerCount({onlyUnassigned:false}),2n);
});

test('Enumeration options preserve arbitrary-width seeds and provider identity',()=>{
 const mol=ck.Molecule.fromSmiles('CC(O)C(F)Cl');
 const source=ck.StereoisomerRandomSource.fromIntegerSeed(-(1n<<90n)+17n);
 const options=new ck.StereoisomerOptions(false,true,1,source,true,false);
 options.unique=false;options.setMaxIsomers(1);
 const first=collect(mol.enumerateStereoisomersWithOptions(options));
 assert.deepEqual(collect(mol.enumerateStereoisomers(options)),first);
 assert.equal(first.length,1);assert.ok(options.randomSource instanceof ck.StereoisomerRandomSource);
 let widths=[];options.randomSource=ck.StereoisomerRandomSource.fromRandomBits(width=>{widths.push(width);return 0n;});
 assert.equal(collect(mol.enumerateStereoisomers(options)).length,1);
 assert.equal(widths.length,1);assert.equal(widths[0],2);
 assert.equal(collect(mol.enumerateStereoisomersWithRandomBits(options,()=>1n)).length,1);
 options.setRandomSource(null);assert.equal(options.randomSource,null);
 assert.throws(()=>ck.StereoisomerRandomSource.fromIntegerSeed(42),TypeError);
 assert.throws(()=>{options.maxIsomers=1.5;},RangeError);assert.equal(options.maxIsomers,1);
});

test('Random provider failures remain original errors and fuse the lazy iterator',()=>{
 const mol=ck.Molecule.fromSmiles('CC(O)C(F)Cl'),options=new ck.StereoisomerOptions();options.maxIsomers=1;
 const original=new Error('random provider failed');
 const iterator=mol.enumerateStereoisomersWithRandomBits(options,()=>{throw original;});
 assert.throws(()=>iterator.next(),e=>e===original);assert.equal(iterator.next(),null);
 assert.equal(iterator.yieldedCount,0);
 const negative=mol.enumerateStereoisomersWithRandomBits(options,()=>-1n);
 assert.throws(()=>negative.next(),RangeError);assert.equal(negative.next(),null);
 const wrong=mol.enumerateStereoisomersWithRandomBits(options,()=>1);
 assert.throws(()=>wrong.next(),TypeError);
 assert.equal(mol.toSmiles(),'CC(O)C(F)Cl');
});
