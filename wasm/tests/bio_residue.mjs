import assert from 'node:assert/strict';import {readFileSync} from 'node:fs';import test from 'node:test';import {pathToFileURL} from 'node:url';
const b=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);b.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
test('BIO residue vocabulary preserves all source codes, kind names and table bounds',()=>{
 const kinds=['UNKNOWN','AA','AAD','PAA','MAA','RNA','DNA','BUF','HOH','PYR','KET','ELS'];for(const [i,name]of kinds.entries()){const kind=b.ResidueInfoKind[name];assert.ok(kind instanceof b.ResidueInfoKind);assert.equal(kind.value,i);assert.equal(kind.name(),name);}
 assert.equal(b.ResidueCode.ALA,0);assert.equal(b.ResidueCode.MSE,17);assert.equal(b.ResidueCode.UNKNOWN,367);
 for(let i=0;i<368;i++){const info=b.residueInfo(i),checked=b.residueInfoChecked(i);assert.ok(info instanceof b.ResidueInfo);assert.equal(info.code(),i);assert.equal(checked.code(),i);assert.equal(typeof b.ResidueCode[i],'string');assert.equal(info.kind().name(),info.kindName());}
 assert.equal(b.residueInfoChecked(368),null);assert.throws(()=>b.residueInfo(368),RangeError);for(const i of [-1,1.5,4294967296])assert.throws(()=>b.residueInfoChecked(i),RangeError);
 for(const name of ['ALA','MSE','HOH','not-a-residue']){const idx=b.findResidueInfoIndex(name),info=b.findResidueInfo(name);assert.equal(b.residueInfo(idx).code(),info.code());assert.equal(b.residueCode(name),info.code());}
});
test('BIO residue information and raw identities retain every Python-level field and predicate',()=>{
 const a=b.findResidueInfo('ALA');assert.deepEqual([a.code(),a.name(),a.kindName(),a.linkingType(),a.oneLetterCode(),a.hydrogenCount(),a.weight()],[b.ResidueCode.ALA,'ALA','AA',1,'A',7,Math.fround(89.0932)]);
 assert.deepEqual([a.found(),a.isWater(),a.isDna(),a.isRna(),a.isNucleicAcid(),a.isAminoAcid()],[true,false,false,false,false,true]);
 assert.equal(a.isAminoAcid(),true);assert.equal(a.isBufferOrWater(),false);assert.equal(a.isStandard(),true);assert.equal(a.isModifiedAminoAcid(),false);assert.equal(a.isPeptideLinking(),true);assert.equal(a.isNaLinking(),false);assert.equal(a.fastaCode(),'A');assert.equal(a.canonicalOneLetterCode(),'A');assert.equal(a.parentStandardCode(),b.ResidueCode.ALA);
 const m=b.findResidueInfo('MSE');assert.equal(m.oneLetterCode(),'m');assert.equal(m.isModifiedAminoAcid(),true);assert.equal(m.isStandard(),false);assert.equal(m.fastaCode(),'X');assert.equal(m.canonicalOneLetterCode(),'M');assert.equal(m.parentStandardCode(),b.ResidueCode.MET);
 const water=b.findResidueInfo('HOH');assert.equal(water.isWater(),true);assert.equal(water.isBufferOrWater(),true);assert.equal(water.canonicalOneLetterCode(),null);assert.equal(water.parentStandardCode(),null);
 for(const identity of [b.ResidueIdentity.new('mse'),new b.ResidueIdentity('mse')]){assert.equal(identity.name(),'mse');assert.equal(identity.code(),b.ResidueCode.MSE);assert.equal(identity.isTabulated(),true);assert.equal(identity.info().name(),'MSE');}const unknown=b.ResidueIdentity.new('  raw unknown  ');assert.equal(unknown.name(),'  raw unknown  ');assert.equal(unknown.isTabulated(),false);assert.equal(unknown.code(),b.ResidueCode.UNKNOWN);
});
test('BIO sequence expansion preserves character/byte distinctions, parenthesized tokens and typed errors',()=>{
 assert.equal(b.expandOneLetter('a',b.ResidueInfoKind.AA),'ALA');assert.equal(b.expandOneLetter('A',b.ResidueInfoKind.DNA),'DA');assert.equal(b.expandOneLetter('A',b.ResidueInfoKind.RNA),'A');assert.equal(b.expandOneLetter('T',b.ResidueInfoKind.RNA),null);assert.equal(b.expandOneLetter('J',b.ResidueInfoKind.AA),null);assert.equal(b.expandOneLetter('😀',b.ResidueInfoKind.AA),null);
 for(const code of ['','AB'])assert.throws(()=>b.expandOneLetter(code,b.ResidueInfoKind.AA),RangeError);
 assert.deepEqual(b.expandOneLetterSequence(' A\tG(MSE)\n()',b.ResidueInfoKind.AA),['ALA','GLY','MSE','']);assert.deepEqual(b.expandOneLetterSequence('',b.ResidueInfoKind.AA),[]);
 assert.throws(()=>b.expandOneLetterSequence('(MSE',b.ResidueInfoKind.AA),e=>e.kind==='UnmatchedParenthesis'&&e.detail instanceof b.ResidueSequenceError);
 assert.throws(()=>b.expandOneLetterSequence('J',b.ResidueInfoKind.AA),e=>e.kind==='UnexpectedLetter'&&e.sequenceKind==='peptide'&&e.letter==='J'&&e.sourceCode===74);
 assert.throws(()=>b.expandOneLetterSequence('é',b.ResidueInfoKind.AA),e=>e.kind==='UnexpectedLetter'&&e.letter==='Ã'&&e.sourceCode===-61);
});
