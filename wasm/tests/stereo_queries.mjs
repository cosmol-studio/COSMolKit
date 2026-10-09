import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import test from 'node:test';
import {pathToFileURL} from 'node:url';
const b = await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
b.initSync({module: readFileSync(process.env.COSMOLKIT_WASM_BINARY)});

test('Modern chiral-center tuples, default filtering and isolated repeated results', () => {
  // Pinned RDKit 2026.03.1 modern FindMolChiralCenters, includeCIP=True.
  for (const [smiles, expected] of [
    ['C', []], ['C[C@H](O)Cl', [[1, 'R']]], ['C[C@@H](O)Cl', [[1, 'S']]],
    ['CC(O)Cl', [[1, '?']]],
    ['C1C[C@H](C)[C@H](C)[C@H](C)C1', [[2, 'S'], [4, 's'], [6, 'R']]],
    ['F[Pt@SP1](Cl)(Br)I', []],
  ]) {
    const source = b.Molecule.fromSmiles(smiles), before = source.toSmiles();
    try {
      const rows = source.findChiralCenters(true);
      assert.deepEqual(rows, expected);
      assert.deepEqual(source.findChiralCenters(), expected.filter(row => row[1] !== '?'));
      assert.deepEqual(source.findChiralCenters(false), source.findChiralCenters());
      rows.push([999, 'modified']);
      assert.deepEqual(source.findChiralCenters(true), expected);
      assert.equal(source.toSmiles(), before);
    } finally { source.free(); }
  }
});

test('Chiral-center host boolean validation and typed source failures', () => {
  const source = b.Molecule.fromSmiles('CC(O)Cl');
  try { for (const value of [null, 1, 'true', {}]) assert.throws(() => source.findChiralCenters(value), TypeError); }
  finally { source.free(); }
  const params = new b.SmilesParseParams(false, true, true, true, false, true);
  const invalid = b.Molecule.fromSmilesWithParams('C[C@H]1CCCC[C@H]1C |atomProp:1._ringStereochemCand.malformed|', params);
  try {
    assert.throws(() => invalid.findChiralCenters(), e => e instanceof Error && e.name === 'StereoReadError' && e.domain === 'stereo' && e.kind === 'PotentialStereo' && e.detail instanceof b.StereoReadError && e.cause instanceof Error);
  } finally { invalid.free(); params.free(); }
});


test('Chiral-center query preserves pre-existing CIP property text', () => {
  const params = new b.SmilesParseParams(false, true, true, true, false, true);
  const source = b.Molecule.fromSmilesWithParams('CC(O)Cl |atomProp:1._CIPCode.foo|', params);
  try {
    assert.deepEqual(source.findChiralCenters(true), [[1, 'foo']]);
    assert.deepEqual(source.findChiralCenters(), []);
  } finally { source.free(); params.free(); }
});
