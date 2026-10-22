import assert from 'node:assert/strict';
import test from 'node:test';
import {pathToFileURL} from 'node:url';
import {readFileSync} from 'node:fs';
const ck = await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
ck.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});
const pdb = 'HETATM    1  C1  LIG A   7       1.000   2.000   3.000  0.50 12.00           C  \nHETATM    2  H1  LIG A   7       1.800   2.000   3.000  1.00  4.00           H  \nCONECT    1    2\nEND\n';
const cif = 'data_lig\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\nHETATM 1 C C1 . LIG A A 7 1.000 2.000 3.000 0.50 12.00\n';
function positions(text) {
    return text.split('\n').filter(l => l.startsWith('ATOM  ') || l.startsWith('HETATM')).map(l => [Number(l.slice(30,38)),Number(l.slice(38,46)),Number(l.slice(46,54))]);
}
test('PDB writer call forms, dimension selection and structured failures', () => {
    const source = ck.molFromSmiles('C');
    assert.deepEqual(positions(source.toPdb()), [[0,0,0]]);
    const builder = source.toBuilder();
    builder.set2dCoordinates([[1.25,-2.5]]);
    builder.add3dConformer([[4,5,6]]);
    const molecule = builder.build();
    assert.throws(() => molecule.toPdb(), e => e.name === 'MolecularIoError' && e.cause.kind === 'AmbiguousCoordinates' && e.cause.detail instanceof ck.PdbWriteError && e.cause.twoD === 1 && e.cause.threeD === 1);
    assert.deepEqual(positions(molecule.toPdb({coordinateSelection: ck.MolCoordinateSelection.TwoD(0)})), [[1.25,-2.5,0]]);
    const params = new ck.PdbWriteParams(ck.MolCoordinateSelection.ThreeD(0));
    assert.deepEqual(positions(molecule.toPdbWithParams(params)), [[4,5,6]]);
    params.coordinateSelection = null;
    params.allConformers = true;
    assert.deepEqual(positions(molecule.toPdb(params)), [[1.25,-2.5,0],[4,5,6]]);
    params.flavor = ck.PdbWriteFlags.MASTER;
    assert.ok(molecule.toPdb(params).includes('MASTER'));
    assert.throws(() => { params.allConformers = 'true'; }, TypeError);
    assert.equal(params.allConformers, true);
    assert.throws(() => molecule.toPdb({coordinateSelection:ck.MolCoordinateSelection.ThreeD(99)}), e => e.cause.kind === 'MissingCoordinate' && e.cause.coordinateId === 99);
    assert.throws(() => molecule.toPdb({allConformers:true,coordinateSelection:ck.MolCoordinateSelection.TwoD(0)}), e => e.cause.kind === 'ConflictingCoordinateSelection');
    assert.deepEqual(Array.from(molecule.coordinates2d()), [1.25,-2.5]);
    const benzene = ck.molFromSmiles('c1ccccc1');
    assert.ok(benzene.toPdb().includes('CONECT    1    2    6    6\n'));
    assert.equal(benzene.toSmiles(), 'c1ccccc1');
});
test('BIO-enabled molecular readers support params/options and generated conveniences', {skip: !ck.BioMoleculeParams}, () => {
    const params = new ck.BioMoleculeParams(true,false,ck.PdbReadFlags.NONE,false);
    for (const read of [ck.Molecule.fromPdb, ck.molFromPdb]) {
        assert.equal(read(pdb).numAtoms(), 1);
        assert.equal(read(pdb, params).numAtoms(), 2);
        assert.equal(read(pdb, {removeHs:false,proximityBonding:false}).numBonds(), 1);
        assert.deepEqual(positions(read(pdb).toPdb()), [[1,2,3]]);
    }
    assert.equal(ck.Molecule.fromPdbWithParams(pdb, params).numAtoms(), 2);
    for (const read of [ck.Molecule.fromMmcif, ck.molFromMmcif]) {
        assert.equal(read(cif, params).numAtoms(), 1);
        assert.deepEqual(positions(read(cif).toPdb()), [[1,2,3]]);
        assert.equal(read(cif, {proximityBonding:false}).numBonds(), 0);
    }
    assert.equal(ck.Molecule.fromMmcifWithParams(cif, params).numAtoms(), 1);
    assert.throws(() => ck.molFromMmcif('not a CIF block'), e => e.kind === 'MmcifRead' && e.cause.detail instanceof ck.BioMmcifReadError);
});
