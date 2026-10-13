import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import test from 'node:test';
import {pathToFileURL} from 'node:url';
const ck=await import(pathToFileURL(process.env.COSMOLKIT_WASM_MODULE).href);
ck.initSync({module:readFileSync(process.env.COSMOLKIT_WASM_BINARY)});

test('Top-level constructors retain the original callable, configuration and errors',()=>{
 assert.equal(ck.molFromSmiles,ck.Molecule.fromSmiles);
 const direct=ck.molFromSmiles('CN(=O)=O',{sanitize:false});
 const original=ck.Molecule.fromSmiles('CN(=O)=O',{sanitize:false});
 assert.equal(direct.toSmiles({canonical:false}),original.toSmiles({canonical:false}));
 const errors=[ck.molFromSmiles,ck.Molecule.fromSmiles].map(factory=>{
  try{factory('C1CC');assert.fail('invalid SMILES must fail');}
  catch(error){return [error.name,error.message,error.kind];}
 });
 assert.deepEqual(errors[0],errors[1]);
 assert.equal('molFromBinary' in ck,false);
 assert.equal(ck.molsFromSmilesList,ck.MoleculeBatch.fromSmilesList);
 const batch=ck.molsFromSmilesList(['C','CC']);
 assert.deepEqual(batch.toSmilesList(),['C','CC']);
 direct.free();original.free();batch.free();
});

test('Detached builder validates topology and returns independent molecule and block snapshots',()=>{
 const carbon=ck.Element.fromSymbol('C'),oxygen=ck.Element.fromSymbol('O');
 const atom=new ck.AtomSpec(carbon),builder=new ck.MoleculeBuilder();
 const first=builder.addAtom(atom),second=builder.addAtom(new ck.AtomSpec(oxygen));
 assert.equal(first,0);assert.equal(second,1);
 assert.equal(builder.addBond(new ck.BondSpec(first,second,ck.BondOrder.Single)),0);
 assert.equal(builder.degree(0),1);assert.deepEqual(builder.neighborBonds(0),[0]);
 assert.equal(builder.bondBetweenAtoms(1,0),0);
 const mol=builder.build().sanitize();assert.equal(mol.toSmiles(),'CO');
 const properties=builder.withName('methanol').withProperty('tracking','42').withSdfDataField('field','value').properties();
 assert.equal(properties.name(),'methanol');assert.equal(builder.properties().name(),null);
 const restored=ck.Molecule.fromParts(mol.topology(),builder.coordinates(),properties).sanitize();
 assert.equal(restored.toSmiles(),'CO');assert.equal(restored.property('tracking').asString(),'42');assert.equal(restored.property('absent'),null);
 const edit=mol.toBuilder();edit.setAtomFormalCharge(1,-1);
 assert.equal(edit.atoms()[1].formalCharge(),-1);assert.equal(mol.atom(1).formalCharge(),0);
 edit.setBondOrder(0,'DOUBLE');assert.equal(edit.bonds()[0].order(),ck.BondOrder.Double);
 assert.throws(()=>edit.addBond(new ck.BondSpec(0,99,ck.BondOrder.Single)),e=>e.name==='OperationError'&&e.detail instanceof ck.OperationError);
 assert.equal(edit.bonds().length,1);
 assert.throws(()=>edit.setAtomFormalCharge(0,256),RangeError);assert.equal(edit.atoms()[0].formalCharge(),0);
 assert.equal(edit.removeBondBetweenAtoms(0,1),true);assert.equal(edit.removeBondBetweenAtoms(0,1),false);
 assert.equal(edit.bondBetweenAtoms(0,1),null);assert.equal(mol.numBonds(),1);
 const isotope=builder.addAtom(atom.withIsotope(13).withExplicitHydrogens(4).withNoImplicit(true).withAtomMap(7).withFormalCharge(1));
 const row=builder.atoms()[isotope];assert.equal(row.isotope(),13);assert.equal(row.atomMap(),7);assert.equal(row.formalCharge(),1);
 mol.free();assert.equal(restored.topology().atoms.length,2);assert.equal(properties.name(),'methanol');
});

test('Builder coordinates and conformers retain separate dimensions, validation and snapshot lifetimes',()=>{
 const mol=ck.Molecule.fromSmiles('CO'),builder=mol.toBuilder();
 builder.set2dCoordinates([[1,2],[3,4]]);
 assert.equal(builder.add3dConformer([[1,2,3],[4,5,6]]),0);
 assert.equal(builder.add2dConformer([[5,6],[7,8]]),1);
 const block=builder.coordinates(),out=builder.build();
 assert.equal(out.has2dCoordinates(),true);
 assert.deepEqual([...out.coordinates2d()],[1,2,3,4]);assert.deepEqual([...out.coordinates3d(0)],[1,2,3,4,5,6]);
 assert.equal(block.conformers2d.length,2);assert.equal(block.conformers3d.length,1);
 assert.deepEqual(block.sourceConformerOrder,[ck.CoordinateDimension.TwoD,ck.CoordinateDimension.ThreeD,ck.CoordinateDimension.TwoD]);
 const conformer=out.conformers3d()[0];assert.equal(conformer.id(),0);assert.equal(conformer.is3d(),true);
 const rows=conformer.coordinates();rows[0][0]=99;assert.equal(conformer.coordinates()[0][0],1);
 assert.throws(()=>builder.set2dCoordinates([[0,0]]),e=>e.name==='OperationError'&&e.kind==='InvalidCoordinates');
 assert.deepEqual(builder.coordinates().conformers2d[0].coordinates(),[[1,2],[3,4]]);
 builder.free();out.free();assert.deepEqual(conformer.coordinates(),[[1,2,3],[4,5,6]]);
 assert.deepEqual(block.conformers2d[1].coordinates(),[[5,6],[7,8]]);
 assert.equal(mol.has2dCoordinates(),false);assert.equal(mol.numConformers3d(),0);
});

test('Manual coordinates use canonical shape policies, exact IDs, COW and atomic in-place operations',()=>{
 const source=ck.Molecule.fromSmiles('CO'),original=source.toSmiles();
 const xy=source.with2dCoordinateBlock([[1,2,0],[3,4,0]],{zPolicy:'require_zero'});
 assert.deepEqual([...xy.coordinates2d()],[1,2,3,4]);assert.equal(source.has2dCoordinates(),false);
 const strict=new ck.Coordinate2DInputParams(ck.CoordinateZPolicy.Error);
 assert.throws(()=>source.with2dCoordinateBlockWithParams([[1,2,0],[3,4,0]],strict),e=>e.name==='OperationError'&&e.cause.kind==='ThreeColumnsForbidden'&&e.cause.detail instanceof ck.CoordinateInputError);
 assert.throws(()=>ck.CoordinateZPolicy.fromName('wrong'),e=>e.name==='CoordinateInputError');
 const xyz=source.withAdded3dConformer([[1,2,3],[4,5,6]],{is3d:false});
 assert.equal(xyz.conformers3d()[0].is3d(),false);
 assert.equal(xyz.add3dConformer_([[7,8,9],[10,11,12]]),1);
 assert.deepEqual([...xyz.coordinates3d(1)],[7,8,9,10,11,12]);
 xyz.set3dCoordinates_([[2,3,4],[5,6,7]],{conformerId:1});
 assert.deepEqual([...xyz.coordinates3d(0)],[1,2,3,4,5,6]);
 const moved=xyz.withAtomPosition(0,[9,9,9],{conformerId:1});
 assert.deepEqual([...moved.coordinates3d(1)],[9,9,9,5,6,7]);
 assert.deepEqual([...xyz.coordinates3d(1)],[2,3,4,5,6,7]);
 moved.setAtomPosition_(1,[8,8,8],new ck.AtomPositionParams(1));
 assert.deepEqual([...moved.coordinates3d(1)],[9,9,9,8,8,8]);
 const snapshot=[...xyz.coordinates3d(1)];
 assert.throws(()=>xyz.set3dCoordinates_([[NaN,0,0],[0,0,0]],{conformerId:1}),e=>e.cause.kind==='NonFinite'&&e.cause.row===0);
 assert.deepEqual([...xyz.coordinates3d(1)],snapshot);
 assert.throws(()=>xyz.coordinates3d(99),e=>e.name==='Coordinate3DReadError'&&e.conformerId===99&&e.count===2);
 assert.throws(()=>source.coordinates3d(),e=>e.name==='Coordinate3DReadError');
 assert.throws(()=>xyz.setAtomPosition_(99,[0,0,0]),e=>e.cause.kind==='AtomOutOfRange'&&e.cause.detail instanceof ck.TransformError);
 assert.deepEqual([...xyz.coordinates3d(1)],snapshot);
 const single=xyz.withOnly3dConformer([[0,0,0],[1,1,1]]);assert.equal(single.numConformers3d(),1);
 const cleared=xyz.withCleared3dConformers();assert.equal(cleared.numConformers3d(),0);assert.equal(xyz.numConformers3d(),2);
 xyz.clear3dConformers_();assert.equal(xyz.numConformers3d(),0);
 assert.equal(source.toSmiles(),original);assert.equal(source.hasValenceViolation(0),false);
});
