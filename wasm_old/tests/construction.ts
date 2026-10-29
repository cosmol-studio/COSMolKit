import {Molecule,MoleculeBuilder,AtomSpec,BondSpec,Element,BondOrder,Coordinate2DInputParams,CoordinateZPolicy,Coordinate3DInputParams,Replace3DCoordinatesParams,AtomPositionParams} from 'cosmolkit-generated';
const element=Element.fromSymbol('C')!;
const builder=new MoleculeBuilder(),atom=new AtomSpec(element);
builder.addAtom(atom.withFormalCharge(1).withExplicitHydrogens(3).withNoImplicit(true));
builder.addAtom(atom);builder.addBond(new BondSpec(0,1,BondOrder.Single));
const molecule:Molecule=builder.build();
Molecule.fromParts(molecule.topology(),builder.coordinates(),molecule.properties());
const coordinates:Float64Array=molecule.coordinates3d();
molecule.with2dCoordinateBlock([[0,0],[1,1]],{zPolicy:'ignore'});
molecule.set2dCoordinatesWithParams_([[0,0],[1,1]],new Coordinate2DInputParams(CoordinateZPolicy.Ignore));
molecule.withAdded3dConformer([[0,0,0],[1,1,1]],new Coordinate3DInputParams());
molecule.set3dCoordinates_([[0,0,0],[1,1,1]],new Replace3DCoordinatesParams(0));
molecule.withAtomPosition(0,[1,2,3],new AtomPositionParams());
// @ts-expect-error formal charge is checked numeric input
atom.withFormalCharge('1');
// @ts-expect-error coordinate row contains a number
molecule.with2dCoordinateBlock([['1',2]]);
void coordinates;

import {molFromSmiles,molsFromSmilesList,MoleculeBatch} from 'cosmolkit-generated';
const convenient:Molecule=molFromSmiles('CCO',{sanitize:false});
const convenientBatch:MoleculeBatch=molsFromSmilesList(['C','CC']);
// @ts-expect-error the alias preserves the original text input type
molFromSmiles(42);
// @ts-expect-error native binary constructors are absent from WASM
import {molFromBinary} from 'cosmolkit-generated';
void [convenient,convenientBatch];
