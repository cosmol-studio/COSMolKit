import * as ck from 'cosmolkit-generated';
const mol=ck.Molecule.fromSmiles('CC(O)Cl');
const count:bigint=mol.stereoisomerCount();
const source=ck.StereoisomerRandomSource.fromIntegerSeed(1n<<80n);
const bits=ck.StereoisomerRandomSource.fromRandomBits(width=>0n);
const options=new ck.StereoisomerOptions();options.randomSource=source;options.maxIsomers=10;
options.setOnlyUnassigned(false);
const iterator:ck.StereoisomerIterator=mol.enumerateStereoisomers({onlyUnassigned:false});
const next:ck.Molecule|null=iterator.next();
const yielded:number=iterator.yieldedCount;
mol.enumerateStereoisomersWithRandomBits(options,width=>1n);
// @ts-expect-error arbitrary-width random providers return bigint
ck.StereoisomerRandomSource.fromRandomBits(()=>1);
// @ts-expect-error configuration booleans are not strings
options.unique='true';
