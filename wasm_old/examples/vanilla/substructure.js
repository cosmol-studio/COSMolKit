import init, { Molecule, QueryGraph } from "@cosmol-studio/cosmolkit";

await init();

const molecule = Molecule.fromSmiles("CCN.CCO");
try {
  const query = QueryGraph.fromSmarts("[#8]");
  try {
    console.log("contains oxygen:", molecule.hasSubstructMatch(query));
  } finally {
    query.free();
  }
  const fragments = molecule.fragments();
  try {
    console.log("fragment SMILES:", fragments.map((fragment) => fragment.toSmiles()));
  } finally {
    for (const fragment of fragments) fragment.free();
  }
  const largest = molecule.largestFragment();
  try {
    console.log("largest fragment:", largest.toSmiles());
  } finally {
    largest.free();
  }
} finally {
  molecule.free();
}
