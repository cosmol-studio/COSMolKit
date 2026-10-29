import { useEffect, useState } from "react";
import init, { Molecule } from "@cosmol-studio/cosmolkit";

export function MoleculePanel({ smiles = "CCO" }) {
  const [ready, setReady] = useState(false);
  const [molecule, setMolecule] = useState(null);
  const [error, setError] = useState("");

  useEffect(() => {
    let active = true;
    init()
      .then(() => active && setReady(true))
      .catch((value) => active && setError(String(value)));
    return () => {
      active = false;
    };
  }, []);

  useEffect(() => {
    if (!ready) return undefined;
    try {
      const next = Molecule.fromSmiles(smiles);
      setMolecule(next);
      setError("");
    } catch (value) {
      setError(String(value));
    }
    return undefined;
  }, [ready, smiles]);

  useEffect(() => () => molecule?.free(), [molecule]);

  if (error) return <pre>{error}</pre>;
  if (!molecule) return <p>Loading COSMolKit...</p>;
  return <output>{molecule.toSmiles()} ({molecule.numAtoms()} atoms)</output>;
}
