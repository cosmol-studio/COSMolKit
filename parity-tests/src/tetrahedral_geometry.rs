//! Historical ordered-ligand geometry contract, using prepared RDKit coordinates.
use crate::{
    Result, tasks,
    value::{self, Value},
};
use cosmolkit::{BatchRecord, LigandRef, Molecule, MoleculeBatch};
use rayon::prelude::*;

fn verify(mol: &Molecule, positions: &[[f64; 3]]) -> Result<Value> {
    if positions.len() != mol.num_atoms() {
        return Err("reference coordinate count differs from molecule".into());
    }
    let mut rows = Vec::new();
    for stereo in mol.tetrahedral_stereo(false)? {
        let center = positions[stereo.center.index()];
        let mut vectors = [[0.0; 3]; 3];
        for (vector, ligand) in vectors.iter_mut().zip(
            stereo
                .ligands
                .ok_or("specified center has no configuration")?,
        ) {
            let LigandRef::Atom(atom) = ligand else {
                return Err("non-atom ligand in first three geometry positions".into());
            };
            for axis in 0..3 {
                vector[axis] = positions[atom.index()][axis] - center[axis];
            }
        }
        let [a, b, c] = vectors;
        let volume = a[0] * (b[1] * c[2] - b[2] * c[1])
            + a[1] * (b[2] * c[0] - b[0] * c[2])
            + a[2] * (b[0] * c[1] - b[1] * c[0]);
        // Preserve 24cfdb1's strict positive-volume assertion, not abs(volume).
        rows.push(Value::List(vec![
            Value::Int(stereo.center.index() as i64),
            Value::Bool(volume > 1.0e-8),
        ]));
    }
    Ok(Value::List(rows))
}

fn observe(input: &str, expected: &[u8]) -> Result<Value> {
    if expected.first() != Some(&0) {
        return Err("RDKit geometry unavailable; comparison unfinished".into());
    }
    let Value::List(witness) = value::decode(&expected[1..])? else {
        return Err("invalid geometry reference".into());
    };
    let Some(Value::List(coordinates)) = witness.first() else {
        return Err("missing reference coordinates".into());
    };
    let positions = coordinates
        .iter()
        .map(|row| {
            let Value::List(values) = row else {
                return Err("invalid coordinate row".into());
            };
            let [Value::Float(x), Value::Float(y), Value::Float(z)] = values.as_slice() else {
                return Err("invalid coordinate shape/type".into());
            };
            Ok([f64::from_bits(*x), f64::from_bits(*y), f64::from_bits(*z)])
        })
        .collect::<Result<Vec<_>>>()?;
    let row: serde_json::Value = serde_json::from_str(input)?;
    let smiles = row["smiles"].as_str().ok_or("missing smiles")?;
    let mol = Molecule::from_smiles(smiles)?;
    let scalar = verify(&mol, &positions)?;
    let batch = MoleculeBatch::from_smiles_list(&[smiles.to_owned(), smiles.to_owned()])?;
    let mut batch_value = None;
    for record in batch.iter() {
        let BatchRecord::Molecule(mol) = record else {
            return Err("batch parse failed".into());
        };
        let observed = verify(mol, &positions)?;
        if observed != scalar {
            return Err("batch differs from scalar geometry".into());
        }
        batch_value = Some(observed);
    }
    // Only the independent coordinates are reused. Centers and orientation
    // flags are calculated from CK's actual public results, never copied.
    Ok(Value::List(vec![
        witness[0].clone(),
        scalar,
        batch_value.ok_or("empty batch")?,
    ]))
}

pub fn evaluate(
    inputs: &[String],
    references: &[Vec<u8>],
    pool: &rayon::ThreadPool,
) -> Result<Vec<Vec<Vec<u8>>>> {
    pool.install(|| {
        inputs
            .par_iter()
            .zip(references)
            .map(|(input, reference)| {
                let expected = tasks::unpack_row(reference, 1)?;
                Ok(vec![match observe(input, expected[0]) {
                    Ok(value) => tasks::success(value),
                    Err(error) => tasks::failed_error(error.as_ref()),
                }])
            })
            .collect()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn geometry_predicate_rejects_opposite_handedness() {
        let positions = [[1., 0., 0.], [0., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
        for (smiles, positive) in [("F[C@H](Cl)Br", true), ("F[C@@H](Cl)Br", false)] {
            let mol = Molecule::from_smiles(smiles).unwrap();
            assert_eq!(
                verify(&mol, &positions).unwrap(),
                Value::List(vec![Value::List(vec![
                    Value::Int(1),
                    Value::Bool(positive)
                ])])
            );
        }
    }
}
