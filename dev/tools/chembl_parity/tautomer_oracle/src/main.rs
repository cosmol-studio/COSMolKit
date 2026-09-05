use std::{
    env,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

use cosmolkit_core::{Molecule, TautomerEnumerator, TautomerOptions, score_tautomer};
use serde_json::{Value, json};

fn error_value(error: &dyn std::error::Error) -> Value {
    json!({"type": std::any::type_name_of_val(error), "message": error.to_string()})
}

fn molecule_state(molecule: &Molecule) -> Result<Value, Box<dyn std::error::Error>> {
    let atoms = molecule
        .atoms()
        .iter()
        .map(|atom| {
            json!({
                "atomic_number": atom.atomic_number(),
                "formal_charge": atom.formal_charge(),
                "explicit_hydrogens": atom.explicit_hydrogens(),
                "no_implicit": atom.no_implicit(),
                "isotope": atom.isotope().unwrap_or(0),
                "radical_electrons": atom.radical_electrons(),
                "aromatic": atom.is_aromatic(),
                "chiral_tag": atom.chiral_tag().rdkit_name(),
                "hybridization": format!("{:?}", atom.hybridization()).to_ascii_uppercase(),
                "cip_code": atom.prop("_CIPCode"),
            })
        })
        .collect::<Vec<_>>();
    let bonds = molecule
        .bonds()
        .iter()
        .map(|bond| {
            let stereo_atoms = bond
                .stereo_atoms()
                .map(|atoms| vec![atoms[0].index(), atoms[1].index()])
                .unwrap_or_default();
            json!({
                "begin": bond.begin().index(),
                "end": bond.end().index(),
                "bond_type": bond.order().rdkit_name(),
                "aromatic": bond.is_aromatic(),
                "conjugated": bond.is_conjugated(),
                "direction": bond.direction().rdkit_name(),
                "stereo": bond.stereo().rdkit_name(),
                "stereo_atoms": stereo_atoms,
            })
        })
        .collect::<Vec<_>>();
    Ok(json!({
        "isomeric_smiles": molecule.to_smiles(true)?,
        "atoms": atoms,
        "bonds": bonds,
    }))
}

fn options(branch: &Value) -> TautomerOptions {
    TautomerOptions::default()
        .with_max_tautomers(branch["max_tautomers"].as_u64().unwrap() as u32)
        .with_max_transforms(branch["max_transforms"].as_u64().unwrap() as u32)
        .with_remove_sp3_stereo(branch["remove_sp3_stereo"].as_bool().unwrap())
        .with_remove_bond_stereo(branch["remove_bond_stereo"].as_bool().unwrap())
        .with_remove_isotopic_hydrogens(branch["remove_isotopic_hydrogens"].as_bool().unwrap())
        .with_reassign_stereo(branch["reassign_stereo"].as_bool().unwrap())
}

fn enumerator(branch: &Value) -> Result<TautomerEnumerator<'static>, Box<dyn std::error::Error>> {
    let configured = options(branch);
    if branch["catalog"] == "v1" {
        let mut value = TautomerEnumerator::v1()?;
        value.set_max_tautomers(configured.max_tautomers());
        value.set_max_transforms(configured.max_transforms());
        value.set_remove_sp3_stereo(configured.remove_sp3_stereo());
        value.set_remove_bond_stereo(configured.remove_bond_stereo());
        value.set_remove_isotopic_hydrogens(configured.remove_isotopic_hydrogens());
        value.set_reassign_stereo(configured.reassign_stereo());
        Ok(value)
    } else {
        Ok(TautomerEnumerator::from_options(configured))
    }
}

fn enumerate_branch(molecule: &Molecule, branch: &Value) -> Result<Value, Box<dyn std::error::Error>> {
    let enumerator = enumerator(branch)?;
    let result = enumerator.enumerate(molecule)?;
    let canonical = enumerator.pick_canonical(&result)?;
    let molecule_states = result.iter().map(molecule_state).collect::<Result<Vec<_>, _>>()?;
    let scores = result
        .iter()
        .map(|molecule| {
            let score = score_tautomer(molecule)?;
            Ok::<_, Box<dyn std::error::Error>>(json!({
                "ring": score.ring(),
                "substructure": score.substructure(),
                "hetero_hydrogen": score.hetero_hydrogen(),
                "total": score.total(),
            }))
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(json!({
        "ok": true,
        "error": null,
        "status": format!("{:?}", result.status()),
        "ordered_smiles": result.canonical_smiles(),
        "modified_atoms": result.modified_atoms().iter().map(|id| id.index()).collect::<Vec<_>>(),
        "modified_bonds": result.modified_bonds().iter().map(|id| id.index()).collect::<Vec<_>>(),
        "molecule_states": molecule_states,
        "scores": scores,
        "canonical_smiles": canonical.to_smiles(true)?,
        "canonical_state": molecule_state(&canonical)?,
    }))
}

fn failed_branch(error: &dyn std::error::Error) -> Value {
    json!({
        "ok": false,
        "error": error_value(error),
        "status": null,
        "ordered_smiles": [],
        "modified_atoms": [],
        "modified_bonds": [],
        "molecule_states": [],
        "scores": [],
        "canonical_smiles": null,
        "canonical_state": null,
    })
}

fn argument(name: &str) -> Result<PathBuf, String> {
    let mut arguments = env::args_os();
    while let Some(value) = arguments.next() {
        if value == name {
            return arguments
                .next()
                .map(PathBuf::from)
                .ok_or_else(|| format!("{name} requires a value"));
        }
    }
    Err(format!("missing {name}"))
}

fn optional_argument(name: &str) -> Result<Option<String>, String> {
    let mut arguments = env::args_os();
    while let Some(value) = arguments.next() {
        if value == name {
            return arguments
                .next()
                .map(|value| Some(value.to_string_lossy().into_owned()))
                .ok_or_else(|| format!("{name} requires a value"));
        }
    }
    Ok(None)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input = argument("--input")?;
    let output = argument("--output")?;
    let profile_path = argument("--profile")?;
    let limit = optional_argument("--limit")?
        .map(|value| value.parse::<usize>())
        .transpose()?;
    let profile: Value = serde_json::from_reader(File::open(profile_path)?)?;
    let selected = profile["chembl_branches"]
        .as_array()
        .ok_or("profile chembl_branches must be an array")?;
    let all = profile["branches"]
        .as_array()
        .ok_or("profile branches must be an array")?;
    let branches = selected
        .iter()
        .map(|name| {
            let name = name.as_str().ok_or("branch name must be a string")?;
            let branch = all
                .iter()
                .find(|branch| branch["name"] == name)
                .ok_or_else(|| format!("unknown branch {name}"))?;
            Ok((name.to_owned(), branch.clone()))
        })
        .collect::<Result<Vec<_>, Box<dyn std::error::Error>>>()?;
    let reader = BufReader::new(File::open(Path::new(&input))?);
    let output: Box<dyn Write> = if output == Path::new("-") {
        Box::new(std::io::stdout())
    } else {
        Box::new(File::create(Path::new(&output))?)
    };
    let mut writer = BufWriter::new(output);
    let mut processed = 0usize;
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        if limit.is_some_and(|limit| processed >= limit) {
            break;
        }
        let record: Value = serde_json::from_str(&line)?;
        let row = record["row"].clone();
        let molecule = record["smiles"]
            .as_str()
            .ok_or("record has no string smiles")
            .and_then(|smiles| Molecule::from_smiles(smiles).map_err(|_| "SMILES parse failed"));
        let parse = match &molecule {
            Ok(_) => json!({"ok": true, "error": null}),
            Err(error) => json!({
                "ok": false,
                "error": {"type": "SmilesParseError", "message": error},
            }),
        };
        serde_json::to_writer(
            &mut writer,
            &json!({"kind": "parse", "row": row.clone(), "parse": parse}),
        )?;
        writer.write_all(b"\n")?;
        if let Ok(molecule) = molecule {
            for (name, branch) in &branches {
                let result = enumerate_branch(&molecule, branch).unwrap_or_else(|error| failed_branch(error.as_ref()));
                serde_json::to_writer(
                    &mut writer,
                    &json!({"kind": "branch", "row": row.clone(), "name": name, "result": result}),
                )?;
                writer.write_all(b"\n")?;
            }
        }
        processed += 1;
    }
    writer.flush()?;
    Ok(())
}
