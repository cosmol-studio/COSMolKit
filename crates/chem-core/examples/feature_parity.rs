use std::collections::VecDeque;
use std::env;
use std::path::PathBuf;
use std::process::Command;

use cosmolkit_chem_core::{
    BondOrder, Molecule, ValenceModel, add_hydrogens_in_place, assign_radicals_rdkit_2025,
    assign_valence,
};

#[derive(Debug, Clone)]
struct AtomLine {
    idx: usize,
    atomic_num: u8,
    degree: usize,
    formal_charge: i8,
    num_hs: usize,
    num_radical_electrons: usize,
    is_aromatic: bool,
    is_in_ring: bool,
    explicit_valence: i32,
    implicit_hs: i32,
    total_valence: i32,
}

impl AtomLine {
    fn as_line(&self) -> String {
        format!(
            "A\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.idx,
            self.atomic_num,
            self.degree,
            self.formal_charge,
            self.num_hs,
            self.num_radical_electrons,
            self.is_aromatic as u8,
            self.is_in_ring as u8,
            self.explicit_valence,
            self.implicit_hs,
            self.total_valence
        )
    }
}

#[derive(Debug, Clone)]
struct BondLine {
    begin_atom: usize,
    end_atom: usize,
    bond_type: String,
}

impl BondLine {
    fn normalize(mut self) -> Self {
        if self.bond_type != "DATIVE" && self.begin_atom > self.end_atom {
            std::mem::swap(&mut self.begin_atom, &mut self.end_atom);
        }
        self
    }

    fn as_line(&self) -> String {
        format!(
            "B\t{}\t{}\t{}",
            self.begin_atom, self.end_atom, self.bond_type
        )
    }
}

fn bond_type_name(order: BondOrder) -> String {
    match order {
        BondOrder::Single => "SINGLE".to_string(),
        BondOrder::Double => "DOUBLE".to_string(),
        BondOrder::Triple => "TRIPLE".to_string(),
        BondOrder::Quadruple => "QUADRUPLE".to_string(),
        BondOrder::Aromatic => "AROMATIC".to_string(),
        BondOrder::Dative => "DATIVE".to_string(),
        BondOrder::Null => "UNSPECIFIED".to_string(),
    }
}

fn compute_ring_flags(mol: &Molecule) -> Vec<bool> {
    let n = mol.atoms.len();
    let mut adj = vec![Vec::<usize>::new(); n];
    for b in &mol.bonds {
        adj[b.begin_atom].push(b.end_atom);
        adj[b.end_atom].push(b.begin_atom);
    }
    let mut in_ring = vec![false; n];

    for a in 0..n {
        if adj[a].len() < 2 {
            continue;
        }
        let neigh = &adj[a];
        let mut found = false;
        for i in 0..neigh.len() {
            for j in (i + 1)..neigh.len() {
                let src = neigh[i];
                let dst = neigh[j];
                let mut q = VecDeque::new();
                let mut seen = vec![false; n];
                seen[a] = true;
                seen[src] = true;
                q.push_back(src);
                while let Some(v) = q.pop_front() {
                    if v == dst {
                        found = true;
                        break;
                    }
                    for &nb in &adj[v] {
                        if !seen[nb] {
                            seen[nb] = true;
                            q.push_back(nb);
                        }
                    }
                }
                if found {
                    break;
                }
            }
            if found {
                break;
            }
        }
        in_ring[a] = found;
    }
    in_ring
}

fn extract_lines(mol: &Molecule) -> Vec<String> {
    let ring_flags = compute_ring_flags(mol);
    let assignment = assign_valence(mol, ValenceModel::RdkitLike)
        .unwrap_or_else(|e| panic!("assign_valence failed: {:?}", e));
    let radicals = assign_radicals_rdkit_2025(mol, &assignment.explicit_valence)
        .unwrap_or_else(|e| panic!("assign_radicals failed: {:?}", e));

    let mut degree = vec![0usize; mol.atoms.len()];
    for b in &mol.bonds {
        degree[b.begin_atom] += 1;
        degree[b.end_atom] += 1;
    }

    let mut out = Vec::new();
    for (i, a) in mol.atoms.iter().enumerate() {
        let implicit_hs = assignment.implicit_hydrogens[i] as i32;
        let explicit_valence = assignment.explicit_valence[i] as i32;
        let num_hs = (a.explicit_hydrogens as i32 + implicit_hs).max(0) as usize;
        let atom = AtomLine {
            idx: i,
            atomic_num: a.atomic_num,
            degree: degree[i] + num_hs,
            formal_charge: a.formal_charge,
            num_hs,
            num_radical_electrons: radicals[i] as usize,
            is_aromatic: a.is_aromatic,
            is_in_ring: ring_flags[i],
            explicit_valence,
            implicit_hs,
            total_valence: explicit_valence + implicit_hs,
        };
        out.push(atom.as_line());
    }

    let mut bonds = Vec::new();
    for b in &mol.bonds {
        bonds.push(
            BondLine {
                begin_atom: b.begin_atom,
                end_atom: b.end_atom,
                bond_type: bond_type_name(b.order),
            }
            .normalize(),
        );
    }
    bonds.sort_by(|l, r| {
        (l.begin_atom, l.end_atom, &l.bond_type).cmp(&(r.begin_atom, r.end_atom, &r.bond_type))
    });
    for b in bonds {
        out.push(b.as_line());
    }

    out
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("crates/")
        .parent()
        .expect("repo root")
        .to_path_buf()
}

fn call_rdkit(smiles: &str) -> String {
    let repo = repo_root();
    let script = repo.join("tests/scripts/rdkit_features_single.py");
    let candidates = [
        env::var("COSMOLKIT_PYTHON").ok(),
        Some(repo.join(".venv/bin/python").display().to_string()),
        Some("python3".to_string()),
    ];

    for py in candidates.iter().flatten() {
        let out = Command::new(py)
            .arg(&script)
            .arg("--smiles")
            .arg(smiles)
            .output();
        match out {
            Ok(o) if o.status.success() => {
                return String::from_utf8(o.stdout)
                    .unwrap_or_else(|_| panic!("rdkit output not utf8 via {}", py));
            }
            Ok(o) => {
                eprintln!(
                    "rdkit call failed via {}: status={} stderr={}",
                    py,
                    o.status,
                    String::from_utf8_lossy(&o.stderr)
                );
            }
            Err(e) => eprintln!("failed to spawn {}: {}", py, e),
        }
    }
    panic!("failed to run rdkit script: {}", script.display());
}

fn parse_sections(output: &str) -> (Vec<String>, Vec<String>) {
    let mut direct = Vec::new();
    let mut with_hs = Vec::new();
    let mut mode = "";
    for raw in output.lines() {
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        if line == "## direct" {
            mode = "direct";
            continue;
        }
        if line == "## with_hs" {
            mode = "with_hs";
            continue;
        }
        if line.starts_with("ERROR\t") {
            panic!("rdkit reported error: {}", line);
        }
        match mode {
            "direct" => direct.push(line.to_string()),
            "with_hs" => with_hs.push(line.to_string()),
            _ => {}
        }
    }
    (direct, with_hs)
}

fn print_diff(label: &str, ours: &[String], rdkit: &[String]) {
    println!("== {} ==", label);
    if ours == rdkit {
        println!("MATCH");
        return;
    }
    println!("DIFF");
    let n = ours.len().max(rdkit.len());
    for i in 0..n {
        let l = ours.get(i).map_or("<none>", String::as_str);
        let r = rdkit.get(i).map_or("<none>", String::as_str);
        if l != r {
            println!("first diff at line {}", i + 1);
            println!("OURS : {}", l);
            println!("RDKIT: {}", r);
            break;
        }
    }
    println!("-- OURS --");
    for line in ours {
        println!("{}", line);
    }
    println!("-- RDKIT --");
    for line in rdkit {
        println!("{}", line);
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!(
            "Usage: cargo run -p cosmolkit-chem-core --example feature_parity -- \"<SMILES>\""
        );
        std::process::exit(2);
    }
    let smiles = &args[1];

    let mol = Molecule::from_smiles(smiles).unwrap_or_else(|e| panic!("parse failed: {}", e));
    let ours_direct = extract_lines(&mol);

    let mut with_h = mol.clone();
    add_hydrogens_in_place(&mut with_h).unwrap_or_else(|e| panic!("AddHs failed: {:?}", e));
    let ours_with_hs = extract_lines(&with_h);

    let rdkit_output = call_rdkit(smiles);
    let (rdkit_direct, rdkit_with_hs) = parse_sections(&rdkit_output);

    println!("{:?}, {:?}", ours_direct, rdkit_direct);

    print_diff("direct", &ours_direct, &rdkit_direct);
    print_diff("with_hs", &ours_with_hs, &rdkit_with_hs);

    if ours_direct != rdkit_direct || ours_with_hs != rdkit_with_hs {
        std::process::exit(1);
    }
}
