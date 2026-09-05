use std::convert::TryFrom;
use std::env;
use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

fn main() {
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR"));
    let out_dir = PathBuf::from(env::var("OUT_DIR").expect("OUT_DIR"));
    let generated_path = out_dir.join("forcefield_generated.rs");
    let uff_generated_path = out_dir.join("uff_defaults_generated.rs");
    let mmff_generated_path = out_dir.join("mmff_defaults_generated.rs");
    let crystalff_generated_path = out_dir.join("crystalff_defaults_generated.rs");
    let isotope_generated_path = out_dir.join("rdkit_isotope_masses_generated.rs");

    let uff_params_path = manifest_dir.join("src/chemistry/forcefield/rdkit/ForceField/UFF/Params.cpp");
    let mmff_params_path = manifest_dir.join("src/chemistry/forcefield/rdkit/ForceField/MMFF/Params.cpp");
    let crystalff_dir = manifest_dir.join("src/chemistry/forcefield/rdkit/GraphMol/ForceFieldHelpers/CrystalFF");
    let rdkit_atomic_data_path = manifest_dir.join("src/chemistry/forcefield/rdkit/GraphMol/atomic_data.cpp");

    println!("cargo:rerun-if-changed={}", uff_params_path.display());
    println!("cargo:rerun-if-changed={}", mmff_params_path.display());
    println!("cargo:rerun-if-changed={}", rdkit_atomic_data_path.display());
    for file_name in [
        "torsionPreferences_v1.in",
        "torsionPreferences_v2.in",
        "torsionPreferences_smallrings.in",
        "torsionPreferences_macrocycles.in",
    ] {
        println!("cargo:rerun-if-changed={}", crystalff_dir.join(file_name).display());
    }

    let uff_source = fs::read_to_string(&uff_params_path).expect("read RDKit UFF Params.cpp");
    let uff_param_data =
        extract_cpp_string(&uff_source, "defaultParamData").expect("extract RDKit UFF defaultParamData");
    let mut uff_atomic_params = parse_uff_atomic_params(&uff_param_data);
    uff_atomic_params.sort_by(|left, right| left.label.cmp(&right.label));

    let mmff_source = fs::read_to_string(&mmff_params_path).expect("read RDKit MMFF Params.cpp");
    let mmff_symbols = [
        "defaultMMFFDef",
        "defaultMMFFProp",
        "defaultMMFFPBCI",
        "defaultMMFFChg",
        "defaultMMFFBond",
        "defaultMMFFBndk",
        "defaultMMFFHerschbachLaurie",
        "defaultMMFFCovRadPauEle",
        "defaultMMFFStbn",
        "defaultMMFFDfsb",
        "defaultMMFFOop",
        "defaultMMFFsOop",
        "defaultMMFFTor",
        "defaultMMFFsTor",
        "defaultMMFFVdW",
    ];
    let mmff_string_defaults: Vec<(&str, String)> = mmff_symbols
        .iter()
        .map(|&symbol| {
            (
                symbol,
                extract_cpp_string(&mmff_source, symbol).unwrap_or_else(|| panic!("extract RDKit MMFF {symbol}")),
            )
        })
        .collect();
    let mmff_angle = extract_cpp_string_array(&mmff_source, "defaultMMFFAngleData", "")
        .expect("extract RDKit MMFF defaultMMFFAngleData");
    let mut mmff_def_rows = parse_mmff_def_rows(mmff_default(&mmff_string_defaults, "defaultMMFFDef"));
    mmff_def_rows.sort_by_key(|row| row.atom_type);
    let mut mmff_prop_rows = parse_mmff_prop_rows(mmff_default(&mmff_string_defaults, "defaultMMFFProp"));
    mmff_prop_rows.sort_by_key(|row| row.atom_type);
    let mut mmff_pbci_rows = parse_mmff_pbci_rows(mmff_default(&mmff_string_defaults, "defaultMMFFPBCI"));
    mmff_pbci_rows.sort_by_key(|row| row.atom_type);
    let mut mmff_chg_rows = parse_mmff_chg_rows(mmff_default(&mmff_string_defaults, "defaultMMFFChg"));
    mmff_chg_rows.sort_by_key(|row| (row.bond_type, row.i_atom_type, row.j_atom_type));
    let mut mmff_bond_rows = parse_mmff_bond_rows(mmff_default(&mmff_string_defaults, "defaultMMFFBond"));
    mmff_bond_rows.sort_by_key(|row| (row.bond_type, row.i_atom_type, row.j_atom_type));
    let mut mmff_bndk_rows = parse_mmff_bndk_rows(mmff_default(&mmff_string_defaults, "defaultMMFFBndk"));
    mmff_bndk_rows.sort_by_key(|row| (row.i_atomic_num, row.j_atomic_num));
    let mut mmff_herschbach_laurie_rows =
        parse_mmff_herschbach_laurie_rows(mmff_default(&mmff_string_defaults, "defaultMMFFHerschbachLaurie"));
    mmff_herschbach_laurie_rows.sort_by_key(|row| (row.i_row, row.j_row));
    let mut mmff_cov_rad_pau_ele_rows =
        parse_mmff_cov_rad_pau_ele_rows(mmff_default(&mmff_string_defaults, "defaultMMFFCovRadPauEle"));
    mmff_cov_rad_pau_ele_rows.sort_by_key(|row| row.atomic_num);
    let mut mmff_angle_rows = parse_mmff_angle_rows(&mmff_angle);
    mmff_angle_rows.sort_by_key(|row| (row.angle_type, row.i_atom_type, row.j_atom_type, row.k_atom_type));
    let mut mmff_stbn_rows = parse_mmff_stbn_rows(mmff_default(&mmff_string_defaults, "defaultMMFFStbn"));
    mmff_stbn_rows.sort_by_key(|row| (row.stretch_bend_type, row.i_atom_type, row.j_atom_type, row.k_atom_type));
    let mut mmff_dfsb_rows = parse_mmff_dfsb_rows(mmff_default(&mmff_string_defaults, "defaultMMFFDfsb"));
    mmff_dfsb_rows.sort_by_key(|row| (row.i_atomic_num, row.j_atomic_num, row.k_atomic_num));
    let mut mmff_oop_rows = parse_mmff_oop_rows(mmff_default(&mmff_string_defaults, "defaultMMFFOop"));
    mmff_oop_rows.sort_by_key(|row| (row.i_atom_type, row.j_atom_type, row.k_atom_type, row.l_atom_type));
    let mut mmffs_oop_rows = parse_mmff_oop_rows(mmff_default(&mmff_string_defaults, "defaultMMFFsOop"));
    mmffs_oop_rows.sort_by_key(|row| (row.i_atom_type, row.j_atom_type, row.k_atom_type, row.l_atom_type));
    let mut mmff_tor_rows = parse_mmff_tor_rows(mmff_default(&mmff_string_defaults, "defaultMMFFTor"));
    mmff_tor_rows.sort_by_key(|row| {
        (
            row.tor_type,
            row.i_atom_type,
            row.j_atom_type,
            row.k_atom_type,
            row.l_atom_type,
        )
    });
    let mut mmffs_tor_rows = parse_mmff_tor_rows(mmff_default(&mmff_string_defaults, "defaultMMFFsTor"));
    mmffs_tor_rows.sort_by_key(|row| {
        (
            row.tor_type,
            row.i_atom_type,
            row.j_atom_type,
            row.k_atom_type,
            row.l_atom_type,
        )
    });
    let mut mmff_vdw_data = parse_mmff_vdw_rows(mmff_default(&mmff_string_defaults, "defaultMMFFVdW"));
    mmff_vdw_data.rows.sort_by_key(|row| row.atom_type);

    let crystalff_defaults = [
        (
            "TORSION_PREFERENCES_V1",
            "torsionPreferences_v1.in",
            "torsionPreferencesV1",
        ),
        (
            "TORSION_PREFERENCES_V2",
            "torsionPreferences_v2.in",
            "torsionPreferencesV2",
        ),
        (
            "TORSION_PREFERENCES_SMALL_RINGS",
            "torsionPreferences_smallrings.in",
            "torsionPreferencesSmallRings",
        ),
        (
            "TORSION_PREFERENCES_MACROCYCLES",
            "torsionPreferences_macrocycles.in",
            "torsionPreferencesMacrocycles",
        ),
    ];
    let crystalff_strings: Vec<(&str, String)> = crystalff_defaults
        .iter()
        .map(|&(rust_name, file_name, cpp_symbol)| {
            let source = fs::read_to_string(crystalff_dir.join(file_name)).expect("read CrystalFF table");
            (
                rust_name,
                extract_cpp_string(&source, cpp_symbol)
                    .unwrap_or_else(|| panic!("extract RDKit CrystalFF {cpp_symbol}")),
            )
        })
        .collect();
    let rdkit_atomic_data_source = fs::read_to_string(&rdkit_atomic_data_path).expect("read RDKit atomic_data.cpp");
    let periodic_table_atom_data = extract_cpp_string(&rdkit_atomic_data_source, "periodicTableAtomData")
        .expect("extract RDKit periodicTableAtomData");
    let isotope_atom_data = extract_cpp_string_array(&rdkit_atomic_data_source, "isotopesAtomData", "\n")
        .expect("extract RDKit isotopesAtomData");
    let mut isotope_masses = parse_rdkit_isotope_masses(&isotope_atom_data);
    isotope_masses.sort_by_key(|row| (row.atomic_number, row.isotope));
    assert_eq!(
        isotope_masses
            .binary_search_by_key(&(20, 47), |row| (row.atomic_number, row.isotope))
            .ok()
            .map(|index| isotope_masses[index].mass),
        Some(46.954546),
        "RDKit isotope array element boundaries must preserve the Ca-47 source row"
    );
    let most_common_isotopes = parse_rdkit_most_common_isotopes(&periodic_table_atom_data);
    let most_common_isotope_masses = parse_rdkit_most_common_isotope_masses(&periodic_table_atom_data);

    let mut generated = String::new();
    writeln!(
        generated,
        "// @generated by crates/cosmolkit-core/build.rs; do not edit by hand."
    )
    .unwrap();
    writeln!(generated).unwrap();
    write_string_const(&mut generated, "UFF_DEFAULT_PARAM_DATA", &uff_param_data);
    writeln!(
        generated,
        "pub(crate) static UFF_DEFAULT_ATOMIC_PARAMS: &[(&str, AtomicParams)] = &["
    )
    .unwrap();
    for row in &uff_atomic_params {
        // RDKit✔️✔️: paramObj.theta0 = boost::lexical_cast<double>(*token);
        // RDKit✔️✔️: paramObj.theta0 = paramObj.theta0 * M_PI / 180.;
        writeln!(
            generated,
            "    ({:?}, AtomicParams {{ r1: {:?}, theta0: {:?} * std::f64::consts::PI / 180.0, x1: {:?}, d1: {:?}, zeta: {:?}, z1: {:?}, v1: {:?}, u1: {:?}, gmp_xi: {:?}, gmp_hardness: {:?}, gmp_radius: {:?} }}),",
            row.label,
            row.r1,
            row.theta0_degrees,
            row.x1,
            row.d1,
            row.zeta,
            row.z1,
            row.v1,
            row.u1,
            row.gmp_xi,
            row.gmp_hardness,
            row.gmp_radius,
        )
        .unwrap();
    }
    writeln!(generated, "];").unwrap();
    writeln!(generated).unwrap();

    let mut uff_generated = generated_header();
    write_string_const(&mut uff_generated, "UFF_DEFAULT_PARAM_DATA", &uff_param_data);
    writeln!(
        uff_generated,
        "pub(crate) static UFF_DEFAULT_ATOMIC_PARAMS: &[(&str, AtomicParams)] = &["
    )
    .unwrap();
    for row in &uff_atomic_params {
        // RDKit✔️✔️: paramObj.theta0 = boost::lexical_cast<double>(*token);
        // RDKit✔️✔️: paramObj.theta0 = paramObj.theta0 * M_PI / 180.;
        writeln!(
            uff_generated,
            "    ({:?}, AtomicParams {{ r1: {:?}, theta0: {:?} * std::f64::consts::PI / 180.0, x1: {:?}, d1: {:?}, zeta: {:?}, z1: {:?}, v1: {:?}, u1: {:?}, gmp_xi: {:?}, gmp_hardness: {:?}, gmp_radius: {:?} }}),",
            row.label,
            row.r1,
            row.theta0_degrees,
            row.x1,
            row.d1,
            row.zeta,
            row.z1,
            row.v1,
            row.u1,
            row.gmp_xi,
            row.gmp_hardness,
            row.gmp_radius,
        )
        .unwrap();
    }
    writeln!(uff_generated, "];").unwrap();
    fs::write(uff_generated_path, uff_generated).expect("write UFF generated Rust");

    let mut mmff_generated = generated_header();
    for (symbol, value) in &mmff_string_defaults {
        let rust_name = mmff_symbol_to_rust(symbol);
        write_string_const(&mut mmff_generated, rust_name, value);
    }
    write_string_const(&mut mmff_generated, "DEFAULT_MMFF_ANGLE", &mmff_angle);
    write_mmff_generated_tables(
        &mut mmff_generated,
        &mmff_def_rows,
        &mmff_prop_rows,
        &mmff_pbci_rows,
        &mmff_chg_rows,
        &mmff_bond_rows,
        &mmff_bndk_rows,
        &mmff_herschbach_laurie_rows,
        &mmff_cov_rad_pau_ele_rows,
        &mmff_angle_rows,
        &mmff_stbn_rows,
        &mmff_dfsb_rows,
        &mmff_oop_rows,
        &mmffs_oop_rows,
        &mmff_tor_rows,
        &mmffs_tor_rows,
        &mmff_vdw_data,
    );
    fs::write(mmff_generated_path, mmff_generated).expect("write MMFF generated Rust");

    let mut crystalff_generated = generated_header();
    for (rust_name, value) in &crystalff_strings {
        write_string_const(&mut crystalff_generated, rust_name, value);
    }
    fs::write(crystalff_generated_path, crystalff_generated).expect("write CrystalFF generated Rust");

    let mut isotope_generated = generated_header();
    writeln!(
        isotope_generated,
        "pub(crate) static RDKIT_ISOTOPE_MASSES: &[((u8, u16), f64)] = &["
    )
    .unwrap();
    for row in &isotope_masses {
        writeln!(
            isotope_generated,
            "    (({}, {}), {:?}),",
            row.atomic_number, row.isotope, row.mass
        )
        .unwrap();
    }
    writeln!(isotope_generated, "];").unwrap();
    writeln!(isotope_generated).unwrap();
    writeln!(
        isotope_generated,
        "pub(crate) static RDKIT_MOST_COMMON_ISOTOPES: &[i16] = &["
    )
    .unwrap();
    for isotope in &most_common_isotopes {
        writeln!(isotope_generated, "    {},", isotope).unwrap();
    }
    writeln!(isotope_generated, "];").unwrap();
    writeln!(isotope_generated).unwrap();
    writeln!(
        isotope_generated,
        "pub(crate) static RDKIT_MOST_COMMON_ISOTOPE_MASSES: &[f64] = &["
    )
    .unwrap();
    for mass in &most_common_isotope_masses {
        writeln!(isotope_generated, "    {:?},", mass).unwrap();
    }
    writeln!(isotope_generated, "];").unwrap();
    fs::write(isotope_generated_path, isotope_generated).expect("write RDKit isotope mass generated Rust");

    fs::write(generated_path, generated).expect("write forcefield generated Rust");
}

fn mmff_default<'a>(defaults: &'a [(&str, String)], symbol: &str) -> &'a str {
    defaults
        .iter()
        .find_map(|(candidate, value)| (*candidate == symbol).then_some(value.as_str()))
        .unwrap_or_else(|| panic!("missing extracted MMFF default {symbol}"))
}

#[derive(Debug)]
struct MmffDefRow {
    atom_type: u8,
    eq_level: [u8; 4],
}

#[derive(Debug)]
struct MmffPropRow {
    atom_type: u8,
    prop: [u8; 8],
}

#[derive(Debug)]
struct MmffPbciRow {
    atom_type: u8,
    pbci: f64,
    fcadj: f64,
}

#[derive(Debug)]
struct MmffChgRow {
    bond_type: u8,
    i_atom_type: u8,
    j_atom_type: u8,
    bci: f64,
}

#[derive(Debug)]
struct MmffBondRow {
    bond_type: u8,
    i_atom_type: u8,
    j_atom_type: u8,
    kb: f64,
    r0: f64,
}

#[derive(Debug)]
struct MmffBndkRow {
    i_atomic_num: u8,
    j_atomic_num: u8,
    r0: f64,
    kb: f64,
}

#[derive(Debug)]
struct MmffHerschbachLaurieRow {
    i_row: u8,
    j_row: u8,
    a_ij: f64,
    d_ij: f64,
    dp_ij: f64,
}

#[derive(Debug)]
struct MmffCovRadPauEleRow {
    atomic_num: u8,
    r0: f64,
    chi: f64,
}

#[derive(Debug)]
struct MmffAngleRow {
    angle_type: u8,
    i_atom_type: u8,
    j_atom_type: u8,
    k_atom_type: u8,
    ka: f64,
    theta0: f64,
}

#[derive(Debug)]
struct MmffStbnRow {
    stretch_bend_type: u8,
    i_atom_type: u8,
    j_atom_type: u8,
    k_atom_type: u8,
    kba_ijk: f64,
    kba_kji: f64,
}

#[derive(Debug)]
struct MmffDfsbRow {
    i_atomic_num: u32,
    j_atomic_num: u32,
    k_atomic_num: u32,
    kba_ijk: f64,
    kba_kji: f64,
}

#[derive(Debug)]
struct MmffOopRow {
    i_atom_type: u8,
    j_atom_type: u8,
    k_atom_type: u8,
    l_atom_type: u8,
    koop: f64,
}

#[derive(Debug)]
struct MmffTorRow {
    tor_type: u8,
    i_atom_type: u8,
    j_atom_type: u8,
    k_atom_type: u8,
    l_atom_type: u8,
    v1: f64,
    v2: f64,
    v3: f64,
}

#[derive(Debug)]
struct MmffVdwRow {
    atom_type: u8,
    alpha_i: f64,
    n_i: f64,
    a_i: f64,
    g_i: f64,
    r_star: f64,
    da: u8,
}

#[derive(Debug)]
struct MmffVdwData {
    power: f64,
    b: f64,
    beta: f64,
    darad: f64,
    daeps: f64,
    rows: Vec<MmffVdwRow>,
}

fn parse_mmff_def_rows(mmff_def: &str) -> Vec<MmffDefRow> {
    let mut rows = Vec::new();
    for (line_idx, line) in mmff_def.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        assert!(columns.len() >= 6, "malformed MMFFDef line {}", line_idx + 1);
        rows.push(MmffDefRow {
            atom_type: parse_u8(columns[1], line),
            eq_level: [
                parse_u8(columns[2], line),
                parse_u8(columns[3], line),
                parse_u8(columns[4], line),
                parse_u8(columns[5], line),
            ],
        });
    }
    rows
}

fn parse_mmff_prop_rows(mmff_prop: &str) -> Vec<MmffPropRow> {
    let mut rows = Vec::new();
    for (line_idx, line) in mmff_prop.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        assert!(columns.len() >= 9, "malformed MMFFProp line {}", line_idx + 1);
        rows.push(MmffPropRow {
            atom_type: parse_u8(columns[0], line),
            prop: [
                parse_u8(columns[1], line),
                parse_u8(columns[2], line),
                parse_u8(columns[3], line),
                parse_u8(columns[4], line),
                parse_u8(columns[5], line),
                parse_u8(columns[6], line),
                parse_u8(columns[7], line),
                parse_u8(columns[8], line),
            ],
        });
    }
    rows
}

fn parse_mmff_pbci_rows(mmff_pbci: &str) -> Vec<MmffPbciRow> {
    let mut rows = Vec::new();
    for (line_idx, line) in mmff_pbci.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        assert!(columns.len() >= 4, "malformed MMFFPBCI line {}", line_idx + 1);
        rows.push(MmffPbciRow {
            atom_type: parse_u8(columns[0], line),
            pbci: parse_f64(columns[2], line),
            fcadj: parse_f64(columns[3], line),
        });
    }
    rows
}

fn parse_mmff_chg_rows(mmff_chg: &str) -> Vec<MmffChgRow> {
    let mut rows = Vec::new();
    for (line_idx, line) in mmff_chg.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        assert!(columns.len() >= 4, "malformed MMFFChg line {}", line_idx + 1);
        rows.push(MmffChgRow {
            bond_type: parse_u8(columns[0], line),
            i_atom_type: parse_u8(columns[1], line),
            j_atom_type: parse_u8(columns[2], line),
            bci: parse_f64(columns[3], line),
        });
    }
    rows
}

fn parse_mmff_bond_rows(mmff_bond: &str) -> Vec<MmffBondRow> {
    let mut rows = Vec::new();
    for (line_idx, line) in mmff_bond.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        assert!(columns.len() >= 5, "malformed MMFFBond line {}", line_idx + 1);
        rows.push(MmffBondRow {
            bond_type: parse_u8(columns[0], line),
            i_atom_type: parse_u8(columns[1], line),
            j_atom_type: parse_u8(columns[2], line),
            kb: parse_f64(columns[3], line),
            r0: parse_f64(columns[4], line),
        });
    }
    rows
}

fn parse_mmff_bndk_rows(mmff_bndk: &str) -> Vec<MmffBndkRow> {
    let mut rows = Vec::new();
    for (line_idx, line) in mmff_bndk.lines().enumerate() {
        if line.trim().is_empty() || line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        assert!(columns.len() >= 4, "malformed MMFFBndk line {}", line_idx + 1);
        rows.push(MmffBndkRow {
            i_atomic_num: parse_u8(columns[0], line),
            j_atomic_num: parse_u8(columns[1], line),
            r0: parse_f64(columns[2], line),
            kb: parse_f64(columns[3], line),
        });
    }
    rows
}

fn parse_mmff_herschbach_laurie_rows(mmff_herschbach_laurie: &str) -> Vec<MmffHerschbachLaurieRow> {
    let mut rows = Vec::new();
    for (line_idx, line) in mmff_herschbach_laurie.lines().enumerate() {
        if line.trim().is_empty() || line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        assert!(
            columns.len() >= 5,
            "malformed MMFFHerschbachLaurie line {}",
            line_idx + 1
        );
        rows.push(MmffHerschbachLaurieRow {
            i_row: parse_u8(columns[0], line),
            j_row: parse_u8(columns[1], line),
            a_ij: parse_f64(columns[2], line),
            d_ij: parse_f64(columns[3], line),
            dp_ij: parse_f64(columns[4], line),
        });
    }
    rows
}

fn parse_mmff_cov_rad_pau_ele_rows(mmff_cov_rad_pau_ele: &str) -> Vec<MmffCovRadPauEleRow> {
    let mut rows = Vec::new();
    for (line_idx, line) in mmff_cov_rad_pau_ele.lines().enumerate() {
        if line.trim().is_empty() || line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        assert!(columns.len() >= 3, "malformed MMFFCovRadPauEle line {}", line_idx + 1);
        rows.push(MmffCovRadPauEleRow {
            atomic_num: parse_u8(columns[0], line),
            r0: parse_f64(columns[1], line),
            chi: parse_f64(columns[2], line),
        });
    }
    rows
}

fn parse_mmff_angle_rows(mmff_angle: &str) -> Vec<MmffAngleRow> {
    let mut rows = Vec::new();
    for (line_idx, line) in mmff_angle.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        if line == "EOS" {
            break;
        }
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        assert!(columns.len() >= 6, "malformed MMFFAngle line {}", line_idx + 1);
        rows.push(MmffAngleRow {
            angle_type: parse_u8(columns[0], line),
            i_atom_type: parse_u8(columns[1], line),
            j_atom_type: parse_u8(columns[2], line),
            k_atom_type: parse_u8(columns[3], line),
            ka: parse_f64(columns[4], line),
            theta0: parse_f64(columns[5], line),
        });
    }
    rows
}

fn parse_mmff_stbn_rows(mmff_stbn: &str) -> Vec<MmffStbnRow> {
    let mut rows = Vec::new();
    for (line_idx, line) in mmff_stbn.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        assert!(columns.len() >= 6, "malformed MMFFStbn line {}", line_idx + 1);
        rows.push(MmffStbnRow {
            stretch_bend_type: parse_u8(columns[0], line),
            i_atom_type: parse_u8(columns[1], line),
            j_atom_type: parse_u8(columns[2], line),
            k_atom_type: parse_u8(columns[3], line),
            kba_ijk: parse_f64(columns[4], line),
            kba_kji: parse_f64(columns[5], line),
        });
    }
    rows
}

fn parse_mmff_dfsb_rows(mmff_dfsb: &str) -> Vec<MmffDfsbRow> {
    let mut rows = Vec::new();
    for (line_idx, line) in mmff_dfsb.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        assert!(columns.len() >= 5, "malformed MMFFDfsb line {}", line_idx + 1);
        rows.push(MmffDfsbRow {
            i_atomic_num: parse_u32(columns[0], line),
            j_atomic_num: parse_u32(columns[1], line),
            k_atomic_num: parse_u32(columns[2], line),
            kba_ijk: parse_f64(columns[3], line),
            kba_kji: parse_f64(columns[4], line),
        });
    }
    rows
}

fn parse_mmff_oop_rows(mmff_oop: &str) -> Vec<MmffOopRow> {
    let mut rows = Vec::new();
    for (line_idx, line) in mmff_oop.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        assert!(columns.len() >= 5, "malformed MMFFOop line {}", line_idx + 1);
        rows.push(MmffOopRow {
            i_atom_type: parse_u8(columns[0], line),
            j_atom_type: parse_u8(columns[1], line),
            k_atom_type: parse_u8(columns[2], line),
            l_atom_type: parse_u8(columns[3], line),
            koop: parse_f64(columns[4], line),
        });
    }
    rows
}

fn parse_mmff_tor_rows(mmff_tor: &str) -> Vec<MmffTorRow> {
    let mut rows = Vec::new();
    for (line_idx, line) in mmff_tor.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        assert!(columns.len() >= 8, "malformed MMFFTor line {}", line_idx + 1);
        rows.push(MmffTorRow {
            tor_type: parse_u8(columns[0], line),
            i_atom_type: parse_u8(columns[1], line),
            j_atom_type: parse_u8(columns[2], line),
            k_atom_type: parse_u8(columns[3], line),
            l_atom_type: parse_u8(columns[4], line),
            v1: parse_f64(columns[5], line),
            v2: parse_f64(columns[6], line),
            v3: parse_f64(columns[7], line),
        });
    }
    rows
}

fn parse_mmff_vdw_rows(mmff_vdw: &str) -> MmffVdwData {
    let mut power = None;
    let mut b = None;
    let mut beta = None;
    let mut darad = None;
    let mut daeps = None;
    let mut rows = Vec::new();
    for (line_idx, line) in mmff_vdw.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with('*') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() == 5 && columns[0].parse::<f64>().is_ok() {
            power = Some(parse_f64(columns[0], line));
            b = Some(parse_f64(columns[1], line));
            beta = Some(parse_f64(columns[2], line));
            darad = Some(parse_f64(columns[3], line));
            daeps = Some(parse_f64(columns[4], line));
            continue;
        }
        assert!(columns.len() >= 6, "malformed MMFFVdW line {}", line_idx + 1);
        let alpha_i = parse_f64(columns[1], line);
        let a_i = parse_f64(columns[3], line);
        rows.push(MmffVdwRow {
            atom_type: parse_u8(columns[0], line),
            alpha_i,
            n_i: parse_f64(columns[2], line),
            a_i,
            g_i: parse_f64(columns[4], line),
            r_star: a_i * alpha_i.powf(power.expect("missing generated MMFF VdW power constant")),
            da: parse_vdw_da(columns[5], line),
        });
    }
    MmffVdwData {
        power: power.expect("missing generated MMFF VdW power constant"),
        b: b.expect("missing generated MMFF VdW B constant"),
        beta: beta.expect("missing generated MMFF VdW Beta constant"),
        darad: darad.expect("missing generated MMFF VdW DARAD constant"),
        daeps: daeps.expect("missing generated MMFF VdW DAEPS constant"),
        rows,
    }
}

fn parse_vdw_da(value: &str, line: &str) -> u8 {
    match value {
        "A" => b'A',
        "D" => b'D',
        "-" => b'-',
        other => panic!("invalid generated VdW DA token {other:?} in {line:?}"),
    }
}

fn write_mmff_generated_tables(
    output: &mut String,
    def_rows: &[MmffDefRow],
    prop_rows: &[MmffPropRow],
    pbci_rows: &[MmffPbciRow],
    chg_rows: &[MmffChgRow],
    bond_rows: &[MmffBondRow],
    bndk_rows: &[MmffBndkRow],
    herschbach_laurie_rows: &[MmffHerschbachLaurieRow],
    cov_rad_pau_ele_rows: &[MmffCovRadPauEleRow],
    angle_rows: &[MmffAngleRow],
    stbn_rows: &[MmffStbnRow],
    dfsb_rows: &[MmffDfsbRow],
    oop_rows: &[MmffOopRow],
    mffs_oop_rows: &[MmffOopRow],
    tor_rows: &[MmffTorRow],
    mffs_tor_rows: &[MmffTorRow],
    vdw_data: &MmffVdwData,
) {
    writeln!(output, "pub(crate) static DEFAULT_MMFF_DEF_ROWS: &[(u8, MmffDef)] = &[").unwrap();
    for row in def_rows {
        writeln!(
            output,
            "    ({}, MmffDef {{ eq_level: [{}, {}, {}, {}] }}),",
            row.atom_type, row.eq_level[0], row.eq_level[1], row.eq_level[2], row.eq_level[3]
        )
        .unwrap();
    }
    writeln!(output, "];\n").unwrap();

    writeln!(
        output,
        "pub(crate) static DEFAULT_MMFF_BNDK_ROWS: &[(u8, u8, MmffBond)] = &["
    )
    .unwrap();
    for row in bndk_rows {
        writeln!(
            output,
            "    ({}, {}, MmffBond {{ kb: {:?}, r0: {:?} }}),",
            row.i_atomic_num, row.j_atomic_num, row.kb, row.r0
        )
        .unwrap();
    }
    writeln!(output, "];\n").unwrap();

    writeln!(
        output,
        "pub(crate) static DEFAULT_MMFF_HERSCHBACH_LAURIE_ROWS: &[(u8, u8, MmffHerschbachLaurie)] = &["
    )
    .unwrap();
    for row in herschbach_laurie_rows {
        writeln!(
            output,
            "    ({}, {}, MmffHerschbachLaurie {{ a_ij: {:?}, d_ij: {:?}, dp_ij: {:?} }}),",
            row.i_row, row.j_row, row.a_ij, row.d_ij, row.dp_ij
        )
        .unwrap();
    }
    writeln!(output, "];\n").unwrap();

    writeln!(
        output,
        "pub(crate) static DEFAULT_MMFF_COV_RAD_PAU_ELE_ROWS: &[(u8, MmffCovRadPauEle)] = &["
    )
    .unwrap();
    for row in cov_rad_pau_ele_rows {
        writeln!(
            output,
            "    ({}, MmffCovRadPauEle {{ r0: {:?}, chi: {:?} }}),",
            row.atomic_num, row.r0, row.chi
        )
        .unwrap();
    }
    writeln!(output, "];\n").unwrap();

    writeln!(
        output,
        "pub(crate) static DEFAULT_MMFF_PROP_ROWS: &[(u8, MmffProp)] = &["
    )
    .unwrap();
    for row in prop_rows {
        writeln!(
            output,
            "    ({}, MmffProp {{ atno: {}, crd: {}, val: {}, pilp: {}, mltb: {}, arom: {}, linh: {}, sbmb: {} }}),",
            row.atom_type,
            row.prop[0],
            row.prop[1],
            row.prop[2],
            row.prop[3],
            row.prop[4],
            row.prop[5],
            row.prop[6],
            row.prop[7],
        )
        .unwrap();
    }
    writeln!(output, "];\n").unwrap();

    writeln!(
        output,
        "pub(crate) static DEFAULT_MMFF_PBCI_ROWS: &[(u8, MmffPbci)] = &["
    )
    .unwrap();
    for row in pbci_rows {
        writeln!(
            output,
            "    ({}, MmffPbci {{ pbci: {:?}, fcadj: {:?} }}),",
            row.atom_type, row.pbci, row.fcadj
        )
        .unwrap();
    }
    writeln!(output, "];\n").unwrap();

    writeln!(
        output,
        "pub(crate) static DEFAULT_MMFF_CHG_ROWS: &[(u8, u8, u8, MmffChg)] = &["
    )
    .unwrap();
    for row in chg_rows {
        writeln!(
            output,
            "    ({}, {}, {}, MmffChg {{ bci: {:?} }}),",
            row.bond_type, row.i_atom_type, row.j_atom_type, row.bci
        )
        .unwrap();
    }
    writeln!(output, "];\n").unwrap();

    writeln!(
        output,
        "pub(crate) static DEFAULT_MMFF_BOND_ROWS: &[(u8, u8, u8, MmffBond)] = &["
    )
    .unwrap();
    for row in bond_rows {
        writeln!(
            output,
            "    ({}, {}, {}, MmffBond {{ kb: {:?}, r0: {:?} }}),",
            row.bond_type, row.i_atom_type, row.j_atom_type, row.kb, row.r0
        )
        .unwrap();
    }
    writeln!(output, "];\n").unwrap();

    writeln!(
        output,
        "pub(crate) static DEFAULT_MMFF_ANGLE_ROWS: &[(u8, u8, u8, u8, MmffAngle)] = &["
    )
    .unwrap();
    for row in angle_rows {
        writeln!(
            output,
            "    ({}, {}, {}, {}, MmffAngle {{ ka: {:?}, theta0: {:?} }}),",
            row.angle_type, row.i_atom_type, row.j_atom_type, row.k_atom_type, row.ka, row.theta0
        )
        .unwrap();
    }
    writeln!(output, "];\n").unwrap();

    writeln!(
        output,
        "pub(crate) static DEFAULT_MMFF_STBN_ROWS: &[(u8, u8, u8, u8, MmffStbn)] = &["
    )
    .unwrap();
    for row in stbn_rows {
        writeln!(
            output,
            "    ({}, {}, {}, {}, MmffStbn {{ kba_ijk: {:?}, kba_kji: {:?} }}),",
            row.stretch_bend_type, row.i_atom_type, row.j_atom_type, row.k_atom_type, row.kba_ijk, row.kba_kji
        )
        .unwrap();
    }
    writeln!(output, "];\n").unwrap();

    writeln!(
        output,
        "pub(crate) static DEFAULT_MMFF_DFSB_ROWS: &[(u32, u32, u32, MmffStbn)] = &["
    )
    .unwrap();
    for row in dfsb_rows {
        writeln!(
            output,
            "    ({}, {}, {}, MmffStbn {{ kba_ijk: {:?}, kba_kji: {:?} }}),",
            row.i_atomic_num, row.j_atomic_num, row.k_atomic_num, row.kba_ijk, row.kba_kji
        )
        .unwrap();
    }
    writeln!(output, "];\n").unwrap();

    writeln!(
        output,
        "pub(crate) static DEFAULT_MMFF_OOP_ROWS: &[(u8, u8, u8, u8, MmffOop)] = &["
    )
    .unwrap();
    for row in oop_rows {
        writeln!(
            output,
            "    ({}, {}, {}, {}, MmffOop {{ koop: {:?} }}),",
            row.i_atom_type, row.j_atom_type, row.k_atom_type, row.l_atom_type, row.koop
        )
        .unwrap();
    }
    writeln!(output, "];\n").unwrap();

    writeln!(
        output,
        "pub(crate) static DEFAULT_MMFFS_OOP_ROWS: &[(u8, u8, u8, u8, MmffOop)] = &["
    )
    .unwrap();
    for row in mffs_oop_rows {
        writeln!(
            output,
            "    ({}, {}, {}, {}, MmffOop {{ koop: {:?} }}),",
            row.i_atom_type, row.j_atom_type, row.k_atom_type, row.l_atom_type, row.koop
        )
        .unwrap();
    }
    writeln!(output, "];\n").unwrap();

    writeln!(
        output,
        "pub(crate) static DEFAULT_MMFF_TOR_ROWS: &[(u8, u8, u8, u8, u8, MmffTor)] = &["
    )
    .unwrap();
    for row in tor_rows {
        writeln!(
            output,
            "    ({}, {}, {}, {}, {}, MmffTor {{ v1: {:?}, v2: {:?}, v3: {:?} }}),",
            row.tor_type, row.i_atom_type, row.j_atom_type, row.k_atom_type, row.l_atom_type, row.v1, row.v2, row.v3
        )
        .unwrap();
    }
    writeln!(output, "];\n").unwrap();

    writeln!(
        output,
        "pub(crate) static DEFAULT_MMFFS_TOR_ROWS: &[(u8, u8, u8, u8, u8, MmffTor)] = &["
    )
    .unwrap();
    for row in mffs_tor_rows {
        writeln!(
            output,
            "    ({}, {}, {}, {}, {}, MmffTor {{ v1: {:?}, v2: {:?}, v3: {:?} }}),",
            row.tor_type, row.i_atom_type, row.j_atom_type, row.k_atom_type, row.l_atom_type, row.v1, row.v2, row.v3
        )
        .unwrap();
    }
    writeln!(output, "];\n").unwrap();

    writeln!(
        output,
        "pub(crate) const DEFAULT_MMFF_VDW_POWER: f64 = {:?};",
        vdw_data.power
    )
    .unwrap();
    writeln!(output, "pub(crate) const DEFAULT_MMFF_VDW_B: f64 = {:?};", vdw_data.b).unwrap();
    writeln!(
        output,
        "pub(crate) const DEFAULT_MMFF_VDW_BETA: f64 = {:?};",
        vdw_data.beta
    )
    .unwrap();
    writeln!(
        output,
        "pub(crate) const DEFAULT_MMFF_VDW_DARAD: f64 = {:?};",
        vdw_data.darad
    )
    .unwrap();
    writeln!(
        output,
        "pub(crate) const DEFAULT_MMFF_VDW_DAEPS: f64 = {:?};\n",
        vdw_data.daeps
    )
    .unwrap();
    writeln!(output, "pub(crate) static DEFAULT_MMFF_VDW_ROWS: &[(u8, MmffVdw)] = &[").unwrap();
    for row in &vdw_data.rows {
        writeln!(
            output,
            "    ({}, MmffVdw {{ alpha_i: {:?}, n_i: {:?}, a_i: {:?}, g_i: {:?}, r_star: {:?}, da: {} }}),",
            row.atom_type, row.alpha_i, row.n_i, row.a_i, row.g_i, row.r_star, row.da
        )
        .unwrap();
    }
    writeln!(output, "];\n").unwrap();
}

#[derive(Debug)]
struct UffAtomicParamRow {
    label: String,
    r1: f64,
    theta0_degrees: f64,
    x1: f64,
    d1: f64,
    zeta: f64,
    z1: f64,
    v1: f64,
    u1: f64,
    gmp_xi: f64,
    gmp_hardness: f64,
    gmp_radius: f64,
}

fn parse_uff_atomic_params(param_data: &str) -> Vec<UffAtomicParamRow> {
    param_data
        .lines()
        .filter(|line| !line.starts_with('#'))
        .enumerate()
        .map(|(idx, line)| {
            let columns: Vec<&str> = line.split('\t').collect();
            assert_eq!(columns.len(), 12, "malformed generated UFF row {}: {line}", idx + 1);
            UffAtomicParamRow {
                label: columns[0].to_owned(),
                r1: parse_f64(columns[1], line),
                theta0_degrees: parse_f64(columns[2], line),
                x1: parse_f64(columns[3], line),
                d1: parse_f64(columns[4], line),
                zeta: parse_f64(columns[5], line),
                z1: parse_f64(columns[6], line),
                v1: parse_f64(columns[7], line),
                u1: parse_f64(columns[8], line),
                gmp_xi: parse_f64(columns[9], line),
                gmp_hardness: parse_f64(columns[10], line),
                gmp_radius: parse_f64(columns[11], line),
            }
        })
        .collect()
}

fn parse_f64(value: &str, line: &str) -> f64 {
    value
        .parse::<f64>()
        .unwrap_or_else(|err| panic!("invalid generated table float {value:?} in {line:?}: {err}"))
}

fn parse_u8(value: &str, line: &str) -> u8 {
    let parsed = value
        .parse::<u32>()
        .unwrap_or_else(|err| panic!("invalid generated table integer {value:?} in {line:?}: {err}"));
    u8::try_from(parsed).unwrap_or_else(|_| panic!("generated table integer {value:?} in {line:?} exceeds u8"))
}

fn parse_u32(value: &str, line: &str) -> u32 {
    value
        .parse::<u32>()
        .unwrap_or_else(|err| panic!("invalid generated table integer {value:?} in {line:?}: {err}"))
}

#[derive(Debug)]
struct RdkitIsotopeMassRow {
    atomic_number: u8,
    isotope: u16,
    mass: f64,
}

fn parse_rdkit_isotope_masses(table: &str) -> Vec<RdkitIsotopeMassRow> {
    table
        .lines()
        .filter_map(|raw_line| {
            let line = raw_line.trim();
            if line.is_empty() || line == "EOS" {
                return None;
            }
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < 4 {
                panic!("invalid RDKit isotope data row {line:?}");
            }
            let atomic_number = parse_u8(fields[0], line);
            let isotope = fields[2]
                .parse::<u16>()
                .unwrap_or_else(|err| panic!("invalid generated isotope integer {:?} in {:?}: {err}", fields[2], line));
            let mass = parse_f64(fields[3], line);
            Some(RdkitIsotopeMassRow {
                atomic_number,
                isotope,
                mass,
            })
        })
        .collect()
}

fn parse_rdkit_most_common_isotopes(table: &str) -> Vec<i16> {
    let mut isotopes = Vec::new();
    for raw_line in table.lines() {
        let line = raw_line.trim();
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 10 {
            panic!("invalid RDKit periodic table row {line:?}");
        }
        let atomic_number = usize::from(parse_u8(fields[0], line));
        if atomic_number != isotopes.len() {
            continue;
        }
        let isotope = fields[8].parse::<i16>().unwrap_or_else(|err| {
            panic!(
                "invalid generated most-common isotope integer {:?} in {:?}: {err}",
                fields[8], line
            )
        });
        isotopes.push(isotope);
    }
    assert_eq!(
        isotopes.len(),
        119,
        "RDKit periodicTableAtomData should produce atomic numbers 0..=118"
    );
    isotopes
}

fn parse_rdkit_most_common_isotope_masses(table: &str) -> Vec<f64> {
    let mut masses = Vec::new();
    for raw_line in table.lines() {
        let line = raw_line.trim();
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 10 {
            panic!("invalid RDKit periodic table row {line:?}");
        }
        let atomic_number = usize::from(parse_u8(fields[0], line));
        if atomic_number != masses.len() {
            continue;
        }
        masses.push(parse_f64(fields[9], line));
    }
    assert_eq!(
        masses.len(),
        119,
        "RDKit periodicTableAtomData should produce atomic numbers 0..=118"
    );
    masses
}

fn mmff_symbol_to_rust(symbol: &str) -> &'static str {
    match symbol {
        "defaultMMFFDef" => "DEFAULT_MMFF_DEF",
        "defaultMMFFProp" => "DEFAULT_MMFF_PROP",
        "defaultMMFFPBCI" => "DEFAULT_MMFF_PBCI",
        "defaultMMFFChg" => "DEFAULT_MMFF_CHG",
        "defaultMMFFBond" => "DEFAULT_MMFF_BOND",
        "defaultMMFFBndk" => "DEFAULT_MMFF_BNDK",
        "defaultMMFFHerschbachLaurie" => "DEFAULT_MMFF_HERSCHBACH_LAURIE",
        "defaultMMFFCovRadPauEle" => "DEFAULT_MMFF_COV_RAD_PAU_ELE",
        "defaultMMFFStbn" => "DEFAULT_MMFF_STBN",
        "defaultMMFFDfsb" => "DEFAULT_MMFF_DFSB",
        "defaultMMFFOop" => "DEFAULT_MMFF_OOP",
        "defaultMMFFsOop" => "DEFAULT_MMFFS_OOP",
        "defaultMMFFTor" => "DEFAULT_MMFF_TOR",
        "defaultMMFFsTor" => "DEFAULT_MMFFS_TOR",
        "defaultMMFFVdW" => "DEFAULT_MMFF_VDW",
        _ => panic!("unknown MMFF default symbol {symbol}"),
    }
}

fn write_string_const(output: &mut String, name: &str, value: &str) {
    writeln!(output, "pub(crate) const {name}: &str = {:?};", value).unwrap();
    writeln!(output).unwrap();
}

fn generated_header() -> String {
    "// @generated by crates/cosmolkit-core/build.rs; do not edit by hand.\n\n".to_owned()
}

fn extract_cpp_string(source: &str, symbol: &str) -> Option<String> {
    let header = format!("const std::string {symbol} =");
    extract_cpp_string_by_header(source, &header)
}

fn extract_cpp_string_array(source: &str, symbol: &str, element_separator: &str) -> Option<String> {
    let header = format!("const std::string {symbol}[] =");
    extract_cpp_string_by_header_with_separator(source, &header, element_separator)
}

fn extract_cpp_string_by_header(source: &str, header: &str) -> Option<String> {
    extract_cpp_string_by_header_with_separator(source, header, "")
}

fn extract_cpp_string_by_header_with_separator(source: &str, header: &str, element_separator: &str) -> Option<String> {
    let start = source.find(header)?;
    let after_equals = &source[start + header.len()..];
    let end = find_initializer_end(after_equals)?;
    Some(decode_cpp_string_literals(&after_equals[..end], element_separator))
}

fn find_initializer_end(source: &str) -> Option<usize> {
    let mut in_string = false;
    let mut escaped = false;
    let mut chars = source.char_indices().peekable();
    while let Some((idx, ch)) = chars.next() {
        if in_string {
            if escaped {
                escaped = false;
                continue;
            }
            match ch {
                '\\' => escaped = true,
                '"' => in_string = false,
                _ => {}
            }
            continue;
        }
        if ch == '/' && chars.peek().is_some_and(|(_, next)| *next == '/') {
            chars.next();
            for (_, comment_ch) in chars.by_ref() {
                if comment_ch == '\n' {
                    break;
                }
            }
            continue;
        }
        if ch == '/' && chars.peek().is_some_and(|(_, next)| *next == '*') {
            chars.next();
            let mut prev = '\0';
            for (_, comment_ch) in chars.by_ref() {
                if prev == '*' && comment_ch == '/' {
                    break;
                }
                prev = comment_ch;
            }
            continue;
        }
        if ch == 'R' && chars.peek().is_some_and(|(_, next)| *next == '"') {
            chars.next();
            let mut delimiter = String::new();
            for (_, delimiter_ch) in chars.by_ref() {
                if delimiter_ch == '(' {
                    break;
                }
                delimiter.push(delimiter_ch);
            }
            let terminator = format!("){delimiter}\"");
            let mut raw_tail = String::new();
            for (_, raw_ch) in chars.by_ref() {
                raw_tail.push(raw_ch);
                if raw_tail.ends_with(&terminator) {
                    break;
                }
            }
            continue;
        }
        match ch {
            '"' => in_string = true,
            ';' => return Some(idx),
            _ => {}
        }
    }
    None
}

fn decode_cpp_string_literals(source: &str, element_separator: &str) -> String {
    let mut decoded = String::new();
    let mut chars = source.chars().peekable();
    while let Some(ch) = chars.next() {
        if ch == '/' && chars.peek() == Some(&'/') {
            for comment_ch in chars.by_ref() {
                if comment_ch == '\n' {
                    break;
                }
            }
            continue;
        }
        if ch == '/' && chars.peek() == Some(&'*') {
            chars.next();
            let mut prev = '\0';
            for comment_ch in chars.by_ref() {
                if prev == '*' && comment_ch == '/' {
                    break;
                }
                prev = comment_ch;
            }
            continue;
        }
        if ch == ',' {
            decoded.push_str(element_separator);
            continue;
        }
        if ch == 'R' && chars.peek() == Some(&'"') {
            chars.next();
            let mut delimiter = String::new();
            for delimiter_ch in chars.by_ref() {
                if delimiter_ch == '(' {
                    break;
                }
                delimiter.push(delimiter_ch);
            }
            let terminator = format!("){delimiter}\"");
            let mut raw = String::new();
            for raw_ch in chars.by_ref() {
                raw.push(raw_ch);
                if raw.ends_with(&terminator) {
                    let content_len = raw.len() - terminator.len();
                    decoded.push_str(&raw[..content_len]);
                    break;
                }
            }
            continue;
        }
        if ch != '"' {
            continue;
        }
        while let Some(literal_ch) = chars.next() {
            match literal_ch {
                '"' => break,
                '\\' => {
                    let escaped = chars.next().expect("unterminated C++ string escape");
                    match escaped {
                        'n' => decoded.push('\n'),
                        't' => decoded.push('\t'),
                        '"' => decoded.push('"'),
                        '\\' => decoded.push('\\'),
                        other => panic!("unsupported C++ string escape \\{other}"),
                    }
                }
                other => decoded.push(other),
            }
        }
    }
    decoded
}
