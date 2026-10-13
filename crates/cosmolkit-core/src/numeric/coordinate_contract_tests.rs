use super::parse_rdkit_atof_prefix;

#[test]
fn v3k_atom_numbers_coordinate_contract_lexical_boundaries() {
    // RDKit ParseV3000AtomBlock uses unscreened atof. Missing fields are a
    // different error at the reader boundary, not a call with an empty token.
    for (text, expected, consumed) in [
        ("-0", -0.0_f64, 2),
        ("+0.000e-12", 0.0, 10),
        ("-0.000E+12", -0.0, 10),
        (" \t\n\r\u{b}\u{c}-1.25tail", -1.25, 11),
        ("+3.125e+2suffix", 312.5, 9),
        ("1e+", 1.0, 1),
        ("-.5E-1", -0.05, 6),
        ("1.25\0ignored", 1.25, 4),
        ("abc", 0.0, 0),
        ("  -.", 0.0, 0),
        ("", 0.0, 0),
    ] {
        let (value, offset) = parse_rdkit_atof_prefix(text.as_bytes());
        assert_eq!(value.to_bits(), expected.to_bits(), "{text:?}");
        assert_eq!(offset, consumed, "{text:?}");
    }
}

#[test]
#[ignore = "explicit live oracle: pinned RDKit wheel and glibc binaries required"]
fn v3k_atom_numbers_coordinate_contract_rdkit_bitwise() {
    use std::io::Write;
    use std::process::{Command, Stdio};

    // Frozen before comparison: 1..=17 mantissa digits, both signs, seven
    // decimal exponents, 20 deterministic samples per combination. This is
    // coverage of ordinary decimal conversion, NOT an input rejection limit
    // or an exhaustive proof for all finite decimal strings.
    let mut tokens: Vec<String> = [
        "0",
        "-0",
        "+0",
        "-0.000E+12",
        "0.1",
        "-0.1",
        ".5",
        "5.",
        "1e+",
        "1e-",
        "1e",
        "12xyz",
        "abc",
        "1,5",
        "1.2345",
        "9999.999",
        "-999.999",
        "99999.9999",
        "-9999.9999",
        "1.0000000000000001",
        "1.0000000000000002",
        "9007199254740993",
        "2.2250738585072014e-308",
        "1.7976931348623157e308",
    ]
    .into_iter()
    .map(str::to_owned)
    .collect();
    let mut state = 2038_u64;
    for digits in 1..=17 {
        for exponent in [-100, -12, -4, 0, 4, 12, 100] {
            for _ in 0..20 {
                let mut mantissa = String::new();
                for index in 0..digits {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    let digit = if index == 0 {
                        1 + (state % 9) as u8
                    } else {
                        (state % 10) as u8
                    };
                    mantissa.push(char::from(b'0' + digit));
                    if index == 0 {
                        mantissa.push('.');
                    }
                }
                for sign in ["+", "-"] {
                    tokens.push(format!("{sign}{mantissa}e{exponent:+}"));
                }
            }
        }
    }
    assert_eq!(tokens.len(), 4784);
    let script = r#"
import ctypes, hashlib, locale, pathlib, struct, sys
from rdkit import Chem, rdBase
assert rdBase.rdkitVersion == '2026.03.1', rdBase.rdkitVersion
locale.setlocale(locale.LC_NUMERIC, 'C')
libc = ctypes.CDLL(None)
assert libc.fegetround() == 0, 'FE_TONEAREST required'
paths = {line.split()[-1] for line in pathlib.Path('/proc/self/maps').read_text().splitlines() if '/' in line}
for suffix, expected in [
    ('/libc.so.6', '85e64f97e348786a8fb4d9f3d52fec289e2fb86bba20f0731dfe61990525e0f7'),
    ('/libRDKitFileParsers-288b044d.so.1', '5005c4b872b61c8bcde1887699a116de571900a2aa6a5d9e3b9c9ead70c6ac1b'),
]:
    found = [p for p in paths if p.endswith(suffix)]
    assert len(found) == 1, (suffix, found)
    assert hashlib.sha256(pathlib.Path(found[0]).read_bytes()).hexdigest() == expected, found[0]
libc.strtod.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_char_p)]
libc.strtod.restype = ctypes.c_double
template = Chem.MolToMolBlock(Chem.MolFromSmiles('C'), forceV3000=True)
tokens = sys.stdin.read().splitlines()
for token in tokens:
    block = '\n'.join(('M  V30 1 C ' + token + ' 0 0 0') if line.startswith('M  V30 1 C ') else line for line in template.split('\n'))
    mol = Chem.MolFromMolBlock(block, sanitize=False, removeHs=False, strictParsing=True)
    assert mol is not None, token
    value = mol.GetConformer().GetAtomPosition(0).x
    buf = ctypes.create_string_buffer(token.encode('ascii'))
    end = ctypes.c_char_p()
    raw = libc.strtod(buf, ctypes.byref(end))
    assert struct.pack('>d', raw) == struct.pack('>d', value), token
    offset = ctypes.cast(end, ctypes.c_void_p).value - ctypes.addressof(buf)
    print(struct.pack('>d', value).hex(), offset)
"#;
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let mut child = Command::new(root.join(".venv/bin/python"))
        .args(["-B", "-c", script])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("start pinned RDKit coordinate oracle");
    {
        let mut input = child.stdin.take().expect("oracle stdin");
        for token in &tokens {
            writeln!(input, "{token}").expect("write token");
        }
    }
    let output = child.wait_with_output().expect("wait for oracle");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).expect("oracle UTF-8");
    let rows: Vec<_> = stdout.lines().collect();
    assert_eq!(rows.len(), tokens.len(), "all oracle rows required");
    for (token, row) in tokens.iter().zip(rows) {
        let (bits, offset) = row.split_once(' ').expect("bits and offset");
        let expected = u64::from_str_radix(bits, 16).expect("binary64 bits");
        let expected_offset: usize = offset.parse().expect("consumed bytes");
        let (actual, actual_offset) = parse_rdkit_atof_prefix(token.as_bytes());
        assert_eq!(actual.to_bits(), expected, "RDKit bit mismatch: {token}");
        assert_eq!(
            actual_offset, expected_offset,
            "RDKit offset mismatch: {token}"
        );
    }
    eprintln!(
        "RDKit coordinate profile: {} rows, exact bits and offsets",
        tokens.len()
    );
}
