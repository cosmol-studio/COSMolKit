use std::io;

use cosmolkit_inchi::inchi_to_inchi_key;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let inchi = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "InChI=1S/CH4/h1H4".to_owned());
    let output = inchi_to_inchi_key(inchi.as_bytes())?;

    for diagnostic in &output.diagnostics {
        eprintln!("{:?}: {}", diagnostic.level, diagnostic.message.trim_end());
    }
    if output.key.is_empty() {
        return Err(io::Error::other("InChI did not produce an InChIKey").into());
    }

    println!("{}", String::from_utf8(output.key)?);
    Ok(())
}
