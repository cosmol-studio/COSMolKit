use pyo3_stub_gen::Result;

fn main() -> Result<()> {
    let stub = cosmolkit::stub_info()?;
    stub.generate()?;
    Ok(())
}
