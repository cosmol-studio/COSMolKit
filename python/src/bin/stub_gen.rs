use pyo3_stub_gen::Result;
use std::fs;
use std::io;
use std::path::Path;

fn main() -> Result<()> {
    let stub = cosmolkit::stub_info()?;
    stub.generate()?;

    let pyi_path = [
        Path::new("./cosmolkit.pyi"),
        Path::new("./python/cosmolkit.pyi"),
    ]
    .into_iter()
    .find(|path| path.exists())
    .ok_or_else(|| io::Error::new(io::ErrorKind::NotFound, "cosmolkit.pyi was not generated"))?;
    let mut text = fs::read_to_string(pyi_path)?;
    let future_line = "from __future__ import annotations\n";
    if !text.contains(future_line) {
        text = format!("{future_line}{text}");
    }
    text = expose_batch_validation_error(text);
    fs::write(pyi_path, text)?;

    Ok(())
}

fn expose_batch_validation_error(mut text: String) -> String {
    let export = "    \"BatchValidationError\",\n";
    if !text.contains(export) {
        text = text.replace(
            "    \"BatchExportReport\",\n",
            &format!("    \"BatchExportReport\",\n{export}"),
        );
    }

    let class_decl = "class BatchValidationError(builtins.ValueError):\n    ...\n\n";
    if !text.contains("class BatchValidationError") {
        text = text.replace(
            "@typing.final\nclass Alignment:",
            &format!("{class_decl}@typing.final\nclass Alignment:"),
        );
    }

    text
}
