use crate::{Result, cli::Options, storage::hex};
use sha2::{Digest, Sha256};
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

pub const MAX_INPUT_BYTES: usize = 1024 * 1024;

#[derive(Clone)]
pub struct Corpus {
    pub path: PathBuf,
    pub sha256: String,
    pub records: u64,
    pub offset: usize,
    pub limit: Option<usize>,
}

impl Corpus {
    pub fn open(options: &Options) -> Result<Self> {
        let path =
            if options.corpus.components().count() == 1 && options.corpus.extension().is_none() {
                Path::new(env!("CARGO_MANIFEST_DIR"))
                    .join("testdata/corpora")
                    .join(&options.corpus)
                    .with_extension("smi")
            } else {
                let path = if options.corpus.is_absolute() || options.corpus.exists() {
                    options.corpus.clone()
                } else {
                    Path::new(env!("CARGO_MANIFEST_DIR"))
                        .parent()
                        .unwrap()
                        .join(&options.corpus)
                };
                std::fs::canonicalize(path)?
            };
        let mut digest = Sha256::new();
        let mut records = 0u64;
        digest.update((options.offset as u64).to_le_bytes());
        for row in Self::rows_at(&path, options.offset, options.limit)? {
            let row = row?;
            digest.update((row.len() as u64).to_le_bytes());
            digest.update(row.as_bytes());
            records += 1;
        }
        if records == 0 {
            return Err("selected corpus is empty".into());
        }
        Ok(Self {
            path,
            sha256: hex(&digest.finalize()),
            records,
            offset: options.offset,
            limit: options.limit,
        })
    }

    pub fn rows(&self) -> Result<impl Iterator<Item = Result<String>>> {
        Self::rows_at(&self.path, self.offset, self.limit)
    }

    pub fn structured(&self) -> bool {
        self.path.extension().is_some_and(|ext| ext == "jsonl")
    }

    fn rows_at(
        path: &Path,
        offset: usize,
        limit: Option<usize>,
    ) -> Result<impl Iterator<Item = Result<String>>> {
        let mut reader = BufReader::new(File::open(path)?);
        let structured = path.extension().is_some_and(|ext| ext == "jsonl");
        let max_bytes = if structured {
            64 * MAX_INPUT_BYTES
        } else {
            MAX_INPUT_BYTES
        };
        let mut line = Vec::new();
        let rows = std::iter::from_fn(move || {
            loop {
                line.clear();
                match (&mut reader)
                    .take(max_bytes as u64 + 1)
                    .read_until(b'\n', &mut line)
                {
                    Ok(0) => return None,
                    Ok(_) if line.len() > max_bytes => {
                        return Some(Err("input line exceeds its bounded row budget".into()));
                    }
                    Ok(_) => {
                        let text = match std::str::from_utf8(&line) {
                            Ok(text) => text.trim(),
                            Err(error) => return Some(Err(error.into())),
                        };
                        if text.is_empty() || text.starts_with('#') {
                            continue;
                        }
                        return Some(Ok(if structured {
                            text.to_owned()
                        } else {
                            text.split_whitespace().next().unwrap().to_owned()
                        }));
                    }
                    Err(error) => return Some(Err(error.into())),
                }
            }
        });
        Ok(rows.skip(offset).take(limit.unwrap_or(usize::MAX)))
    }
}

use std::io::Read;

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn structured_rows_preserve_spaces_and_hash_the_actual_text() {
        let mut file = tempfile::Builder::new()
            .suffix(".jsonl")
            .tempfile()
            .unwrap();
        writeln!(file, "{{\"text\": \"ATOM has spaces\"}}").unwrap();
        let options =
            Options::parse(["--corpus".to_owned(), file.path().display().to_string()]).unwrap();
        let first = Corpus::open(&options).unwrap();
        assert!(first.structured());
        assert_eq!(
            first.rows().unwrap().next().unwrap().unwrap(),
            "{\"text\": \"ATOM has spaces\"}"
        );
        writeln!(file, "{{\"text\": \"second row\"}}").unwrap();
        assert_ne!(first.sha256, Corpus::open(&options).unwrap().sha256);
    }
}
