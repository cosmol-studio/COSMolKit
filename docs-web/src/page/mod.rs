mod home;
mod placeholders;
mod python;

pub(crate) use home::Home;
pub(crate) use placeholders::{Benchmarks, JavaScript, Validation};
pub(crate) use python::{
    Api, Batch, Confseq, Descriptors, Fingerprints, Genindex, Installation, Io, Molecule, Protein,
    PyModindex, Python, Quickstart, SearchPage,
};
