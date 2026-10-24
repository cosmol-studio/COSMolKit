//! Exact shared fingerprint error graph, retaining non-Clone native causes.
use crate::alignment_errors::ThreadCountError;
use crate::allocation_errors::TryReserveError;
use crate::canonical_errors::PropertyStringError;
use crate::cip_atom_code_errors::AtomCodeError;
use crate::error_framework::{Source, error_projection};
use crate::errors::FingerprintError;
use crate::fingerprint_json::FingerprintJsonError;
use crate::io_source_errors::{ErrorRef, IoError};
use crate::match_errors::{QueryMatchContextError, SubstructMatchError};
use crate::matrix_errors::MatrixError;
use crate::model_edit_errors::MoleculePropertyError;
use crate::model_errors::TopologyValidationError;
use crate::path_errors::{PathError, PeriodicTableError};
use crate::retained_error::RetainedError;
use crate::ring_errors::RingFindingError;
use crate::smarts_errors::SmartsParseError;
use crate::smarts_writer::SmartsWriteError;
use crate::stereo_errors::LegacyStereoError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

// A real source chain can alternate Morgan and atom-pair errors. Each wrapper
// borrows its exact node from one retained native root, never rebuilding it.
fn project(node: ErrorRef) -> PyErr {
    macro_rules! retained {
        ($($name:ident),+ $(,)?) => { $(
            if node.error().is::<ck::$name>() { return $name::from_ref(node); }
        )+ };
    }
    retained!(
        MorganReadError,
        AtomPairReadError,
        TopologicalTorsionReadError,
        FingerprintPreparationError,
        MorganError,
        AtomPairError,
        TopologicalTorsionError,
        FingerprintWorkerError,
        FingerprintJsonError,
        SmartsWriteError,
        ThreadCountError
    );
    if node.error().is::<std::io::Error>() {
        return IoError::from_ref(node);
    }
    if let Some(value) = node
        .error()
        .downcast_ref::<std::collections::TryReserveError>()
    {
        return TryReserveError::from_rust(value.clone());
    }
    macro_rules! value {
        ($($name:ident),+ $(,)?) => { $(
            if let Some(value) = node.error().downcast_ref::<ck::$name>() {
                return $name::from_rust(value.clone());
            }
        )+ };
    }
    value!(
        FingerprintError,
        MatrixError,
        AtomCodeError,
        TopologyValidationError,
        ValenceError,
        PeriodicTableError,
        LegacyStereoError,
        MoleculePropertyError,
        PropertyStringError,
        SmartsParseError,
        QueryMatchContextError,
        SubstructMatchError,
        PathError,
        RingFindingError
    );
    unreachable!("unmodeled concrete source in the native fingerprint error graph")
}

error_projection!(retained FingerprintPreparationError, FingerprintPreparationErrorKind, ck::FingerprintPreparationError,"fingerprint",
    { MissingPreparedValence, RingPreparation(_) },
    cause=|inner: &RetainedError<ck::FingerprintPreparationError>| inner.node().source().map(project), source=Option<crate::error_framework::ErrorObject<RingFindingError>>);
error_projection!(retained FingerprintWorkerError, FingerprintWorkerErrorKind, ck::FingerprintWorkerError,"fingerprint",
    { ThreadCount(_), ThreadSpawn(_), Panic, Protocol },
    cause=|inner: &RetainedError<ck::FingerprintWorkerError>| inner.node().source().map(project), source=Option<Source<crate::error_framework::ErrorObject<ThreadCountError>,crate::error_framework::ErrorObject<IoError>>>);
type MorganSource = Source<
    Source<
        crate::error_framework::ErrorObject<FingerprintError>,
        Source<
            crate::error_framework::ErrorObject<FingerprintJsonError>,
            Source<
                crate::error_framework::ErrorObject<SmartsWriteError>,
                crate::error_framework::ErrorObject<AtomPairError>,
            >,
        >,
    >,
    Source<
        Source<
            crate::error_framework::ErrorObject<FingerprintWorkerError>,
            Source<
                crate::error_framework::ErrorObject<MatrixError>,
                Source<
                    crate::error_framework::ErrorObject<ValenceError>,
                    crate::error_framework::ErrorObject<PeriodicTableError>,
                >,
            >,
        >,
        Source<
            Source<
                crate::error_framework::ErrorObject<LegacyStereoError>,
                Source<
                    crate::error_framework::ErrorObject<MoleculePropertyError>,
                    Source<
                        crate::error_framework::ErrorObject<PropertyStringError>,
                        crate::error_framework::ErrorObject<SmartsParseError>,
                    >,
                >,
            >,
            Source<
                crate::error_framework::ErrorObject<QueryMatchContextError>,
                crate::error_framework::ErrorObject<SubstructMatchError>,
            >,
        >,
    >,
>;
error_projection!(retained MorganError, MorganErrorKind, ck::MorganError,"fingerprint", {
    Fingerprint(_), Json(_), SmartsWrite(_), AtomPair(_), Worker(_), StatePoisoned,
    Matrix(_), Valence(_), PeriodicTable(_), LegacyStereo(_), MoleculeProperty(_),
    PropertyString(_), SmartsParse(_), QueryMatchContext(_), SubstructMatch(_),
}, cause=|inner: &RetainedError<ck::MorganError>| inner.node().source().map(project), source=Option<MorganSource>);
type AtomPairSource = Source<
    crate::error_framework::ErrorObject<FingerprintError>,
    Source<
        crate::error_framework::ErrorObject<AtomCodeError>,
        Source<
            crate::error_framework::ErrorObject<MatrixError>,
            Source<
                crate::error_framework::ErrorObject<TopologyValidationError>,
                crate::error_framework::ErrorObject<MorganError>,
            >,
        >,
    >,
>;
error_projection!(retained AtomPairError, AtomPairErrorKind, ck::AtomPairError,"fingerprint", {
    Fingerprint(_), AtomCode(_), Matrix(_), Topology(_), Preparation(_),
    AtomInvariantLength { length, atom_index },
}, cause=|inner: &RetainedError<ck::AtomPairError>| inner.node().source().map(project), source=Option<AtomPairSource>, methods={
    #[getter]
    fn length(&self) -> Option<usize> { match self.inner.as_ref() {
        ck::AtomPairError::AtomInvariantLength { length, .. } => Some(*length), _ => None,
    } }
    #[getter]
    fn atom_index(&self) -> Option<usize> { match self.inner.as_ref() {
        ck::AtomPairError::AtomInvariantLength { atom_index, .. } => Some(*atom_index), _ => None,
    } }
});
type TorsionSource = Source<
    Source<
        crate::error_framework::ErrorObject<FingerprintError>,
        Source<
            crate::error_framework::ErrorObject<AtomPairError>,
            Source<
                crate::error_framework::ErrorObject<PathError>,
                crate::error_framework::ErrorObject<MorganError>,
            >,
        >,
    >,
    Source<
        crate::error_framework::ErrorObject<FingerprintJsonError>,
        Source<
            crate::error_framework::ErrorObject<ThreadCountError>,
            Source<
                crate::error_framework::ErrorObject<IoError>,
                crate::error_framework::ErrorObject<TryReserveError>,
            >,
        >,
    >,
>;
error_projection!(retained TopologicalTorsionError, TopologicalTorsionErrorKind, ck::TopologicalTorsionError,"fingerprint", {
    Fingerprint(_), AtomInvariants(_), Path(_), Preparation(_), Json(_), StatePoisoned,
    ThreadCount(_), ThreadSpawn(_), WorkerPanic, WorkerProtocol, OutputAllocation(_),
}, cause=|inner: &RetainedError<ck::TopologicalTorsionError>| inner.node().source().map(project), source=Option<TorsionSource>);
error_projection!(retained MorganReadError, MorganReadErrorKind, ck::MorganReadError,"fingerprint", { Preparation(_), Generator(_) },
    cause=|inner: &RetainedError<ck::MorganReadError>| inner.node().source().map(project), source=Option<Source<crate::error_framework::ErrorObject<FingerprintPreparationError>,crate::error_framework::ErrorObject<MorganError>>>);
error_projection!(retained AtomPairReadError, AtomPairReadErrorKind, ck::AtomPairReadError,"fingerprint", { Preparation(_), Generator(_) },
    cause=|inner: &RetainedError<ck::AtomPairReadError>| inner.node().source().map(project), source=Option<Source<crate::error_framework::ErrorObject<FingerprintPreparationError>,crate::error_framework::ErrorObject<AtomPairError>>>);
error_projection!(retained TopologicalTorsionReadError, TopologicalTorsionReadErrorKind, ck::TopologicalTorsionReadError,"fingerprint", { Preparation(_), Generator(_) },
    cause=|inner: &RetainedError<ck::TopologicalTorsionReadError>| inner.node().source().map(project), source=Option<Source<crate::error_framework::ErrorObject<FingerprintPreparationError>,crate::error_framework::ErrorObject<TopologicalTorsionError>>>);
