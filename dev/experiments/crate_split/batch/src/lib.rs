use cosmolkit_experiment_model::{CoordinateBlock, MoleculeProperties, TopologyBlock};
use rayon::prelude::*;

#[derive(Clone, Debug, PartialEq)]
pub struct BatchRecord {
    pub index: usize,
    pub topology: TopologyBlock,
    pub coordinates: CoordinateBlock,
    pub properties: MoleculeProperties,
}

pub type BatchProgress<'a> = Option<&'a (dyn Fn() + Sync)>;

pub fn map_parsed<S, R, E>(
    sources: impl IntoIterator<Item = S>,
    n_jobs: Option<usize>,
    progress: BatchProgress<'_>,
    parser: impl Fn(S) -> Result<R, E> + Send + Sync,
) -> Result<Vec<BatchRecord>, String>
where
    S: AsRef<str> + Send,
    R: Into<(TopologyBlock, CoordinateBlock, MoleculeProperties)>,
    E: ToString + Send,
{
    map_ordered(
        sources.into_iter().enumerate(),
        n_jobs,
        progress,
        |(index, source)| {
            let (topology, coordinates, properties) = parser(source)?.into();
            Ok(BatchRecord {
                index,
                topology,
                coordinates,
                properties,
            })
        },
    )
    .map_err(|error: E| error.to_string())
}

pub fn map_ordered<T, R, E>(
    values: impl IntoIterator<Item = T>,
    n_jobs: Option<usize>,
    progress: BatchProgress<'_>,
    operation: impl Fn(T) -> Result<R, E> + Send + Sync,
) -> Result<Vec<R>, E>
where
    T: Send,
    R: Send,
    E: Send,
{
    let values: Vec<T> = values.into_iter().collect();
    let run = move || {
        values
            .into_par_iter()
            .map(|value| {
                let result = operation(value);
                if let Some(progress) = progress {
                    progress();
                }
                result
            })
            .collect()
    };
    match n_jobs.map(|value| value.max(1)) {
        Some(n_jobs) => rayon::ThreadPoolBuilder::new()
            .num_threads(n_jobs)
            .build()
            .expect("experiment batch thread pool build must succeed")
            .install(run),
        None => run(),
    }
}
