use cosmolkit_parity_tests_fixed::testing;

macro_rules! define_tests {
    (composition: $($key:ident),*) => {
        $(#[cfg(parity_corpus_smiles)] #[test] fn $key() { testing::run_corpus(stringify!($key)).unwrap(); })*
    };
    ($(($key:ident, $operation:expr, CorpusType::$corpus:ident, $generator:literal)),* $(,)?) => {
        $(one_test!($corpus, $key);)*
    };
}
macro_rules! one_test {
    (Smiles, $key:ident) => {
        #[cfg(parity_corpus_smiles)]
        #[test]
        fn $key() {
            testing::run_corpus(stringify!($key)).unwrap();
        }
    };
    (Pdb, $key:ident) => {
        #[cfg(parity_corpus_bio)]
        #[test]
        fn $key() {
            testing::run_corpus(stringify!($key)).unwrap();
        }
    };
    (Cif, $key:ident) => {
        #[cfg(parity_corpus_bio)]
        #[test]
        fn $key() {
            testing::run_corpus(stringify!($key)).unwrap();
        }
    };
}
cosmolkit_parity_tests_fixed::corpus_tasks!(define_tests);
