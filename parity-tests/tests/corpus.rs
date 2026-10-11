fn main() {
    if let Err(error) = cosmolkit_parity_tests::runner::main() {
        eprintln!("{error}");
        std::process::exit(1);
    }
}
