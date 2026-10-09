//! Export the actual compiled API contract for development binding checks.
#[path = "support/binding_contract_manifest.rs"]
mod contract;

fn main() {
    println!("{}", contract::manifest());
}
