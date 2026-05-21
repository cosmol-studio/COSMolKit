use cosmolkit_core::Molecule;

fn main() {
    let smiles = "[#6]-[#8]-[#6]-1-[#6](-[#8])-[#6]-[#6]C2([#6]-[#8]2)[#6]-1[C@@]1([#6])[#8]-[#6]1-[#6]\\[#6]=[#6](\\[#6])-[#6]";
    let mol = Molecule::from_smiles(smiles).expect("parse row 123 smiles");
    let svg = mol.to_svg(300, 300).expect("draw row 123 svg");
    println!("{svg}");
}
