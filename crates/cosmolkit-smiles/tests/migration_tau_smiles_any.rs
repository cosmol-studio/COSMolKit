use cosmolkit_smiles::{parse_smiles, write_smiles};
use cosmolkit_types::BondStereo;
#[test]
fn source_stereo_any_serializes_the_writer_copy_without_changing_the_input() {
    for text in [
        "CC=C(O)C",
        "C1=CCCCC1",
        "CC=N",
        "CC=CC=O",
        "O=C(O)C(CO)=NC(=O)C1=C(O)C(O)CC=C1",
        "O=C(O)C(CO)=NC(=O)C1=CCCC(O)=C1O",
    ] {
        let mut record = parse_smiles(text, &Default::default()).unwrap();
        let expected = write_smiles(&record).unwrap();
        for bond in &mut record.topology.bonds {
            if bond.order() == cosmolkit_types::BondOrder::Double {
                bond.set_stereo(BondStereo::Any).unwrap();
            }
        }
        let before = record.clone();
        assert_eq!(write_smiles(&record).unwrap(), expected, "{text}");
        assert_eq!(record, before, "{text}: input changed");
    }
}
