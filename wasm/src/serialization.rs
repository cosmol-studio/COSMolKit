//! Binary projection delegates to the public molecule archive boundary.
use crate::Molecule;
use cosmolkit as ck;
use std::cell::RefCell;
impl Molecule {
    pub fn to_binary(&self) -> Result<Vec<u8>, ck::PickleError> {
        // COSMolKit❗✔️: self.inner.to_binary()
        self.inner.borrow().to_binary()
    }
    pub fn from_binary(data: &[u8]) -> Result<Self, ck::PickleError> {
        // COSMolKit❗✔️: ck::Molecule::from_binary(data)
        ck::Molecule::from_binary(data).map(|inner| Self {
            inner: RefCell::new(inner),
        })
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn archive_bytes_and_restored_state_match_owner_without_mutation() {
        for text in [
            "",
            "CCO",
            "c1ccccc1",
            "F[C@H](Cl)Br",
            "[13CH3:7][NH3+]",
            "C.O",
        ] {
            let m = Molecule::from_smiles(text).unwrap();
            let owner = ck::Molecule::from_smiles(text).unwrap();
            let before = m.inner.borrow().clone();
            let bytes = m.to_binary().unwrap();
            assert_eq!(bytes, owner.to_binary().unwrap());
            let restored = Molecule::from_binary(&bytes).unwrap();
            assert_eq!(
                *restored.inner.borrow(),
                ck::Molecule::from_binary(&bytes).unwrap()
            );
            assert_eq!(*restored.inner.borrow(), before);
            assert_eq!(restored.to_binary().unwrap(), bytes);
            assert_eq!(*m.inner.borrow(), before);
            for end in 0..bytes.len() {
                let data = &bytes[..end];
                assert_eq!(
                    Molecule::from_binary(data).err(),
                    ck::Molecule::from_binary(data).err(),
                    "cut {end}"
                );
                assert!(Molecule::from_binary(data).is_err());
            }
            for (offset, value) in [(8, 99), (16, 99), (19, 255), (18, 0), (14, 99)] {
                let mut damaged = bytes.clone();
                damaged[offset] = value;
                assert_eq!(
                    Molecule::from_binary(&damaged).err(),
                    ck::Molecule::from_binary(&damaged).err()
                );
                assert!(Molecule::from_binary(&damaged).is_err());
            }
        }
        let owner =
            ck::Molecule::from_xyz_block("2\narchive\nC -0.0 1.25 -2.5\nO 3 4 5\n").unwrap();
        let m = Molecule {
            inner: RefCell::new(owner.clone()),
        };
        let bytes = m.to_binary().unwrap();
        assert_eq!(bytes, owner.to_binary().unwrap());
        let restored = Molecule::from_binary(&bytes).unwrap();
        assert_eq!(*restored.inner.borrow(), owner);
        assert_eq!(restored.to_binary().unwrap(), bytes);
        assert_eq!(
            Molecule::from_binary(&[]).err(),
            Some(ck::PickleError::UnexpectedEof)
        );
        assert_eq!(
            Molecule::from_binary(&[255]).err(),
            Some(ck::PickleError::UnsupportedVersion(255))
        );
    }
}
