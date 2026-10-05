/// Completion state of a tautomer-enumeration run.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TautomerEnumerationStatus {
    #[default]
    Completed,
    MaxTautomersReached,
    MaxTransformsReached,
    Canceled,
}

/// Source-compatible limits and stereochemistry policies for enumeration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TautomerParams {
    max_tautomers: u32,
    max_transforms: u32,
    remove_sp3_stereo: bool,
    remove_bond_stereo: bool,
    remove_isotopic_hydrogens: bool,
    reassign_stereo: bool,
}

impl Default for TautomerParams {
    fn default() -> Self {
        Self {
            max_tautomers: 1000,
            max_transforms: 1000,
            remove_sp3_stereo: true,
            remove_bond_stereo: true,
            remove_isotopic_hydrogens: true,
            reassign_stereo: true,
        }
    }
}

impl TautomerParams {
    #[must_use]
    pub const fn max_tautomers(self) -> u32 {
        self.max_tautomers
    }

    pub fn set_max_tautomers(&mut self, value: u32) {
        self.max_tautomers = value;
    }

    #[must_use]
    pub const fn with_max_tautomers(mut self, value: u32) -> Self {
        self.max_tautomers = value;
        self
    }

    #[must_use]
    pub const fn max_transforms(self) -> u32 {
        self.max_transforms
    }

    pub fn set_max_transforms(&mut self, value: u32) {
        self.max_transforms = value;
    }

    #[must_use]
    pub const fn with_max_transforms(mut self, value: u32) -> Self {
        self.max_transforms = value;
        self
    }

    #[must_use]
    pub const fn remove_sp3_stereo(self) -> bool {
        self.remove_sp3_stereo
    }

    pub fn set_remove_sp3_stereo(&mut self, value: bool) {
        self.remove_sp3_stereo = value;
    }

    #[must_use]
    pub const fn with_remove_sp3_stereo(mut self, value: bool) -> Self {
        self.remove_sp3_stereo = value;
        self
    }

    #[must_use]
    pub const fn remove_bond_stereo(self) -> bool {
        self.remove_bond_stereo
    }

    pub fn set_remove_bond_stereo(&mut self, value: bool) {
        self.remove_bond_stereo = value;
    }

    #[must_use]
    pub const fn with_remove_bond_stereo(mut self, value: bool) -> Self {
        self.remove_bond_stereo = value;
        self
    }

    #[must_use]
    pub const fn remove_isotopic_hydrogens(self) -> bool {
        self.remove_isotopic_hydrogens
    }

    pub fn set_remove_isotopic_hydrogens(&mut self, value: bool) {
        self.remove_isotopic_hydrogens = value;
    }

    #[must_use]
    pub const fn with_remove_isotopic_hydrogens(mut self, value: bool) -> Self {
        self.remove_isotopic_hydrogens = value;
        self
    }

    #[must_use]
    pub const fn reassign_stereo(self) -> bool {
        self.reassign_stereo
    }

    pub fn set_reassign_stereo(&mut self, value: bool) {
        self.reassign_stereo = value;
    }

    #[must_use]
    pub const fn with_reassign_stereo(mut self, value: bool) -> Self {
        self.reassign_stereo = value;
        self
    }
}
