# RDKit MMFFFormalChargeCollection Absence

## Source Verdict

`ForceFields::MMFF::MMFFFormalChargeCollection` is not a symbol in the
vendored RDKit MMFF parameter source.

Checked source files:

- `third_party/rdkit/Code/ForceField/MMFF/Params.h`
- `third_party/rdkit/Code/ForceField/MMFF/Params.cpp`
- `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h`
- `third_party/rdkit/Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp`

Observed RDKit formal-charge handling is attached to MMFF atom properties in
`AtomTyper`, not to a parameter collection in `Params.h` or `Params.cpp`.

Relevant source-backed symbols:

- `MMFFMolProperties::getMMFFFormalCharge`
- `MMFFMolProperties::setMMFFFormalCharge`
- `MMFFAtomProperties::mmffFormalCharge`

## Plan Correction

The original plan entries that requested porting
`ForceFields::MMFF::MMFFFormalChargeCollection` were invalid because the
referenced RDKit collection does not exist. The corrected plan records this
absence explicitly and adds a source-inventory regression test so future porting
work does not fabricate an unsupported parameter collection.

Future MMFF formal-charge porting must use `AtomTyper.h` and `AtomTyper.cpp`
as the source anchors.
