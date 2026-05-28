# RDKit UFF TorsionAngle Missing `calcTorsionEnergy`

## Scope

This report covers the planned port target:

```text
ForceFields::UFF::TorsionAngleContrib::Utils::calcTorsionEnergy
```

The target was checked against:

```text
third_party/rdkit/Code/ForceField/UFF/TorsionAngle.h
third_party/rdkit/Code/ForceField/UFF/TorsionAngle.cpp
third_party/rdkit/Code/ForceField/UFF/
```

## Finding

RDKit UFF has no `calcTorsionEnergy` source symbol.

The UFF torsion energy behavior is implemented directly in:

```text
TorsionAngleContrib::getEnergy
```

The actual UFF TorsionAngle source-backed symbols are:

```text
ForceFields::UFF::Utils::calculateCosTorsion
ForceFields::UFF::Utils::calcTorsionGrad
ForceFields::UFF::Utils::equation17
ForceFields::UFF::Utils::isInGroup6
ForceFields::UFF::TorsionAngleContrib::TorsionAngleContrib
ForceFields::UFF::TorsionAngleContrib::calcTorsionParams
ForceFields::UFF::TorsionAngleContrib::getEnergy
ForceFields::UFF::TorsionAngleContrib::getGrad
ForceFields::UFF::TorsionAngleContrib::getThetaDeriv
```

## Required Handling

Do not create a Rust `calc_torsion_energy` helper as a substitute for the
missing RDKit symbol.

Do not substitute similarly shaped helpers from another force-field family or
from old COSMolKit code.

Future torsion energy changes must remain anchored to
`TorsionAngleContrib::getEnergy` unless the pinned RDKit source gains a real
UFF `calcTorsionEnergy` symbol.
