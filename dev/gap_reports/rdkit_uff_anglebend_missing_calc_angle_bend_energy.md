# RDKit UFF AngleBend Missing `calcAngleBendEnergy` Audit

## Scope

This audit covers the previously planned port target:

```text
ForceFields::UFF::AngleBendContrib::Utils::calcAngleBendEnergy
```

The audited source files are:

```text
third_party/rdkit/Code/ForceField/UFF/AngleBend.h
third_party/rdkit/Code/ForceField/UFF/AngleBend.cpp
```

## Finding

`ForceFields::UFF::AngleBendContrib::Utils::calcAngleBendEnergy` does not
exist in the audited UFF source files. No Rust helper with this name may be
invented for UFF, and no similarly named helper from another force-field family
may be substituted.

The UFF AngleBend source-backed symbols present in the audited files are:

```text
ForceFields::UFF::AngleBendContrib::AngleBendContrib
ForceFields::UFF::AngleBendContrib::getEnergy
ForceFields::UFF::AngleBendContrib::getGrad
ForceFields::UFF::AngleBendContrib::getEnergyTerm
ForceFields::UFF::AngleBendContrib::getThetaDeriv
ForceFields::UFF::Utils::calcAngleForceConstant
ForceFields::UFF::Utils::calcAngleBendGrad
```

## Source Evidence

`AngleBend.h` declares the private UFF energy-shape helpers:

```text
double getEnergyTerm(double cosTheta, double sinThetaSq) const;
double getThetaDeriv(double cosTheta, double sinTheta) const;
```

`AngleBend.h` declares the UFF `Utils` functions:

```text
RDKIT_FORCEFIELD_EXPORT double calcAngleForceConstant(
RDKIT_FORCEFIELD_EXPORT void calcAngleBendGrad(RDGeom::Point3D *r, double *dist,
```

`AngleBend.cpp` defines the corresponding UFF implementations:

```text
double calcAngleForceConstant(double theta0, double bondOrder12,
void calcAngleBendGrad(RDGeom::Point3D *r, double *dist, double **g,
double AngleBendContrib::getEnergyTerm(double cosTheta,
double AngleBendContrib::getThetaDeriv(double cosTheta, double sinTheta) const {
```

There is no UFF declaration or definition for:

```text
calcAngleBendEnergy
```

## Non-Substitution Rule

`calcAngleBendEnergy` exists under RDKit MMFF AngleBend, not UFF. That MMFF
symbol is outside this UFF audit scope and must not be used as a substitute for
the missing UFF symbol.

The UFF port must continue through source-backed UFF behavior only:

- `calcAngleForceConstant`
- `calcAngleBendGrad`
- `AngleBendContrib::getEnergy`
- `AngleBendContrib::getGrad`
- `AngleBendContrib::getEnergyTerm`
- `AngleBendContrib::getThetaDeriv`

Any future plan item referring to UFF `calcAngleBendEnergy` is a plan error and
must be corrected before implementation.
