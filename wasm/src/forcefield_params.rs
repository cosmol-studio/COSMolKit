//! Native force-field input configuration; executor state is never reconstructed.
use cosmolkit as ck;
use wasm_bindgen::prelude::*;
crate::enums::enum_value!(MmffVariant,"MmffVariant"{Mmff94="MMFF94",Mmff94s="MMFF94s"});
ck::__ck_configuration_ForceFieldDistanceConstraintParams!(
    crate::project_parameters,
    ck::ForceFieldDistanceConstraintParams,
    {},
    { ForceFieldDistanceConstraintParams }
);
ck::__ck_configuration_ForceFieldAngleConstraintParams!(
    crate::project_parameters,
    ck::ForceFieldAngleConstraintParams,
    {},
    { ForceFieldAngleConstraintParams }
);
ck::__ck_configuration_ForceFieldTorsionConstraintParams!(
    crate::project_parameters,
    ck::ForceFieldTorsionConstraintParams,
    {},
    { ForceFieldTorsionConstraintParams }
);
ck::__ck_configuration_ForceFieldPositionConstraintParams!(
    crate::project_parameters,
    ck::ForceFieldPositionConstraintParams,
    {},
    { ForceFieldPositionConstraintParams }
);
ck::__ck_configuration_ForceFieldMinimizeParams!(
    crate::project_parameters,
    ck::ForceFieldMinimizeParams,
    {},
    { ForceFieldMinimizeParams }
);
ck::__ck_configuration_MmffForceFieldParams!(
    crate::project_parameters,
    ck::MmffForceFieldParams,
    {},
    { MmffForceFieldParams }
);
ck::__ck_configuration_UffForceFieldParams!(
    crate::project_parameters,
    ck::UffForceFieldParams,
    {},
    { UffForceFieldParams }
);
ck::__ck_configuration_MmffPropertiesParams!(
    crate::project_parameters,
    ck::MmffPropertiesParams,
    {},
    { MmffPropertiesParams }
);
