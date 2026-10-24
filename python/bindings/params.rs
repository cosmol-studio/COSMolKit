use cosmolkit as ck;

ck::project_configuration!(ck::MmffPropertiesParams);

ck::project_configuration!(ck::SmartsWriteParams);

ck::project_configuration!(ck::ForceFieldDistanceConstraintParams);

ck::project_configuration!(ck::ForceFieldAngleConstraintParams);

ck::project_configuration!(ck::ForceFieldTorsionConstraintParams);

ck::project_configuration!(ck::MaccsFingerprintParams);

ck::project_configuration!(ck::McsAtomCompareParameters);
ck::project_configuration!(ck::McsBondCompareParameters);
ck::project_configuration!(ck::MmffEvaluationParams);
ck::project_configuration!(ck::MmffOptimizationParams);
ck::project_configuration!(ck::MmffConformerOptimizationParams);
ck::project_configuration!(ck::Coordinate3DInputParams);
ck::project_configuration!(ck::Replace3DCoordinatesParams);
ck::project_configuration!(ck::MorganCallParams);
ck::project_configuration!(ck::TopologicalTorsionCallParams);

ck::project_configuration!(ck::RemoveHsParams);
ck::project_configuration!(ck::KekulizeParams);
ck::project_configuration!(ck::RingSearchParams);
ck::project_configuration!(ck::StructureTagParams);
ck::project_configuration!(ck::PotentialStereoParams);
ck::project_configuration!(ck::DistanceMatrixParams);
ck::project_configuration!(ck::DistanceMatrix3dParams);
ck::project_configuration!(ck::AtomPositionParams);
ck::project_configuration!(ck::RandomSmilesWriteParams);
ck::project_configuration!(ck::PatternFingerprintParams);
ck::project_configuration!(ck::TopologicalFingerprintOutputRequest);
ck::project_configuration!(ck::UffEvaluationParams);
ck::project_configuration!(ck::UffOptimizationParams);
ck::project_configuration!(ck::UffConformerOptimizationParams);
ck::project_configuration!(ck::XyzWriteParams);
ck::project_configuration!(ck::InchiReadParams);
ck::project_configuration!(ck::InchiWriteParams);
ck::project_configuration!(ck::ReactionParseParams);
ck::project_configuration!(ck::ReactionValidationParams);
ck::project_configuration!(ck::ReactionApplyParams);
ck::project_configuration!(ck::ReactionTemplateRemovalParams);
ck::project_configuration!(ck::BioPdbReadParams);
ck::project_configuration!(ck::BioPdbWriteParams);
ck::project_configuration!(ck::BioMoleculeParams);

/// Required inputs stay required; Rust's canonical constructor owns the value.
ck::project_configuration!(ck::ForceFieldPositionConstraintParams);
