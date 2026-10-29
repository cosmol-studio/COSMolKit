// Use the full package, or change the import to @cosmol-studio/cosmolkit-core-3d.
import init, * as ck from "@cosmol-studio/cosmolkit";

await init();

const input = ck.molFromSmiles("CCO");
const withHydrogens = input.withHydrogens();
const molecule = withHydrogens.with3dConformer({ randomSeed: 42 });
withHydrogens.free();
input.free();

const original = molecule.coordinates3d(0).flat();
try {
    for (const model of ["UFF", ck.MmffVariant.Mmff94, ck.MmffVariant.Mmff94s]) {
        const field = model === "UFF"
            ? molecule.uffForceField()
            : molecule.mmffForceField({ mmffVariant: model });
        // These values are explicit example choices, not API defaults.
        const distance = new ck.ForceFieldDistanceConstraintParams({relative: false, minDistance: 1.4, maxDistance: 1.6, forceConstant: 50});
        let outcome;
        let initial;
        let relaxed;
        try {
            // The same method accepts either a parameter instance or a plain options object.
            distance.forceConstant = 100;
            field.addDistanceConstraint_(0, 1, distance);
            field.addAngleConstraint_(0, 1, 2, {
                relative: true, minAngleDeg: -5, maxAngleDeg: 5, forceConstant: 1,
            });
            field.addTorsionConstraint_(3, 0, 1, 2, {
                relative: true, minDihedralDeg: -10, maxDihedralDeg: 10, forceConstant: 1,
            });
            // Distances/displacements are in Angstroms; angles/dihedrals are in degrees.
            // Relative limits and the soft O position anchor are captured at addition.
            field.addPositionConstraint_(2, { maxDisplacement: 0.25, forceConstant: 10 });
            field.setFixedAtoms_([0]); // A fixed atom cannot move during minimization.
            const pinned = field.position(0);
            const oxygen = field.position(2);
            field.setPosition_(2, [oxygen[0] + 0.3, oxygen[1] + 0.1, oxygen[2] + 0.2]);
            initial = field.energyGradient();
            outcome = field.minimize_();
            if (!Number.isFinite(outcome.energy) || outcome.energy > initial.energy) {
                throw new Error(`${model}: minimization did not lower the finite energy`);
            }
            if (!field.position(0).every((value, i) => Object.is(value, pinned[i]))) {
                throw new Error(`${model}: fixed atom moved`);
            }
            // Force-field edits never write back. Apply the snapshot explicitly to a new value.
            const positions = field.positions;
            const flattened = positions.flat();
            relaxed = molecule.with3dCoordinates(positions);
            if (!molecule.coordinates3d(0).flat().every((value, i) => Object.is(value, original[i]))) {
                throw new Error(`${model}: source molecule changed`);
            }
            if (!relaxed.coordinates3d(0).flat().every((value, i) => Object.is(value, flattened[i]))) {
                throw new Error(`${model}: result coordinates were not installed`);
            }
            console.log(
                `${model} restrained: energy=${initial.energy.toFixed(6)} -> ${outcome.energy.toFixed(6)}, `
                + `converged=${outcome.converged}, iterations=${outcome.iterations}`,
            );
        } finally {
            relaxed?.free();
            outcome?.free();
            initial?.free();
            distance.free();
            field.free();
        }
    }
} finally {
    molecule.free();
}
