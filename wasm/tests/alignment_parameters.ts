import { AlignmentAtomMap, AlignmentParameters, BestAlignmentParameters, CoordinateRmsdParameters, AllConformerRmsdParameters, ConformerAlignmentParameters } from "cosmolkit-generated";

const map: AlignmentAtomMap = AlignmentAtomMap.new(1, 2);
const params: AlignmentParameters = new AlignmentParameters();
const mapCopy: AlignmentAtomMap[] | null = params.atomMap;
const weights: number[] | null = params.weights;
params.atomMap = [map];
params.weights = new Float64Array([1]);
params.probeConformerId = -1;
const best: BestAlignmentParameters = BestAlignmentParameters.new(-1, -1, [[map]]);
const coordinate: CoordinateRmsdParameters = new CoordinateRmsdParameters();
const all: AllConformerRmsdParameters = AllConformerRmsdParameters.new();
const conformers: ConformerAlignmentParameters = new ConformerAlignmentParameters([1], new Uint32Array([0]));
const maps: AlignmentAtomMap[][] = best.atomMaps;
const indices: number[] | null = conformers.atomIndices;
conformers.conformerIds = null;
void [mapCopy, weights, coordinate, all, maps, indices];
