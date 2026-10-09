import {Molecule, QueryGraph, McsResult, McsParameters, McsAtomComparator, maximumCommonSubstructure} from 'cosmolkit-generated';
const inputs = [Molecule.fromSmiles('CCO'), Molecule.fromSmiles('CCN')];
const p = new McsParameters();
p.timeout = 10;
p.atomComparator = McsAtomComparator.Any;
p.atomCompareParameters.matchFormalCharge = true;
const result: McsResult = maximumCommonSubstructure(inputs, p);
const other: McsResult = maximumCommonSubstructure(inputs, {timeout: 5, atomCompareParameters: {matchFormalCharge: true}});
const query: QueryGraph | null = result.query;
const alternatives: Map<string, QueryGraph> = other.degenerate;
void query; void alternatives;
// @ts-expect-error unknown configuration fields must not be accepted
maximumCommonSubstructure(inputs, {unknownOption: true});
