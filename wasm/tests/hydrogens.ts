import {Molecule,AddHsParams,RemoveHsParams,HydrogenError} from 'cosmolkit-generated';const m=Molecule.fromSmiles('CCO'),a=new AddHsParams(),r=new RemoveHsParams();const v1:Molecule=m.withHydrogens(),v2:Molecule=m.withHydrogensWithParams(a),v3:Molecule=m.withoutHydrogens(),v4:Molecule=m.withoutHydrogensWithParams(r),i1:void=m.addHydrogens(),i2:void=m.addHydrogensWithParams(a),i3:void=m.removeHydrogens(),i4:void=m.removeHydrogensWithParams(r);declare const e:HydrogenError;const domain:string=e.domain,cause:Error|null=e.cause;
// @ts-expect-error typed addition parameters
m.withHydrogensWithParams(r);
// @ts-expect-error in-place returns void
const output:Molecule=m.addHydrogens();
