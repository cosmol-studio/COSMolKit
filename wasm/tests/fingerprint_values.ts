import {Fingerprint,SparseBitFingerprint,SparseCountFingerprint,SparseCountFingerprint32,FingerprintError} from 'cosmolkit-generated';
const fp=Fingerprint.fromOnBits(8,new Uint32Array([1,7])),bits:number[]=fp.onBits(),n:number=fp.nBits(),similarity:number=fp.tanimoto(fp);
const wide=SparseCountFingerprint.new(18446744073709551615n),small=SparseCountFingerprint32.new(8);
wide.setValue(9007199254740993n,-2);small.setValue(1,2);
const length:bigint=wide.length(),shortLength:number=small.length(),map:Map<bigint,number>=wide.nonzeroElements(),shortMap:Map<number,number>=small.nonzeroElements();
for(const value of [wide.fuzzyAnd(wide),wide.fuzzyOr(wide),wide.withAdded(wide),wide.withSubtracted(wide),wide.withAddedScalar(1),wide.withSubtractedScalar(1),wide.withMultipliedScalar(2),wide.withDividedScalar(2)]){const total:number=value.totalValue();const at:number=value.value(1n);void [total,at];}
declare const sparse:SparseBitFingerprint;const signed:number[]=sparse.onBits();
declare const error:FingerprintError;const cause:null=error.cause,domain:string=error.domain;
// @ts-expect-error u64 requires exact bigint
SparseCountFingerprint.new(8);
// @ts-expect-error u32 requires number
SparseCountFingerprint32.new(8n);
// @ts-expect-error boolean option
wide.totalValue('false');
void [bits,n,similarity,length,shortLength,map,shortMap,signed,cause,domain];
