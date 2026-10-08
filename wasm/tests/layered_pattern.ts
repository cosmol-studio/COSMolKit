import {Molecule,QueryGraph,Fingerprint,LayeredFingerprintParams,LayeredFingerprintResult,PatternFingerprintParams,layeredQueryFingerprintWithParams,layeredQueryFingerprintWithOutputWithParams,patternQueryFingerprint,patternQueryFingerprintWithParams} from 'cosmolkit-generated';
const m=Molecule.fromSmiles('CCO'),q=QueryGraph.fromSmarts('[#8]'),p=new LayeredFingerprintParams(),pattern=new PatternFingerprintParams();
const a:Fingerprint=m.layeredFingerprint(),b:Fingerprint=m.layeredFingerprintWithParams(p),c:LayeredFingerprintResult=m.layeredFingerprintWithOutput(),d:LayeredFingerprintResult=m.layeredFingerprintWithOutputWithParams(p),e:Fingerprint=m.patternFingerprint(),f:Fingerprint=m.patternFingerprintWithParams(pattern);
const g:Fingerprint=layeredQueryFingerprintWithParams(q,p),h:LayeredFingerprintResult=layeredQueryFingerprintWithOutputWithParams(q,p),i:Fingerprint=patternQueryFingerprint(q),j:Fingerprint=patternQueryFingerprintWithParams(q,pattern);
const counts:number[]|null=h.atomCounts();
// @ts-expect-error queries remain typed query values
patternQueryFingerprint(m);
// @ts-expect-error required explicit query configuration
layeredQueryFingerprintWithParams(q);
