import {Molecule,QueryGraph,Fingerprint,LayeredFingerprintParams,LayeredFingerprintResult,PatternFingerprintParams,fingerprintLayeredQueryWithParams,fingerprintLayeredQueryWithOutputWithParams,fingerprintPatternQuery,fingerprintPatternQueryWithParams} from 'cosmolkit-generated';
const m=Molecule.fromSmiles('CCO'),q=QueryGraph.fromSmarts('[#8]'),p=new LayeredFingerprintParams(),pattern=new PatternFingerprintParams();
const a:Fingerprint=m.fingerprintLayered(),b:Fingerprint=m.fingerprintLayeredWithParams(p),c:LayeredFingerprintResult=m.fingerprintLayeredWithOutput(),d:LayeredFingerprintResult=m.fingerprintLayeredWithOutputWithParams(p),e:Fingerprint=m.fingerprintPattern(),f:Fingerprint=m.fingerprintPatternWithParams(pattern);
const g:Fingerprint=fingerprintLayeredQueryWithParams(q,p),h:LayeredFingerprintResult=fingerprintLayeredQueryWithOutputWithParams(q,p),i:Fingerprint=fingerprintPatternQuery(q),j:Fingerprint=fingerprintPatternQueryWithParams(q,pattern);
const counts:number[]|null=h.atomCounts();
const configured:Fingerprint=m.fingerprintLayered({fpSize:128}),withParams:LayeredFingerprintResult=m.fingerprintLayeredWithOutput(p),withOptions:LayeredFingerprintResult=m.fingerprintLayeredWithOutput({fpSize:128,atomCounts:[0,0,0]});
// @ts-expect-error options reject unknown fields
m.fingerprintLayeredWithOutput({unknown:true});
// @ts-expect-error configured fingerprint size remains a number
m.fingerprintLayered({fpSize:'128'});
// @ts-expect-error queries remain typed query values
fingerprintPatternQuery(m);
// @ts-expect-error required explicit query configuration
fingerprintLayeredQueryWithParams(q);
