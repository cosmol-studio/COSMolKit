import {ResidueCode,ResidueInfoKind,ResidueInfo,ResidueIdentity,ResidueCodeParseError,ResidueSequenceError,residueInfo,residueInfoChecked,findResidueInfo,findResidueInfoIndex,residueCode,expandOneLetter,expandOneLetterSequence} from 'cosmolkit-generated';
const a:ResidueInfo=findResidueInfo('ALA'),code:ResidueCode=a.code(),parent:ResidueCode|null=a.parentStandardCode(),kind:ResidueInfoKind=ResidueInfoKind.AA,name:string=kind.name(),num:number=kind.value,maybe:ResidueInfo|null=residueInfoChecked(368),identity=ResidueIdentity.new('MSE'),letter:string|null=expandOneLetter('A',kind),seq:string[]=expandOneLetterSequence('AG',kind);const index:number=findResidueInfoIndex('MSE'),fixed:ResidueInfo=residueInfo(index),value:ResidueCode=residueCode(identity.name());declare const error:ResidueCodeParseError;const input:string=error.input();declare const sequenceError:ResidueSequenceError;const errorKind:string=sequenceError.kind;
// @ts-expect-error enum kind class is required
expandOneLetter('A',1);
// @ts-expect-error nullable result
const bad:ResidueInfo=residueInfoChecked(368);
// @ts-expect-error readonly kind value
kind.value=3;
void [code,parent,name,num,maybe,letter,seq,fixed,value,input,errorKind,bad];
