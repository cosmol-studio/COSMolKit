import { Element, ElementInfo, elementInfo } from "cosmolkit-generated";

const fromSymbol: (symbol: string) => Element | null = Element.fromSymbol;
const fromAtomicNumber: (number: number) => Element | null = Element.fromAtomicNumber;
declare const element: Element;
const atomicNumber: number = element.atomicNumber();
const symbol: string = element.symbol();
const info: ElementInfo = elementInfo(element);
const copiedElement: Element = info.element();
const period: number = info.period();
const outerElectrons: number = info.outerElectrons();
const valences: Int32Array = info.valences();
const radius: number = info.rb0();
const weight: number = info.atomicWeight();
void [fromSymbol, fromAtomicNumber, atomicNumber, symbol, copiedElement, period, outerElectrons, valences, radius, weight];
