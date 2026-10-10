import {Molecule, Atom, Bond, AtomMetadata, ChiralTag, Hybridization, BondDirection, BondStereo, BondOrder, TemplateAttachment, TemplateAttachmentOrder} from 'cosmolkit-generated';
const molecule = Molecule.fromSmiles('CCO');
const atoms: Atom[] = molecule.atoms();
const bonds: Bond[] = molecule.bonds();
const atom: Atom | null = molecule.atom(0);
const bond: Bond | null = molecule.bond(0);
const metadata: AtomMetadata[] = molecule.atomMetadata(false);
const chirality: ChiralTag = atoms[0].chiralTag();
const hybridization: Hybridization = atoms[0].hybridization();
const direction: BondDirection = bonds[0].direction();
const stereo: BondStereo = bonds[0].stereo();
const order: BondOrder = bonds[0].order();
const attachments = new TemplateAttachmentOrder([new TemplateAttachment(0, 'A')]);
const rows: TemplateAttachment[] = attachments.entries();
// @ts-expect-error atom selection is a checked number
molecule.atom('0');
// @ts-expect-error a bond has no mutable endpoint API
bonds[0].setBegin(2);
void [atom, bond, metadata, chirality, hybridization, direction, stereo, order, rows];
