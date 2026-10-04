// Fixed source-only property-presence regression preparation.
// RDKit source pin351f8f378f8ad6bbd517980c38896e66bf907af8c; BSD-3-Clause.
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <cstdint>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

namespace {

std::uint64_t bits(double value) {
  std::uint64_t result = 0;
  static_assert(sizeof(result) == sizeof(value));
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

void probe(const std::string &label, const std::string &smiles, int state, bool orient) {
  std::clog << "D2LAST\tCASE\t" << label << '\n';
  RDKit::Chirality::setUseLegacyStereoPerception(true);
  RDKit::SmilesParserParams parser;
  parser.sanitize = true;
  parser.removeHs = true;
  parser.allowCXSMILES = true;
  parser.strictCXSMILES = true;
  parser.parseName = true;
  parser.skipCleanup = false;
  parser.debugParse = 0;

  std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(smiles, parser));
  if (!mol) {
    std::cout << label << "\tERROR\tSmilesToMol returned null\n";
    return;
  }

  mol->clearProp("_StereochemDone");
  switch (state) {
    case 0: break;
    case 1: mol->setProp("_StereochemDone", 0, true); break;
    case 2: mol->setProp("_StereochemDone", 1, true); break;
    case 3: mol->setProp("_StereochemDone", 0, false); break;
    case 4: mol->setProp("_StereochemDone", std::string(""), true); break;
    case 5: mol->setProp("_StereochemDone", std::string("0"), true); break;
    case 6: mol->setProp("_StereochemDone", std::string("false"), true); break;
    case 7: mol->setProp("_stereochemDone", 1); break;
    default: throw std::runtime_error("invalid frozen state");
  }
  RDDepict::preferCoordGen = false;
  RDGeom::INT_POINT2D_MAP coord_map;
  RDDepict::Compute2DCoordParameters params;
  params.canonOrient = orient;
  params.clearConfs = true;
  params.coordMap = &coord_map;
  params.nFlipsPerSample = 0;
  params.nSamples = 0;
  params.sampleSeed = 0;
  params.permuteDeg4Nodes = false;
  params.forceRDKit = false;
  params.useRingTemplates = false;
  std::cout << label << "\tPROFILE\tparser=1,1,1,1,1,0,0\tgeometry=" << orient << ",1,0,0,0,0,0,0\tlegacy=1\tpreferCoordGen=0\n";
  const auto conf_id = RDDepict::compute2DCoords(*mol, params);
  const auto &conf = mol->getConformer(conf_id);

  std::cout << label << "\tMETA\t" << mol->getNumAtoms() << '\t'
            << mol->getNumBonds() << '\t' << conf_id << '\n';
  for (const auto *atom : mol->atoms()) {
    std::cout << label << "\tATOM\t" << atom->getIdx() << '\t'
              << atom->getAtomicNum() << '\t' << atom->getIsotope() << '\t'
              << atom->getFormalCharge() << '\t' << atom->getNumExplicitHs()
              << '\t' << atom->getNoImplicit() << '\t'
              << atom->getNumRadicalElectrons() << '\t' << atom->getIsAromatic()
              << '\t' << static_cast<int>(atom->getHybridization()) << '\t'
              << static_cast<int>(atom->getChiralTag()) << '\t'
              << atom->getAtomMapNum() << '\n';
  }
  for (const auto *bond : mol->bonds()) {
    std::cout << label << "\tBOND\t" << bond->getIdx() << '\t'
              << bond->getBeginAtomIdx() << '\t' << bond->getEndAtomIdx()
              << '\t' << static_cast<int>(bond->getBondType()) << '\t'
              << bond->getIsAromatic() << '\t' << bond->getIsConjugated()
              << '\t' << static_cast<int>(bond->getBondDir()) << '\t'
              << static_cast<int>(bond->getStereo()) << '\t';
    const auto &stereo_atoms = bond->getStereoAtoms();
    for (std::size_t i = 0; i < stereo_atoms.size(); ++i) {
      if (i) {
        std::cout << ',';
      }
      std::cout << stereo_atoms[i];
    }
    std::cout << '\n';
  }
  for (unsigned int atom_idx = 0; atom_idx < mol->getNumAtoms(); ++atom_idx) {
    const auto &point = conf.getAtomPos(atom_idx);
    std::cout << label << "\tXY\t" << atom_idx << '\t' << bits(point.x)
              << '\t' << bits(point.y) << '\n';
  }
}

}  // namespace

int main() {
  std::string line;
  while (std::getline(std::cin, line)) {
    const auto separator = line.find('\t');
    if (separator == std::string::npos) {
      std::cerr << "input row lacks label/SMILES tab separator\n";
      return 2;
    }
    try {
      const auto second = line.find('\t', separator + 1);
      const auto third = line.find('\t', second + 1);
      if (second == std::string::npos || third == std::string::npos) return 2;
      const int state = std::stoi(line.substr(separator + 1, second - separator - 1));
      const int orient = std::stoi(line.substr(second + 1, third - second - 1));
      if (orient != 0 && orient != 1) return 2;
      probe(line.substr(0, separator), line.substr(third + 1), state, orient != 0);
    } catch (const std::exception &error) {
      std::cout << line.substr(0, separator) << "\tERROR\t" << error.what()
                << '\n';
    }
  }
  return 0;
}
