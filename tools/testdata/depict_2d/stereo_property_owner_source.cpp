// Source-only observation driver; RDKit pin351f8f378f8ad6bbd517980c38896e66bf907af8c.
// No chemistry implementation: the executable observation seam forwards each
// native assignAtomChiralCodes invocation exactly once to the unchanged library.
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <dlfcn.h>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
std::string active;
unsigned calls = 0;

std::uint64_t bits(double value) {
  std::uint64_t result;
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

void properties(const std::string &phase, const std::string &owner,
                const RDKit::RDProps &props) {
  std::vector<std::string> computed;
  props.getPropIfPresent("__computedProps", computed);
  for (const auto &entry : props.getDict()) {
    std::cout << active << '\t' << phase << "\tPROP\t" << owner << '\t'
              << std::quoted(entry.key) << '\t' << entry.val.getTag() << '\t'
              << (std::find(computed.begin(), computed.end(), entry.key) != computed.end()) << '\t';
    switch (entry.val.getTag()) {
      case RDKit::RDTypeTag::IntTag: std::cout << RDKit::rdvalue_cast<int>(entry.val); break;
      case RDKit::RDTypeTag::UnsignedIntTag: std::cout << RDKit::rdvalue_cast<unsigned>(entry.val); break;
      case RDKit::RDTypeTag::DoubleTag: std::cout << bits(RDKit::rdvalue_cast<double>(entry.val)); break;
      case RDKit::RDTypeTag::BoolTag: std::cout << RDKit::rdvalue_cast<bool>(entry.val); break;
      case RDKit::RDTypeTag::StringTag: std::cout << std::quoted(RDKit::rdvalue_cast<std::string>(entry.val)); break;
      case RDKit::RDTypeTag::VecStringTag:
        for (const auto &value : RDKit::rdvalue_cast<std::vector<std::string>>(entry.val)) std::cout << std::quoted(value) << ',';
        break;
      case RDKit::RDTypeTag::VecIntTag:
        for (auto value : RDKit::rdvalue_cast<std::vector<int>>(entry.val)) std::cout << value << ',';
        break;
      case RDKit::RDTypeTag::VecUnsignedIntTag:
        for (auto value : RDKit::rdvalue_cast<std::vector<unsigned>>(entry.val)) std::cout << value << ',';
        break;
      default: throw std::runtime_error("unhandled native property tag at " + owner + ":" + entry.key);
    }
    std::cout << '\n';
  }
}

void state(const std::string &phase, const RDKit::ROMol &mol) {
  std::cout << active << '\t' << phase << "\tSIZE\t" << mol.getNumAtoms() << '\t' << mol.getNumBonds() << '\n';
  properties(phase, "MOL", mol);
  for (const auto *atom : mol.atoms()) {
    std::cout << active << '\t' << phase << "\tATOM\t" << atom->getIdx() << '\t'
              << atom->getAtomicNum() << '\t' << atom->getIsotope() << '\t'
              << atom->getFormalCharge() << '\t' << atom->getNumExplicitHs() << '\t'
              << atom->getNoImplicit() << '\t' << atom->getNumRadicalElectrons() << '\t'
              << atom->getIsAromatic() << '\t' << static_cast<int>(atom->getHybridization()) << '\t'
              << static_cast<int>(atom->getChiralTag()) << '\t' << atom->getAtomMapNum() << '\n';
    properties(phase, "A" + std::to_string(atom->getIdx()), *atom);
  }
  for (const auto *bond : mol.bonds()) {
    std::cout << active << '\t' << phase << "\tBOND\t" << bond->getIdx() << '\t'
              << bond->getBeginAtomIdx() << '\t' << bond->getEndAtomIdx() << '\t'
              << static_cast<int>(bond->getBondType()) << '\t' << bond->getIsAromatic() << '\t'
              << bond->getIsConjugated() << '\t' << static_cast<int>(bond->getBondDir()) << '\t'
              << static_cast<int>(bond->getStereo()) << '\t';
    for (auto index : bond->getStereoAtoms()) std::cout << index << ',';
    std::cout << '\n';
    properties(phase, "B" + std::to_string(bond->getIdx()), *bond);
  }
}
}  // namespace

namespace RDKit::Chirality {
void legacyStereoPerception(ROMol &, bool, bool);
std::pair<bool, bool> assignAtomChiralCodes(ROMol &mol, std::vector<unsigned> &ranks, bool possible) {
  using Owner = std::pair<bool, bool> (*)(ROMol &, std::vector<unsigned> &, bool);
  static const auto owner = reinterpret_cast<Owner>(dlsym(RTLD_NEXT,
      "_ZN5RDKit9Chirality21assignAtomChiralCodesERNS_5ROMolERSt6vectorIjSaIjEEb"));
  if (!owner) throw std::runtime_error("actual native assignAtomChiralCodes binding unavailable");
  const bool observe = !active.empty();
  const unsigned ordinal = observe ? calls++ : 0;
  if (observe) state("INNER_PRE_" + std::to_string(ordinal), mol);
  const auto result = owner(mol, ranks, possible);
  if (observe) {
    state("INNER_POST_" + std::to_string(ordinal), mol);
    std::cout << active << "\tRETURN\t" << ordinal << '\t' << possible << '\t'
              << result.first << '\t' << result.second << '\t';
    for (auto rank : ranks) std::cout << rank << ',';
    std::cout << '\n';
  }
  return result;
}
}  // namespace RDKit::Chirality

int main() {
  RDKit::Chirality::setUseLegacyStereoPerception(true);
  std::string line;
  unsigned outer = 0;
  while (std::getline(std::cin, line)) {
    auto tab = line.find('\t');
    if (tab == std::string::npos) return 2;
    const auto label = line.substr(0, tab);
    const auto smiles = line.substr(tab + 1);
    for (int family = 0; family < 2; ++family) {
      for (int first = 0; first < 2; ++first) for (int second = 0; second < 2; ++second) {
        active.clear();
        RDKit::SmilesParserParams params;
        params.sanitize = family == 0 ? true : second;
        params.removeHs = family == 0 ? true : first;
        params.allowCXSMILES = true; params.strictCXSMILES = true;
        params.parseName = true; params.skipCleanup = false; params.debugParse = 0;
        std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(smiles, params));
        if (!mol) throw std::runtime_error("source parser returned null");
        active = label + (family == 0 ? "_OWNER_" : "_LIFECYCLE_") + std::to_string(first) + std::to_string(second);
        calls = 0;
        if (family == 0) {
          state("OUTER_PRE", *mol);
          RDKit::Chirality::legacyStereoPerception(*mol, first, second);
          state("OUTER_POST", *mol);
        } else state("CONSTRUCTOR", *mol);
        std::cout << active << "\tEND\t" << calls << '\n';
        active.clear(); ++outer;
      }
    }
  }
  std::cout << "CENSUS\t" << outer << '\n';
  return outer == 24 ? 0 : 2;
}
