#include <GraphMol/MolOps.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Canon.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {

using namespace RDKit;

void print_bond_state(const ROMol &mol, unsigned int bondIdx) {
  const auto *bond = mol.getBondWithIdx(bondIdx);
  const auto *ri = mol.getRingInfo();
  std::cout << "bond " << bondIdx << " "
            << bond->getBeginAtomIdx() << "-" << bond->getEndAtomIdx()
            << " dir=" << static_cast<int>(bond->getBondDir())
            << " stereo=" << static_cast<int>(bond->getStereo());
  const auto &stereoAtoms = bond->getStereoAtoms();
  if (stereoAtoms.size() == 2) {
    std::cout << " stereo_atoms=" << stereoAtoms[0] << "," << stereoAtoms[1];
  } else {
    std::cout << " stereo_atoms=-";
  }
  std::cout << " ring_count=" << ri->numBondRings(bondIdx)
            << " min_ring=" << ri->minBondRingSize(bondIdx) << "\n";
}

void print_checkpoint(const std::string &name, const ROMol &mol,
                      const std::vector<unsigned int> &focusBonds) {
  const auto *ri = mol.getRingInfo();
  std::cout << "checkpoint=" << name
            << " rings_initialized=" << ri->isInitialized()
            << " is_symm_sssr=" << ri->isSymmSssr()
            << " stereo_done=" << mol.hasProp(common_properties::_StereochemDone)
            << "\n";
  for (auto bondIdx : focusBonds) {
    print_bond_state(mol, bondIdx);
  }
}

std::vector<unsigned int> rank_atoms_for_canonical_smiles(ROMol &mol) {
  std::vector<unsigned int> ranks(mol.getNumAtoms());
  const bool breakTies = true;
  const bool includeChiralPresence = false;
  const bool includeIsotopes = false;
  const bool includeChirality = false;
  const bool includeStereoGroups = false;
  const bool useNonStereoRanks = false;
  const bool includeAtomMaps = true;
  Canon::rankMolAtoms(mol, ranks, breakTies, includeChirality, includeIsotopes,
                      includeAtomMaps, includeChiralPresence,
                      includeStereoGroups, useNonStereoRanks);
  return ranks;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "usage: rdkit_writer_probe <smiles>\n";
    return 2;
  }

  const std::string smiles = argv[1];
  bool sanitize = false;
  SmilesParserParams parseParams;
  parseParams.sanitize = sanitize;
  std::unique_ptr<ROMol> mol{SmilesToMol(smiles, parseParams)};
  if (!mol) {
    std::cerr << "failed to parse smiles\n";
    return 1;
  }

  const std::vector<unsigned int> focusBonds = {1u, 3u, 4u};
  print_checkpoint("raw_parse", *mol, focusBonds);

  for (auto atom : mol->atoms()) {
    atom->updatePropertyCache(false);
  }
  print_checkpoint("post_update_property_cache", *mol, focusBonds);

  MolOps::symmetrizeSSSR(*mol);
  print_checkpoint("post_symmetrize_sssr", *mol, focusBonds);

  MolOps::findPotentialStereoBonds(*mol, false);
  print_checkpoint("post_find_potential_stereo_bonds", *mol, focusBonds);

  MolOps::assignStereochemistry(*mol, false);
  print_checkpoint("post_assign_stereochemistry", *mol, focusBonds);

  std::vector<Canon::AtomColors> colors(mol->getNumAtoms(), Canon::WHITE_NODE);
  auto ranks = rank_atoms_for_canonical_smiles(*mol);
  int startAtomIdx = -1;
  unsigned int nextRank = mol->getNumAtoms() + 1;
  for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
    if (colors[i] == Canon::WHITE_NODE && ranks[i] < nextRank) {
      nextRank = ranks[i];
      startAtomIdx = static_cast<int>(i);
    }
  }
  Canon::MolStack stack;
  Canon::canonicalizeFragment(*mol, startAtomIdx, colors, ranks, stack, nullptr,
                              nullptr, false, false, true);
  print_checkpoint("post_canonicalize_fragment", *mol, focusBonds);

  SmilesWriteParams writeParams;
  writeParams.canonical = true;
  writeParams.cleanStereo = false;
  writeParams.doIsomericSmiles = true;
  const auto out = MolToSmiles(*mol, writeParams);
  std::cout << "checkpoint=final_output smiles=" << out << "\n";
  return 0;
}
