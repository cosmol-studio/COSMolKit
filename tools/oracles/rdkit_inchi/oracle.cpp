#include <boost/tuple/tuple.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <exception>
#include <functional>
#include <iomanip>
#include <iostream>
#include <list>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <inchi_api.h>

struct ScriptedInchiOutput {
  int return_code = inchi_Ret_OKAY;
  std::vector<inchi_Atom> atoms;
  std::vector<inchi_Stereo0D> stereo0d;
  bool has_message = false;
  bool has_log = false;
  std::string message;
  std::string log;
  bool throw_on_get = false;
  bool throw_on_free = false;
};

static ScriptedInchiOutput scripted_inchi_output;
static std::string scripted_seen_inchi;
static std::string scripted_seen_options;
static unsigned int scripted_get_count = 0;
static unsigned int scripted_free_count = 0;
static unsigned int scripted_outstanding_outputs = 0;

struct ScriptedMolToInchiOutput {
  int return_code = inchi_Ret_OKAY;
  bool has_inchi = false;
  bool has_message = false;
  bool has_log = false;
  bool has_aux_info = false;
  std::string inchi;
  std::string message;
  std::string log;
  std::string aux_info;
  bool throw_on_get = false;
  bool throw_on_free = false;
};

static ScriptedMolToInchiOutput scripted_mol_to_inchi_output;
static std::vector<inchi_Atom> scripted_seen_atoms;
static std::vector<inchi_Stereo0D> scripted_seen_stereo0d;
static bool scripted_seen_null_options = true;
static std::string scripted_seen_generation_options;
static unsigned int scripted_generation_get_count = 0;
static unsigned int scripted_generation_free_count = 0;
static unsigned int scripted_generation_outstanding_outputs = 0;
static std::vector<std::string> scripted_generation_calls;

static int scripted_inchi_key_status = INCHIKEY_OK;
static std::vector<unsigned char> scripted_inchi_key_buffer;
static bool scripted_inchi_key_throw = false;
static std::string scripted_seen_inchi_key_input;
static int scripted_seen_inchi_key_xtra1 = -1;
static int scripted_seen_inchi_key_xtra2 = -1;
static unsigned int scripted_inchi_key_call_count = 0;
static std::vector<std::string> scripted_adapter_calls;

static char *copy_c_string(const std::string &value) {
  auto *result = new char[value.size() + 1];
  std::copy(value.begin(), value.end(), result);
  result[value.size()] = '\0';
  return result;
}

extern "C" int INCHI_DECL GetStructFromINCHI(inchi_InputINCHI *input,
                                               inchi_OutputStruct *output) {
  ++scripted_get_count;
  scripted_seen_inchi = input->szInChI == nullptr ? "" : input->szInChI;
  scripted_seen_options = input->szOptions == nullptr ? "" : input->szOptions;
  if (scripted_inchi_output.throw_on_get) {
    throw std::runtime_error("GetStructFromINCHI");
  }
  std::memset(output, 0, sizeof(*output));
  output->num_atoms = static_cast<AT_NUM>(scripted_inchi_output.atoms.size());
  output->num_stereo0D =
      static_cast<AT_NUM>(scripted_inchi_output.stereo0d.size());
  if (!scripted_inchi_output.atoms.empty()) {
    output->atom = new inchi_Atom[scripted_inchi_output.atoms.size()];
    std::copy(scripted_inchi_output.atoms.begin(),
              scripted_inchi_output.atoms.end(), output->atom);
  }
  if (!scripted_inchi_output.stereo0d.empty()) {
    output->stereo0D =
        new inchi_Stereo0D[scripted_inchi_output.stereo0d.size()];
    std::copy(scripted_inchi_output.stereo0d.begin(),
              scripted_inchi_output.stereo0d.end(), output->stereo0D);
  }
  if (scripted_inchi_output.has_message) {
    output->szMessage = copy_c_string(scripted_inchi_output.message);
  }
  if (scripted_inchi_output.has_log) {
    output->szLog = copy_c_string(scripted_inchi_output.log);
  }
  ++scripted_outstanding_outputs;
  return scripted_inchi_output.return_code;
}

extern "C" void INCHI_DECL FreeStructFromINCHI(inchi_OutputStruct *output) {
  ++scripted_free_count;
  delete[] output->atom;
  delete[] output->stereo0D;
  delete[] output->szMessage;
  delete[] output->szLog;
  output->atom = nullptr;
  output->stereo0D = nullptr;
  output->szMessage = nullptr;
  output->szLog = nullptr;
  if (scripted_outstanding_outputs != 0) {
    --scripted_outstanding_outputs;
  }
  if (scripted_inchi_output.throw_on_free) {
    throw std::runtime_error("FreeStructFromINCHI");
  }
}

extern "C" int INCHI_DECL GetINCHI(inchi_Input *input,
                                     inchi_Output *output) {
  scripted_adapter_calls.emplace_back("GetINCHI");
  scripted_generation_calls.emplace_back("GetINCHI");
  ++scripted_generation_get_count;
  scripted_seen_atoms.clear();
  if (input->num_atoms != 0) {
    scripted_seen_atoms.assign(input->atom,
                               input->atom + input->num_atoms);
  }
  scripted_seen_stereo0d.clear();
  if (input->num_stereo0D != 0) {
    scripted_seen_stereo0d.assign(input->stereo0D,
                                  input->stereo0D + input->num_stereo0D);
  }
  scripted_seen_null_options = input->szOptions == nullptr;
  scripted_seen_generation_options =
      input->szOptions == nullptr ? "" : input->szOptions;
  if (scripted_mol_to_inchi_output.throw_on_get) {
    throw std::runtime_error("GetINCHI");
  }
  std::memset(output, 0, sizeof(*output));
  if (scripted_mol_to_inchi_output.has_inchi) {
    output->szInChI = copy_c_string(scripted_mol_to_inchi_output.inchi);
  }
  if (scripted_mol_to_inchi_output.has_message) {
    output->szMessage = copy_c_string(scripted_mol_to_inchi_output.message);
  }
  if (scripted_mol_to_inchi_output.has_log) {
    output->szLog = copy_c_string(scripted_mol_to_inchi_output.log);
  }
  if (scripted_mol_to_inchi_output.has_aux_info) {
    output->szAuxInfo = copy_c_string(scripted_mol_to_inchi_output.aux_info);
  }
  ++scripted_generation_outstanding_outputs;
  return scripted_mol_to_inchi_output.return_code;
}

extern "C" void INCHI_DECL FreeINCHI(inchi_Output *output) {
  scripted_adapter_calls.emplace_back("FreeINCHI");
  scripted_generation_calls.emplace_back("FreeINCHI");
  ++scripted_generation_free_count;
  delete[] output->szInChI;
  delete[] output->szMessage;
  delete[] output->szLog;
  delete[] output->szAuxInfo;
  output->szInChI = nullptr;
  output->szMessage = nullptr;
  output->szLog = nullptr;
  output->szAuxInfo = nullptr;
  if (scripted_generation_outstanding_outputs != 0) {
    --scripted_generation_outstanding_outputs;
  }
  if (scripted_mol_to_inchi_output.throw_on_free) {
    throw std::runtime_error("FreeINCHI");
  }
}

extern "C" int INCHI_DECL GetINCHIKeyFromINCHI(
    const char *inchi, const int xtra1, const int xtra2, char *inchi_key,
    char *extra1, char *extra2) {
  scripted_adapter_calls.emplace_back("GetINCHIKeyFromINCHI");
  ++scripted_inchi_key_call_count;
  scripted_seen_inchi_key_input = inchi == nullptr ? "" : inchi;
  scripted_seen_inchi_key_xtra1 = xtra1;
  scripted_seen_inchi_key_xtra2 = xtra2;
  if (scripted_inchi_key_throw) {
    throw std::runtime_error("GetINCHIKeyFromINCHI");
  }
  std::memset(inchi_key, 0xA5, 29);
  std::memset(extra1, 0xB6, 65);
  std::memset(extra2, 0xC7, 65);
  if (scripted_inchi_key_buffer.size() > 29) {
    throw std::runtime_error("scripted InChIKey exceeds source buffer");
  }
  std::copy(scripted_inchi_key_buffer.begin(),
            scripted_inchi_key_buffer.end(), inchi_key);
  return scripted_inchi_key_status;
}

namespace RDKit {

class RWMol;

class RangeError final : public std::exception {
 public:
  RangeError(unsigned int index, unsigned int upper_bound)
      : index(index), upper_bound(upper_bound) {}

  const char *what() const noexcept override { return "Range Error"; }

  unsigned int index;
  unsigned int upper_bound;
};

class ValueErrorException final : public std::exception {
 public:
  ValueErrorException(unsigned int bond_index, int bond_type)
      : bond_index(bond_index), bond_type(bond_type) {}

  const char *what() const noexcept override { return "Bad bond type"; }

  unsigned int bond_index;
  int bond_type;
};

class InvariantException final : public std::exception {
 public:
  InvariantException(const char *kind, const char *expression,
                     std::string message)
      : kind(kind), expression(expression), message(std::move(message)) {}

  const char *what() const noexcept override { return message.c_str(); }

  std::string kind;
  std::string expression;
  std::string message;
};

#define PRECONDITION(expr, message)                                      \
  if (!(expr)) {                                                         \
    throw InvariantException("Pre-condition Violation", #expr, message); \
  }

#define CHECK_INVARIANT(expr, message)                              \
  if (!(expr)) {                                                     \
    throw InvariantException("Invariant Violation", #expr, message); \
  }

static std::ostringstream rdWarningLog;
static std::ostringstream rdErrorLog;
#define BOOST_LOG(stream) stream

using INT_PAIR_VECT = std::vector<std::pair<int, int>>;
using UINT_VECT = std::vector<unsigned int>;

class MolSanitizeException final : public std::exception {
 public:
  explicit MolSanitizeException(std::string message)
      : message_(std::move(message)) {}
  const char *what() const noexcept override { return message_.c_str(); }

 private:
  std::string message_;
};

struct InchiToMolControl {
  std::vector<std::string> calls;
  std::vector<unsigned int> ranks;
  std::string fail_on;
};

static InchiToMolControl inchi_to_mol_control;

struct MolToInchiControl {
  bool active = false;
  bool needs_update_property_cache = false;
  std::vector<std::string> calls;
  std::string fail_on;
};

static MolToInchiControl mol_to_inchi_control;

static void record_mol_to_inchi_call(const char *call) {
  mol_to_inchi_control.calls.emplace_back(call);
  if (mol_to_inchi_control.fail_on == call) {
    throw MolSanitizeException(call);
  }
}

static void record_inchi_to_mol_call(const char *call) {
  inchi_to_mol_control.calls.emplace_back(call);
  if (inchi_to_mol_control.fail_on == call) {
    throw MolSanitizeException(call);
  }
}

class Atom {
 public:
  typedef enum {
    CHI_UNSPECIFIED = 0,
    CHI_TETRAHEDRAL_CW,
    CHI_TETRAHEDRAL_CCW,
    CHI_OTHER,
  } ChiralType;

  Atom() : Atom(0) {}

  explicit Atom(int atomic_number)
      : index_(0),
        atomic_number_(atomic_number),
        formal_charge_(0),
        num_explicit_hydrogens_(0),
        is_aromatic_(false),
        isotope_(0),
        num_radical_electrons_(0),
        no_implicit_(false),
        chiral_tag_(CHI_UNSPECIFIED),
        has_cip_rank_(false),
        cip_rank_(0),
        total_num_hydrogens_(0),
        total_degree_override_(-1),
        owning_mol_(nullptr) {}

  Atom(unsigned int index, int atomic_number, int formal_charge)
      : index_(index),
        atomic_number_(atomic_number),
        formal_charge_(formal_charge),
        num_explicit_hydrogens_(0),
        is_aromatic_(false),
        isotope_(0),
        num_radical_electrons_(0),
        no_implicit_(false),
        chiral_tag_(CHI_UNSPECIFIED),
        has_cip_rank_(false),
        cip_rank_(0),
        total_num_hydrogens_(0),
        total_degree_override_(-1),
        owning_mol_(nullptr) {}

  Atom(unsigned int index, int atomic_number, int formal_charge,
       unsigned int num_explicit_hydrogens)
      : index_(index),
        atomic_number_(atomic_number),
        formal_charge_(formal_charge),
        num_explicit_hydrogens_(num_explicit_hydrogens),
        is_aromatic_(false),
        isotope_(0),
        num_radical_electrons_(0),
        no_implicit_(false),
        chiral_tag_(CHI_UNSPECIFIED),
        has_cip_rank_(false),
        cip_rank_(0),
        total_num_hydrogens_(num_explicit_hydrogens),
        total_degree_override_(-1),
        owning_mol_(nullptr) {}

  unsigned int getIdx() const { return index_; }
  int getAtomicNum() const { return atomic_number_; }
  int getFormalCharge() const { return formal_charge_; }
  unsigned int getNumExplicitHs() const { return num_explicit_hydrogens_; }
  bool getIsAromatic() const { return is_aromatic_; }
  unsigned int getIsotope() const { return isotope_; }
  unsigned int getNumRadicalElectrons() const {
    return num_radical_electrons_;
  }
  unsigned int getTotalNumHs() const {
    record_mol_to_inchi_call("total_num_hydrogens");
    return total_num_hydrogens_;
  }
  unsigned int getTotalDegree() const;
  bool getNoImplicit() const { return no_implicit_; }
  ChiralType getChiralTag() const { return chiral_tag_; }
  bool hasCIPRank() const { return has_cip_rank_; }
  unsigned int getCIPRank() const { return cip_rank_; }
  double getMass() const;
  unsigned int getDegree() const;
  void setAtomicNum(int atomic_number) { atomic_number_ = atomic_number; }
  void setFormalCharge(int formal_charge) { formal_charge_ = formal_charge; }
  void setNumExplicitHs(unsigned int count) {
    num_explicit_hydrogens_ = count;
  }
  void setIsAromatic(bool is_aromatic) { is_aromatic_ = is_aromatic; }
  void setIsotope(unsigned int isotope) { isotope_ = isotope; }
  void setNumRadicalElectrons(unsigned int count) {
    num_radical_electrons_ = count;
  }
  void setTotalNumHs(unsigned int count) { total_num_hydrogens_ = count; }
  void setTotalDegree(unsigned int count) {
    total_degree_override_ = static_cast<int>(count);
  }
  void setNoImplicit(bool no_implicit) { no_implicit_ = no_implicit; }
  void setChiralTag(ChiralType chiral_tag) { chiral_tag_ = chiral_tag; }
  void setCIPRank(unsigned int rank) {
    has_cip_rank_ = true;
    cip_rank_ = rank;
  }
  bool invertChirality();
  int getPerturbationOrder(const std::list<int> &probe) const;
  int calcExplicitValence(bool strict);
  int calcImplicitValence() {
    record_mol_to_inchi_call("calc_implicit_valence");
    return 0;
  }
  void setIdx(unsigned int index) { index_ = index; }
  void setOwningMol(RWMol *molecule) { owning_mol_ = molecule; }

 private:
  unsigned int index_;
  int atomic_number_;
  int formal_charge_;
  unsigned int num_explicit_hydrogens_;
  bool is_aromatic_;
  unsigned int isotope_;
  unsigned int num_radical_electrons_;
  bool no_implicit_;
  ChiralType chiral_tag_;
  bool has_cip_rank_;
  unsigned int cip_rank_;
  unsigned int total_num_hydrogens_;
  int total_degree_override_;
  RWMol *owning_mol_;
};

class Bond {
 public:
  typedef enum {
    UNSPECIFIED = 0,
    SINGLE,
    DOUBLE,
    TRIPLE,
    QUADRUPLE,
    QUINTUPLE,
    HEXTUPLE,
    ONEANDAHALF,
    TWOANDAHALF,
    THREEANDAHALF,
    FOURANDAHALF,
    FIVEANDAHALF,
    AROMATIC,
    IONIC,
    HYDROGEN,
    THREECENTER,
    DATIVEONE,
    DATIVE,
    DATIVEL,
    DATIVER,
    OTHER,
    ZERO,
  } BondType;

  typedef enum {
    NONE = 0,
    BEGINWEDGE,
    BEGINDASH,
    ENDDOWNRIGHT,
    ENDUPRIGHT,
    EITHERDOUBLE,
    UNKNOWN,
  } BondDir;

  typedef enum {
    STEREONONE = 0,
    STEREOANY,
    STEREOZ,
    STEREOE,
    STEREOCIS,
    STEREOTRANS,
    STEREOATROPCW,
    STEREOATROPCCW,
  } BondStereo;

  Bond() : Bond(UNSPECIFIED) {}

  explicit Bond(BondType bond_type)
      : index_(0),
        begin_atom_index_(0),
        end_atom_index_(0),
        bond_type_(bond_type),
        direction_(NONE),
        is_aromatic_(false),
        stereo_(STEREONONE) {}

  Bond(unsigned int index, BondDir direction)
      : index_(index),
        begin_atom_index_(0),
        end_atom_index_(0),
        bond_type_(UNSPECIFIED),
        direction_(direction),
        is_aromatic_(false),
        stereo_(STEREONONE) {}

  Bond(unsigned int index, unsigned int begin_atom_index,
       unsigned int end_atom_index, BondType bond_type)
      : index_(index),
        begin_atom_index_(begin_atom_index),
        end_atom_index_(end_atom_index),
        bond_type_(bond_type),
        direction_(NONE),
        is_aromatic_(false),
        stereo_(STEREONONE) {}

  BondDir getBondDir() const { return direction_; }
  void setBondDir(BondDir direction) { direction_ = direction; }
  unsigned int getIdx() const { return index_; }
  void setIdx(unsigned int index) { index_ = index; }
  unsigned int getBeginAtomIdx() const { return begin_atom_index_; }
  unsigned int getEndAtomIdx() const { return end_atom_index_; }
  BondType getBondType() const { return bond_type_; }
  void setBondType(BondType bond_type) { bond_type_ = bond_type; }
  void setBeginAtomIdx(unsigned int index) { begin_atom_index_ = index; }
  void setEndAtomIdx(unsigned int index) { end_atom_index_ = index; }
  unsigned int getOtherAtomIdx(unsigned int index) const {
    return index == begin_atom_index_ ? end_atom_index_ : begin_atom_index_;
  }
  bool getIsAromatic() const { return is_aromatic_; }
  void setIsAromatic(bool aromatic) { is_aromatic_ = aromatic; }
  BondStereo getStereo() const { return stereo_; }
  void setStereo(BondStereo stereo) { stereo_ = stereo; }
  std::vector<unsigned int> &getStereoAtoms() { return stereo_atoms_; }
  const std::vector<unsigned int> &getStereoAtoms() const {
    return stereo_atoms_;
  }

 private:
  unsigned int index_;
  unsigned int begin_atom_index_;
  unsigned int end_atom_index_;
  BondType bond_type_;
  BondDir direction_;
  bool is_aromatic_;
  BondStereo stereo_;
  std::vector<unsigned int> stereo_atoms_;
};

struct GraphBond {
  unsigned int begin_atom_index;
  unsigned int end_atom_index;
  Bond::BondType bond_type;
  Bond::BondDir direction = Bond::NONE;
};

struct GraphAtom {
  int atomic_number;
  int formal_charge;
  unsigned int num_explicit_hydrogens;
  bool is_aromatic = false;
};

namespace RDGeom {
class Point3D {
 public:
  Point3D() : values_{0.0, 0.0, 0.0} {}
  Point3D(double x, double y, double z) : values_{x, y, z} {}
  double operator[](std::size_t index) const { return values_[index]; }

 private:
  double values_[3];
};
}  // namespace RDGeom

class Conformer {
 public:
  explicit Conformer(std::vector<RDGeom::Point3D> positions)
      : positions_(std::move(positions)) {}
  RDGeom::Point3D getAtomPos(unsigned int index) const {
    return positions_.at(index);
  }
  const std::vector<RDGeom::Point3D> &positions() const { return positions_; }

 private:
  std::vector<RDGeom::Point3D> positions_;
};

class RWMol {
 public:
  using ADJ_ITER = std::vector<unsigned int>::const_iterator;
  class AtomIterator {
   public:
    AtomIterator() = default;
    explicit AtomIterator(std::vector<Atom>::iterator iterator)
        : iterator_(iterator) {}
    Atom *operator*() const { return &*iterator_; }
    AtomIterator &operator++() {
      ++iterator_;
      return *this;
    }
    bool operator!=(const AtomIterator &other) const {
      return iterator_ != other.iterator_;
    }

   private:
    std::vector<Atom>::iterator iterator_;
  };

  RWMol() = default;

  RWMol(const RWMol &other)
      : atoms_(other.atoms_),
        bonds_(other.bonds_),
        adjacency_(other.adjacency_),
        conformers_(other.conformers_),
        needs_update_property_cache_(other.needs_update_property_cache_) {
    for (auto &atom : atoms_) atom.setOwningMol(this);
  }

  explicit RWMol(const std::vector<Bond::BondDir> &directions) {
    bonds_.reserve(directions.size());
    for (std::size_t index = 0; index < directions.size(); ++index) {
      bonds_.emplace_back(static_cast<unsigned int>(index), directions[index]);
    }
  }

  RWMol(const std::vector<std::pair<int, int>> &atoms,
        const std::vector<GraphBond> &bonds)
      : adjacency_(atoms.size()) {
    atoms_.reserve(atoms.size());
    for (std::size_t index = 0; index < atoms.size(); ++index) {
      atoms_.emplace_back(static_cast<unsigned int>(index), atoms[index].first,
                          atoms[index].second);
      atoms_.back().setOwningMol(this);
    }
    bonds_.reserve(bonds.size());
    for (std::size_t index = 0; index < bonds.size(); ++index) {
      const auto &bond = bonds[index];
      bonds_.emplace_back(static_cast<unsigned int>(index),
                          bond.begin_atom_index, bond.end_atom_index,
                          bond.bond_type);
      bonds_.back().setBondDir(bond.direction);
      adjacency_.at(bond.begin_atom_index).push_back(bond.end_atom_index);
      adjacency_.at(bond.end_atom_index).push_back(bond.begin_atom_index);
    }
  }

  RWMol(const std::vector<GraphAtom> &atoms,
        const std::vector<GraphBond> &bonds)
      : adjacency_(atoms.size()) {
    atoms_.reserve(atoms.size());
    for (std::size_t index = 0; index < atoms.size(); ++index) {
      atoms_.emplace_back(static_cast<unsigned int>(index),
                          atoms[index].atomic_number,
                          atoms[index].formal_charge,
                          atoms[index].num_explicit_hydrogens);
      atoms_.back().setIsAromatic(atoms[index].is_aromatic);
      atoms_.back().setOwningMol(this);
    }
    bonds_.reserve(bonds.size());
    for (std::size_t index = 0; index < bonds.size(); ++index) {
      const auto &bond = bonds[index];
      bonds_.emplace_back(static_cast<unsigned int>(index),
                          bond.begin_atom_index, bond.end_atom_index,
                          bond.bond_type);
      bonds_.back().setBondDir(bond.direction);
      adjacency_.at(bond.begin_atom_index).push_back(bond.end_atom_index);
      adjacency_.at(bond.end_atom_index).push_back(bond.begin_atom_index);
    }
  }

  unsigned int addAtom(Atom *atom, bool, bool take_ownership) {
    const auto index = static_cast<unsigned int>(atoms_.size());
    atoms_.push_back(*atom);
    atoms_.back().setIdx(index);
    atoms_.back().setOwningMol(this);
    adjacency_.emplace_back();
    if (take_ownership) {
      delete atom;
    }
    return index;
  }

  void addBond(unsigned int begin_atom_index, unsigned int end_atom_index,
               Bond::BondType bond_type) {
    const auto index = static_cast<unsigned int>(bonds_.size());
    bonds_.emplace_back(index, begin_atom_index, end_atom_index, bond_type);
    adjacency_.at(begin_atom_index).push_back(end_atom_index);
    adjacency_.at(end_atom_index).push_back(begin_atom_index);
  }

  unsigned int addBond(Bond *bond, bool take_ownership) {
    const auto index = static_cast<unsigned int>(bonds_.size());
    bonds_.push_back(*bond);
    bonds_.back().setIdx(index);
    adjacency_.at(bond->getBeginAtomIdx()).push_back(bond->getEndAtomIdx());
    adjacency_.at(bond->getEndAtomIdx()).push_back(bond->getBeginAtomIdx());
    if (take_ownership) {
      delete bond;
    }
    return index;
  }

  Bond *getBondWithIdx(unsigned int index) {
    if (index >= bonds_.size()) {
      throw RangeError(index, static_cast<unsigned int>(bonds_.size()));
    }
    return &bonds_[index];
  }

  Atom *getAtomWithIdx(unsigned int index) { return &atoms_.at(index); }

  Bond *getBondBetweenAtoms(unsigned int first, unsigned int second) {
    for (auto &bond : bonds_) {
      if ((bond.getBeginAtomIdx() == first && bond.getEndAtomIdx() == second) ||
          (bond.getBeginAtomIdx() == second && bond.getEndAtomIdx() == first)) {
        return &bond;
      }
    }
    return nullptr;
  }

  const Bond *getBondBetweenAtoms(unsigned int first,
                                  unsigned int second) const {
    for (const auto &bond : bonds_) {
      if ((bond.getBeginAtomIdx() == first &&
           bond.getEndAtomIdx() == second) ||
          (bond.getBeginAtomIdx() == second &&
           bond.getEndAtomIdx() == first)) {
        return &bond;
      }
    }
    return nullptr;
  }

  boost::tuple<ADJ_ITER, ADJ_ITER> getAtomNeighbors(const Atom *atom) const {
    const auto &neighbors = adjacency_.at(atom->getIdx());
    return boost::make_tuple(neighbors.begin(), neighbors.end());
  }

  AtomIterator beginAtoms() { return AtomIterator(atoms_.begin()); }
  AtomIterator endAtoms() { return AtomIterator(atoms_.end()); }
  std::vector<Bond *> atomBonds(const Atom *atom) {
    std::vector<Bond *> result;
    for (auto &bond : bonds_) {
      if (bond.getBeginAtomIdx() == atom->getIdx() ||
          bond.getEndAtomIdx() == atom->getIdx()) {
        result.push_back(&bond);
      }
    }
    return result;
  }
  void updatePropertyCache(bool strict) {
    if (strict) {
      throw MolSanitizeException("strict property cache not supported");
    }
    if (mol_to_inchi_control.active) {
      record_mol_to_inchi_call("update_property_cache");
    } else {
      record_inchi_to_mol_call("update_property_cache");
    }
  }
  bool needsUpdatePropertyCache() const {
    if (mol_to_inchi_control.active) {
      record_mol_to_inchi_call("needs_update_property_cache");
    }
    return needs_update_property_cache_;
  }
  void setNeedsUpdatePropertyCache(bool value) {
    needs_update_property_cache_ = value;
  }
  void setConformers(
      const std::vector<std::vector<RDGeom::Point3D>> &conformers) {
    conformers_.clear();
    for (const auto &positions : conformers) {
      conformers_.push_back(std::make_shared<Conformer>(positions));
    }
  }
  unsigned int getNumConformers() const {
    return static_cast<unsigned int>(conformers_.size());
  }
  auto beginConformers() { return conformers_.begin(); }
  unsigned int atomDegree(unsigned int index) const {
    return static_cast<unsigned int>(adjacency_.at(index).size());
  }

  const std::vector<Atom> &atoms() const { return atoms_; }
  const std::vector<Bond> &bonds() const { return bonds_; }
  unsigned int getNumAtoms() const {
    return static_cast<unsigned int>(atoms_.size());
  }
  unsigned int getNumBonds() const {
    return static_cast<unsigned int>(bonds_.size());
  }
  const std::vector<std::shared_ptr<Conformer>> &conformers() const {
    return conformers_;
  }

 private:
  std::vector<Atom> atoms_;
  std::vector<Bond> bonds_;
  std::vector<std::vector<unsigned int>> adjacency_;
  std::vector<std::shared_ptr<Conformer>> conformers_;
  bool needs_update_property_cache_ = false;
};

double Atom::getMass() const {
  record_inchi_to_mol_call("average_atomic_weight");
  switch (atomic_number_) {
    case 1:
      return 1.008;
    case 6:
      return 12.011;
    case 7:
      return 14.007;
    case 8:
      return 15.999;
    case 16:
      return 32.06;
    case 17:
      return 35.45;
    case 35:
      return 79.904;
    default:
      throw InvariantException("Pre-condition Violation", "atomic_number",
                               "Atomic number not found");
  }
}

bool Atom::invertChirality() {
  if (chiral_tag_ == CHI_TETRAHEDRAL_CW) {
    chiral_tag_ = CHI_TETRAHEDRAL_CCW;
    return true;
  }
  if (chiral_tag_ == CHI_TETRAHEDRAL_CCW) {
    chiral_tag_ = CHI_TETRAHEDRAL_CW;
    return true;
  }
  return false;
}

int Atom::getPerturbationOrder(const std::list<int> &probe) const {
  std::vector<int> reference;
  for (const auto &bond : owning_mol_->bonds()) {
    if (bond.getBeginAtomIdx() == index_ || bond.getEndAtomIdx() == index_) {
      reference.push_back(static_cast<int>(bond.getIdx()));
    }
  }
  std::vector<int> permutation(probe.begin(), probe.end());
  int swaps = 0;
  for (std::size_t position = 0; position < reference.size(); ++position) {
    if (permutation.at(position) == reference[position]) {
      continue;
    }
    auto found = std::find(permutation.begin() + position, permutation.end(),
                           reference[position]);
    CHECK_INVARIANT(found != permutation.end(), "invalid perturbation order");
    std::iter_swap(permutation.begin() + position, found);
    ++swaps;
  }
  return swaps;
}

class PeriodicTable {
 public:
  static PeriodicTable *getTable() {
    static PeriodicTable table;
    return &table;
  }

  int getAtomicNumber(const char *element) const {
    record_inchi_to_mol_call("atomic_number");
    const std::string value(element);
    if (value == "H") return 1;
    if (value == "C") return 6;
    if (value == "N") return 7;
    if (value == "O") return 8;
    if (value == "S") return 16;
    if (value == "Cl") return 17;
    if (value == "Br") return 35;
    throw InvariantException("Post-condition Violation", "atomic_number",
                             "Element '" + value + "' not found");
  }

  std::string getElementSymbol(unsigned int atomic_number) const {
    record_mol_to_inchi_call("element_symbol");
    switch (atomic_number) {
      case 1:
        return "H";
      case 6:
        return "C";
      case 7:
        return "N";
      case 8:
        return "O";
      case 9:
        return "F";
      case 15:
        return "P";
      case 17:
        return "Cl";
      case 35:
        return "Br";
      case 53:
        return "I";
      default:
        throw InvariantException("Pre-condition Violation", "atomic_number",
                                 "Atomic number not found");
    }
  }

  double getAtomicWeight(unsigned int atomic_number) const {
    record_mol_to_inchi_call("atomic_weight");
    switch (atomic_number) {
      case 1:
        return 1.008;
      case 6:
        return 12.011;
      case 7:
        return 14.007;
      case 8:
        return 15.999;
      case 9:
        return 18.998;
      case 15:
        return 30.5;
      case 17:
        return 35.45;
      case 35:
        return 79.904;
      case 53:
        return 126.904;
      default:
        throw InvariantException("Pre-condition Violation", "atomic_number",
                                 "Atomic number not found");
    }
  }
};

namespace Chirality {
static void assignAtomCIPRanks(RWMol &molecule,
                               std::vector<unsigned int> &ranks) {
  record_inchi_to_mol_call("assign_atom_cip_ranks");
  if (inchi_to_mol_control.ranks.empty()) {
    ranks.resize(molecule.getNumAtoms());
    for (unsigned int index = 0; index < ranks.size(); ++index) {
      ranks[index] = index;
    }
  } else {
    ranks = inchi_to_mol_control.ranks;
  }
  for (unsigned int index = 0;
       index < molecule.getNumAtoms() && index < ranks.size(); ++index) {
    molecule.getAtomWithIdx(index)->setCIPRank(ranks[index]);
  }
}
}  // namespace Chirality

namespace MolOps {
static void Kekulize(RWMol &, bool mark_atoms_bonds) {
  CHECK_INVARIANT(!mark_atoms_bonds, "markAtomsBonds must be false");
  record_mol_to_inchi_call("kekulize");
}
static void removeHs(RWMol &) { record_inchi_to_mol_call("remove_hydrogens"); }
static void sanitizeMol(RWMol &) {
  record_inchi_to_mol_call("sanitize_molecule");
}
static void assignStereochemistry(RWMol &, bool clean_it, bool force) {
  CHECK_INVARIANT(clean_it, "cleanIt must be true");
  CHECK_INVARIANT(force, "force must be true");
  record_inchi_to_mol_call("assign_stereochemistry");
}
}  // namespace MolOps

struct ExtraInchiReturnValues {
  int returnCode = 0;
  std::string messagePtr;
  std::string logPtr;
  std::string auxInfoPtr;
};

int Atom::calcExplicitValence(bool) {
  double valence = static_cast<double>(num_explicit_hydrogens_);
  for (const auto &bond : owning_mol_->bonds()) {
    if (bond.getBeginAtomIdx() != index_ && bond.getEndAtomIdx() != index_) {
      continue;
    }
    switch (bond.getBondType()) {
      case Bond::UNSPECIFIED:
      case Bond::IONIC:
      case Bond::HYDROGEN:
      case Bond::ZERO:
        break;
      case Bond::SINGLE:
        valence += 1.0;
        break;
      case Bond::DOUBLE:
        valence += 2.0;
        break;
      case Bond::TRIPLE:
        valence += 3.0;
        break;
      case Bond::QUADRUPLE:
        valence += 4.0;
        break;
      case Bond::QUINTUPLE:
        valence += 5.0;
        break;
      case Bond::HEXTUPLE:
        valence += 6.0;
        break;
      case Bond::ONEANDAHALF:
      case Bond::AROMATIC:
        valence += 1.5;
        break;
      case Bond::TWOANDAHALF:
        valence += 2.5;
        break;
      case Bond::THREEANDAHALF:
        valence += 3.5;
        break;
      case Bond::FOURANDAHALF:
        valence += 4.5;
        break;
      case Bond::FIVEANDAHALF:
        valence += 5.5;
        break;
      case Bond::DATIVEONE:
      case Bond::DATIVE:
        if (bond.getEndAtomIdx() == index_) {
          valence += 1.0;
        }
        break;
      case Bond::THREECENTER:
      case Bond::DATIVEL:
      case Bond::DATIVER:
      case Bond::OTHER:
        throw ValueErrorException(bond.getIdx(),
                                  static_cast<int>(bond.getBondType()));
    }
  }
  return static_cast<int>(std::round(valence + 0.1));
}

unsigned int Atom::getDegree() const { return owning_mol_->atomDegree(index_); }

unsigned int Atom::getTotalDegree() const {
  record_mol_to_inchi_call("total_degree");
  if (total_degree_override_ >= 0) {
    return static_cast<unsigned int>(total_degree_override_);
  }
  return getDegree() + total_num_hydrogens_;
}

using MatchVectType = std::vector<std::pair<int, int>>;

static void SubstructMatch(RWMol &molecule, const RWMol &query,
                           std::vector<MatchVectType> &matches) {
  matches.clear();
  std::set<std::vector<unsigned int>> seen_atom_sets;
  std::vector<unsigned int> mapping(query.getNumAtoms(), 0);
  std::vector<bool> used(molecule.getNumAtoms(), false);

  std::function<void(unsigned int)> visit = [&](unsigned int query_index) {
    if (matches.size() == 1000) {
      return;
    }
    if (query_index == query.getNumAtoms()) {
      auto atom_set = mapping;
      std::sort(atom_set.begin(), atom_set.end());
      if (!seen_atom_sets.insert(atom_set).second) {
        return;
      }
      MatchVectType match;
      match.reserve(mapping.size());
      for (unsigned int index = 0; index < mapping.size(); ++index) {
        match.emplace_back(static_cast<int>(index),
                           static_cast<int>(mapping[index]));
      }
      matches.push_back(std::move(match));
      return;
    }

    const auto &query_atom = query.atoms()[query_index];
    for (unsigned int atom_index = 0; atom_index < molecule.getNumAtoms();
         ++atom_index) {
      if (used[atom_index]) {
        continue;
      }
      const auto &atom = molecule.atoms()[atom_index];
      if (query_atom.getAtomicNum() != atom.getAtomicNum() ||
          (query_atom.getFormalCharge() != 0 &&
           query_atom.getFormalCharge() != atom.getFormalCharge())) {
        continue;
      }

      bool compatible = true;
      for (unsigned int previous = 0; previous < query_index; ++previous) {
        const auto *query_bond =
            query.getBondBetweenAtoms(query_index, previous);
        if (query_bond == nullptr) {
          continue;
        }
        const auto *molecule_bond =
            molecule.getBondBetweenAtoms(atom_index, mapping[previous]);
        if (molecule_bond == nullptr ||
            (query_bond->getBondType() != Bond::UNSPECIFIED &&
             molecule_bond->getBondType() != query_bond->getBondType())) {
          compatible = false;
          break;
        }
      }
      if (!compatible) {
        continue;
      }

      mapping[query_index] = atom_index;
      used[atom_index] = true;
      visit(query_index + 1);
      used[atom_index] = false;
    }
  };
  visit(0);
}

using ROMol = RWMol;

static RWMol *SmilesToMol(const std::string &smiles) {
  CHECK_INVARIANT(smiles == "[O-][Cl+3]([O-])([O-])O",
                  "unexpected rCleanUp query");
  return new RWMol(
      std::vector<std::pair<int, int>>{{8, -1}, {17, 3}, {8, -1},
                                      {8, -1}, {8, 0}},
      std::vector<GraphBond>{{0, 1, Bond::SINGLE}, {1, 2, Bond::SINGLE},
                             {1, 3, Bond::SINGLE}, {1, 4, Bond::SINGLE}});
}

namespace {
#include ASSIGN_BOND_DIRS_SOURCE
#include FIND_ALTERNATING_BONDS_SOURCE
#include NEIGHBORING_SI_SOURCE
#include VALENCE4N_CLEANUP1_SOURCE
#include VALENCE4N_CLEANUP2_SOURCE
#include VALENCE5N_CLEANUP1_SOURCE
#include VALENCE5N_CLEANUP2_SOURCE
#include VALENCE5N_CLEANUP3_SOURCE
#include VALENCE5N_CLEANUP4_SOURCE
#include VALENCE5N_CLEANUP5_SOURCE
#include VALENCE5N_CLEANUP6_SOURCE
#pragma GCC diagnostic push
// GCC -O2 cannot prove that the official substructure match initializes map.
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include VALENCE5N_CLEANUP7_SOURCE
#pragma GCC diagnostic pop
#include VALENCE5N_CLEANUP8_SOURCE
#include VALENCE5N_CLEANUP9_SOURCE
#include VALENCE5N_CLEANUPA_SOURCE
#include VALENCE5N_CLEANUPB_SOURCE
#include VALENCE7S_CLEANUP1_SOURCE
#include VALENCE7S_CLEANUP2_SOURCE
#include VALENCE7S_CLEANUP3_SOURCE
#include VALENCE8S_CLEANUP1_SOURCE
#include VALENCE8CL_CLEANUP1_SOURCE
#include VALENCE5CL_CLEANUP1_SOURCE
#include VALENCE3CL_CLEANUP1_SOURCE
#include CLEAN_UP_SOURCE
}  // namespace

#include INCHI_TO_MOL_SOURCE
#include FIX_OPTION_SYMBOL_SOURCE
#include R_CLEAN_UP_SOURCE
#include MOL_TO_INCHI_SOURCE
#include INCHI_TO_INCHI_KEY_SOURCE
#include MOL_TO_INCHI_KEY_SOURCE

}  // namespace RDKit

struct AssignCase {
  const char *case_id;
  std::vector<RDKit::Bond::BondDir> initial_directions;
  std::vector<std::pair<int, int>> z_pairs;
  std::vector<std::pair<int, int>> e_pairs;
};

struct SearchCase {
  const char *case_id;
  std::vector<std::pair<int, int>> atoms;
  std::vector<RDKit::GraphBond> bonds;
  unsigned int current_atom_index;
  int desired_atomic_number;
  int desired_atom_charge;
  RDKit::Bond::BondType desired_next_bond_type;
  RDKit::Bond::BondType desired_ending_bond_type;
  unsigned int current_path_length;
  unsigned int max_path_length;
  int last_bond_index;
  std::vector<unsigned int> initial_path;
  std::set<int> initial_visited;
};

struct NeighboringSiCase {
  const char *case_id;
  std::vector<std::pair<int, int>> atoms;
  std::vector<RDKit::GraphBond> bonds;
  unsigned int atom_index;
};

struct Valence4NCleanup1Case {
  const char *case_id;
  std::vector<RDKit::GraphAtom> atoms;
  std::vector<RDKit::GraphBond> bonds;
  unsigned int atom_index;
};

struct Valence4NCleanup2Case {
  const char *case_id;
  std::vector<RDKit::GraphAtom> atoms;
  std::vector<RDKit::GraphBond> bonds;
  unsigned int atom_index;
};

using Valence5NCleanup1Case = Valence4NCleanup2Case;
using Valence5NCleanup2Case = Valence4NCleanup2Case;
using Valence5NCleanup3Case = Valence4NCleanup2Case;
using Valence5NCleanup4Case = Valence4NCleanup2Case;

struct Valence5NCleanup5Case {
  const char *case_id;
  std::vector<RDKit::GraphAtom> atoms;
  std::vector<RDKit::GraphBond> bonds;
  unsigned int atom_index;
  int atomic_number;
};

struct Valence5NCleanup6Case {
  std::string case_id;
  std::vector<RDKit::GraphAtom> atoms;
  std::vector<RDKit::GraphBond> bonds;
  unsigned int atom_index;
};

using Valence5NCleanup7Case = Valence5NCleanup6Case;
using Valence5NCleanup8Case = Valence5NCleanup6Case;
using Valence5NCleanup9Case = Valence5NCleanup6Case;
using Valence5NCleanupACase = Valence5NCleanup6Case;
using Valence5NCleanupBCase = Valence5NCleanup6Case;
using Valence7SCleanup1Case = Valence5NCleanup6Case;
using Valence7SCleanup2Case = Valence5NCleanup6Case;
using Valence7SCleanup3Case = Valence5NCleanup6Case;
using Valence8SCleanup1Case = Valence5NCleanup6Case;
using Valence8ClCleanup1Case = Valence5NCleanup6Case;
using Valence5ClCleanup1Case = Valence5NCleanup6Case;
using Valence3ClCleanup1Case = Valence5NCleanup6Case;

struct CleanUpCase {
  std::string case_id;
  std::vector<RDKit::GraphAtom> atoms;
  std::vector<RDKit::GraphBond> bonds;
};

template <typename T>
static void print_numbers(const std::vector<T> &values) {
  std::cout << '[';
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (index != 0) {
      std::cout << ',';
    }
    std::cout << values[index];
  }
  std::cout << ']';
}

static void print_set(const std::set<int> &values) {
  std::cout << '[';
  std::size_t index = 0;
  for (int value : values) {
    if (index++ != 0) {
      std::cout << ',';
    }
    std::cout << value;
  }
  std::cout << ']';
}

static void print_directions(const std::vector<RDKit::Bond::BondDir> &values) {
  std::cout << '[';
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (index != 0) {
      std::cout << ',';
    }
    std::cout << static_cast<int>(values[index]);
  }
  std::cout << ']';
}

static void print_pairs(const std::vector<std::pair<int, int>> &pairs) {
  std::cout << '[';
  for (std::size_t index = 0; index < pairs.size(); ++index) {
    if (index != 0) {
      std::cout << ',';
    }
    std::cout << '[' << pairs[index].first << ',' << pairs[index].second << ']';
  }
  std::cout << ']';
}

static std::vector<RDKit::Bond::BondDir> bond_directions(
    const RDKit::RWMol &molecule) {
  std::vector<RDKit::Bond::BondDir> result;
  result.reserve(molecule.bonds().size());
  for (const auto &bond : molecule.bonds()) {
    result.push_back(bond.getBondDir());
  }
  return result;
}

static void print_assign_bond_fields(const RDKit::RWMol &molecule) {
  std::cout << '[';
  for (std::size_t index = 0; index < molecule.bonds().size(); ++index) {
    if (index != 0) {
      std::cout << ',';
    }
    const auto &bond = molecule.bonds()[index];
    std::cout << "{\"index\":" << bond.getIdx() << ",\"direction\":"
              << static_cast<int>(bond.getBondDir()) << '}';
  }
  std::cout << ']';
}

static void print_atom_fields(const RDKit::RWMol &molecule) {
  std::cout << '[';
  for (std::size_t index = 0; index < molecule.atoms().size(); ++index) {
    if (index != 0) {
      std::cout << ',';
    }
    const auto &atom = molecule.atoms()[index];
    std::cout << "{\"index\":" << atom.getIdx()
              << ",\"atomic_number\":" << atom.getAtomicNum()
              << ",\"formal_charge\":" << atom.getFormalCharge()
              << ",\"num_explicit_hydrogens\":"
              << atom.getNumExplicitHs() << ",\"is_aromatic\":"
              << (atom.getIsAromatic() ? "true" : "false") << '}';
  }
  std::cout << ']';
}

static void print_graph_bond_fields(const RDKit::RWMol &molecule) {
  std::cout << '[';
  for (std::size_t index = 0; index < molecule.bonds().size(); ++index) {
    if (index != 0) {
      std::cout << ',';
    }
    const auto &bond = molecule.bonds()[index];
    std::cout << "{\"index\":" << bond.getIdx()
              << ",\"begin_atom_index\":" << bond.getBeginAtomIdx()
              << ",\"end_atom_index\":" << bond.getEndAtomIdx()
              << ",\"bond_type\":" << static_cast<int>(bond.getBondType())
              << ",\"direction\":" << static_cast<int>(bond.getBondDir())
              << '}';
  }
  std::cout << ']';
}

static auto atom_snapshot(const RDKit::RWMol &molecule) {
  std::vector<std::tuple<unsigned int, int, int, unsigned int, bool>> result;
  for (const auto &atom : molecule.atoms()) {
    result.emplace_back(atom.getIdx(), atom.getAtomicNum(),
                        atom.getFormalCharge(), atom.getNumExplicitHs(),
                        atom.getIsAromatic());
  }
  return result;
}

static auto bond_snapshot(const RDKit::RWMol &molecule) {
  std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, int>>
      result;
  for (const auto &bond : molecule.bonds()) {
    result.emplace_back(bond.getIdx(), bond.getBeginAtomIdx(),
                        bond.getEndAtomIdx(),
                        static_cast<int>(bond.getBondType()),
                        static_cast<int>(bond.getBondDir()));
  }
  return result;
}

static std::vector<unsigned int> stack_bottom_to_top(
    std::stack<RDKit::Bond *> path) {
  std::vector<unsigned int> result;
  while (!path.empty()) {
    result.insert(result.begin(), path.top()->getIdx());
    path.pop();
  }
  return result;
}

static void print_assign_case(const AssignCase &test_case) {
  RDKit::RWMol molecule(test_case.initial_directions);
  auto z_pairs = test_case.z_pairs;
  auto e_pairs = test_case.e_pairs;
  const auto z_before = z_pairs;
  const auto e_before = e_pairs;

  bool returned = false;
  bool result = false;
  bool threw = false;
  unsigned int error_index = 0;
  unsigned int error_upper_bound = 0;
  try {
    result = RDKit::assignBondDirs(molecule, z_pairs, e_pairs);
    returned = true;
  } catch (const RDKit::RangeError &error) {
    threw = true;
    error_index = error.index;
    error_upper_bound = error.upper_bound;
  }

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"operation\":\"assignBondDirs\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"initial_directions\":";
  print_directions(test_case.initial_directions);
  std::cout << ",\"z_pairs\":";
  print_pairs(test_case.z_pairs);
  std::cout << ",\"e_pairs\":";
  print_pairs(test_case.e_pairs);
  std::cout << "},\"output\":{\"status\":\""
            << (returned ? "return" : "exception") << "\",\"result\":";
  if (returned) {
    std::cout << (result ? "true" : "false");
  } else {
    std::cout << "null";
  }
  std::cout << ",\"exception\":";
  if (threw) {
    std::cout << "{\"kind\":\"Range Error\",\"expression\":\"idx\",\"index\":"
              << error_index << ",\"upper_bound\":" << error_upper_bound
              << ",\"detail\":\"" << error_index << " < " << error_upper_bound
              << "\"}";
  } else {
    std::cout << "null";
  }
  std::cout << ",\"diagnostics\":";
  if (threw) {
    std::cout << "[{\"kind\":\"Range Error\",\"expression\":\"idx\",\"detail\":\""
              << error_index << " < " << error_upper_bound << "\"}]";
  } else {
    std::cout << "[]";
  }
  std::cout << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"bond_fields\":";
  print_assign_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"atom_fields\":[],\"properties\":[],\"z_pairs_unchanged\":"
            << (z_pairs == z_before ? "true" : "false")
            << ",\"e_pairs_unchanged\":"
            << (e_pairs == e_before ? "true" : "false") << "}}\n";
}

static void print_search_case(const SearchCase &test_case) {
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);
  std::stack<RDKit::Bond *> path;
  for (unsigned int bond_index : test_case.initial_path) {
    path.push(molecule.getBondWithIdx(bond_index));
  }
  auto visited = test_case.initial_visited;
  RDKit::Bond *last_bond = test_case.last_bond_index < 0
                              ? nullptr
                              : molecule.getBondWithIdx(static_cast<unsigned int>(
                                    test_case.last_bond_index));
  RDKit::Atom *target = RDKit::findAlternatingBonds(
      molecule, molecule.getAtomWithIdx(test_case.current_atom_index),
      test_case.desired_atomic_number, test_case.desired_atom_charge,
      test_case.desired_next_bond_type, test_case.desired_ending_bond_type,
      test_case.current_path_length, test_case.max_path_length, last_bond, path,
      visited);
  const auto output_path = stack_bottom_to_top(path);
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"ef0a05d3a018d6928d66af0d5238dea5bb22742e233f1fe8229f62578fe7d003\","
               "\"operation\":\"findAlternatingBonds\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"current_atom_index\":" << test_case.current_atom_index
            << ",\"desired_atomic_number\":"
            << test_case.desired_atomic_number
            << ",\"desired_atom_charge\":" << test_case.desired_atom_charge
            << ",\"desired_next_bond_type\":"
            << static_cast<int>(test_case.desired_next_bond_type)
            << ",\"desired_ending_bond_type\":"
            << static_cast<int>(test_case.desired_ending_bond_type)
            << ",\"current_path_length\":" << test_case.current_path_length
            << ",\"max_path_length\":" << test_case.max_path_length
            << ",\"last_bond_index\":";
  if (test_case.last_bond_index < 0) {
    std::cout << "null";
  } else {
    std::cout << test_case.last_bond_index;
  }
  std::cout << ",\"initial_path\":";
  print_numbers(test_case.initial_path);
  std::cout << ",\"initial_visited\":";
  print_set(test_case.initial_visited);
  std::cout << "},\"output\":{\"status\":\"return\",\"target_atom_index\":";
  if (target == nullptr) {
    std::cout << "null";
  } else {
    std::cout << target->getIdx();
  }
  std::cout << ",\"path\":";
  print_numbers(output_path);
  std::cout << ",\"visited\":";
  print_set(visited);
  std::cout << ",\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[],\"diagnostics\":[],\"exception\":null}}\n";
}

static int assign_bond_dirs_records() {
  using D = RDKit::Bond::BondDir;
  const std::vector<AssignCase> cases = {
      {"no-rules-all-directions",
       {D::NONE, D::BEGINWEDGE, D::BEGINDASH, D::ENDDOWNRIGHT,
        D::ENDUPRIGHT, D::EITHERDOUBLE, D::UNKNOWN},
       {},
       {}},
      {"z-reversed-rules", {D::NONE, D::NONE, D::NONE}, {{2, 0}, {0, 1}}, {}},
      {"e-reversed-rules", {D::NONE, D::NONE, D::NONE}, {}, {{2, 0}, {0, 1}}},
      {"disconnected-smallest-seed",
       {D::NONE, D::NONE, D::NONE, D::NONE, D::NONE, D::NONE},
       {{5, 4}, {3, 2}},
       {{1, 0}}},
      {"duplicate-rules", {D::NONE, D::NONE}, {{0, 1}, {1, 0}, {0, 1}}, {}},
      {"self-pair", {D::NONE, D::NONE, D::NONE}, {{1, 1}}, {}},
      {"same-pair-z-and-e-conflict", {D::NONE, D::NONE}, {{0, 1}}, {{0, 1}}},
      {"delayed-cycle-conflict",
       {D::NONE, D::NONE, D::NONE},
       {{0, 1}, {1, 2}},
       {{0, 2}}},
      {"existing-seed-direction-conflict", {D::BEGINWEDGE}, {{0, 0}}, {}},
      {"existing-other-z-conflict", {D::NONE, D::ENDDOWNRIGHT}, {{0, 1}}, {}},
      {"existing-other-e-then-seed-conflict",
       {D::NONE, D::ENDDOWNRIGHT},
       {},
       {{0, 1}}},
      {"negative-index-exception", {D::NONE}, {{-1, -1}}, {}},
      {"partial-mutation-index-exception", {D::NONE, D::NONE}, {{0, 3}}, {}},
  };
  for (const auto &test_case : cases) {
    print_assign_case(test_case);
  }
  return 0;
}

static int find_alternating_bonds_records() {
  using T = RDKit::Bond::BondType;
  const unsigned int umax = std::numeric_limits<unsigned int>::max();
  const std::vector<SearchCase> cases = {
      {"root-clears-state-at-cutoff", {{6, 0}, {6, 0}}, {{0, 1, T::SINGLE}}, 0, 7, 0, T::SINGLE, T::SINGLE, 0, 0, -1, {0, 0}, {7, 8}},
      {"direct-target-replaces-longer-path", {{6, 0}, {7, 0}}, {{0, 1, T::SINGLE}}, 1, 7, 0, T::DOUBLE, T::SINGLE, 1, 0, 0, {0, 0}, {0}},
      {"direct-target-rejects-equal-path", {{6, 0}, {7, 0}}, {{0, 1, T::SINGLE}}, 1, 7, 0, T::DOUBLE, T::SINGLE, 1, 9, 0, {0}, {0}},
      {"atomic-number-miss", {{6, 0}, {6, 0}}, {{0, 1, T::SINGLE}}, 1, 7, 0, T::DOUBLE, T::SINGLE, 1, 1, 0, {}, {}},
      {"formal-charge-miss", {{6, 0}, {7, 1}}, {{0, 1, T::SINGLE}}, 1, 7, 0, T::DOUBLE, T::SINGLE, 1, 1, 0, {}, {}},
      {"ending-bond-type-miss", {{6, 0}, {7, 0}}, {{0, 1, T::DOUBLE}}, 1, 7, 0, T::SINGLE, T::SINGLE, 1, 1, 0, {}, {}},
      {"single-double-single-alternation", {{6, 0}, {6, 0}, {6, 0}, {7, 0}}, {{0, 1, T::SINGLE}, {1, 2, T::DOUBLE}, {2, 3, T::SINGLE}}, 0, 7, 0, T::SINGLE, T::SINGLE, 0, 4, -1, {}, {}},
      {"double-start-defaults-next-to-single", {{6, 0}, {6, 0}, {7, 0}}, {{0, 1, T::DOUBLE}, {1, 2, T::SINGLE}}, 0, 7, 0, T::DOUBLE, T::SINGLE, 0, 2, -1, {}, {}},
      {"triple-then-single-special-case", {{7, 0}, {6, -1}, {7, -1}}, {{0, 1, T::TRIPLE}, {1, 2, T::SINGLE}}, 0, 7, -1, T::TRIPLE, T::SINGLE, 0, 2, -1, {}, {}},
      {"nonalternating-ending-success", {{6, 0}, {7, 0}}, {{0, 1, T::TRIPLE}}, 0, 7, 0, T::DOUBLE, T::TRIPLE, 0, 3, -1, {}, {}},
      {"nonalternating-ending-target-miss", {{6, 0}, {6, 0}, {7, 0}}, {{0, 1, T::TRIPLE}, {1, 2, T::SINGLE}}, 0, 7, 0, T::DOUBLE, T::TRIPLE, 0, 3, -1, {}, {}},
      {"single-double-ending-disables-special-branch", {{6, 0}, {7, 0}}, {{0, 1, T::DOUBLE}}, 0, 7, 0, T::SINGLE, T::DOUBLE, 0, 2, -1, {}, {}},
      {"cycle-global-visited", {{6, 0}, {6, 0}, {6, 0}}, {{0, 1, T::SINGLE}, {1, 2, T::DOUBLE}, {2, 0, T::SINGLE}}, 0, 7, 0, T::SINGLE, T::SINGLE, 0, 9, -1, {}, {}},
      {"shorter-path-found-later", {{6, 0}, {6, 0}, {6, 0}, {7, 0}, {7, 0}}, {{0, 1, T::SINGLE}, {1, 2, T::DOUBLE}, {2, 3, T::SINGLE}, {0, 4, T::SINGLE}}, 0, 7, 0, T::SINGLE, T::SINGLE, 0, 4, -1, {}, {}},
      {"equal-path-keeps-first-insertion", {{6, 0}, {7, 0}, {7, 0}}, {{0, 1, T::SINGLE}, {0, 2, T::SINGLE}}, 0, 7, 0, T::SINGLE, T::SINGLE, 0, 1, -1, {}, {}},
      {"equal-path-reversed-insertion", {{6, 0}, {7, 0}, {7, 0}}, {{0, 2, T::SINGLE}, {0, 1, T::SINGLE}}, 0, 7, 0, T::SINGLE, T::SINGLE, 0, 1, -1, {}, {}},
      {"nonroot-preserves-visited-and-skips-neighbor", {{6, 0}, {7, 0}, {6, 0}}, {{0, 1, T::SINGLE}, {0, 2, T::DOUBLE}}, 0, 7, 0, T::SINGLE, T::SINGLE, 0, 1, 1, {1}, {1}},
      {"negative-atomic-number", {{6, 0}, {-1, 0}}, {{0, 1, T::SINGLE}}, 0, -1, 0, T::SINGLE, T::SINGLE, 0, 1, -1, {}, {}},
      {"near-unsigned-path-limit", {{6, 0}, {7, 0}}, {{0, 1, T::SINGLE}}, 0, 7, 0, T::SINGLE, T::SINGLE, umax - 1, umax, -1, {}, {}},
      {"unmatched-edge-type", {{6, 0}, {7, 0}}, {{0, 1, T::IONIC}}, 0, 7, 0, T::SINGLE, T::SINGLE, 0, 3, -1, {}, {}},
  };
  for (const auto &test_case : cases) {
    print_search_case(test_case);
  }
  return 0;
}

static void print_neighboring_si_case(const NeighboringSiCase &test_case) {
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);
  const int count = RDKit::getNumDoubleBondedNegativelyChargedNeighboringSi(
      molecule, molecule.getAtomWithIdx(test_case.atom_index));
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"a5d70ac89adb467c691ba0a8a16878e951833d598c7f7d15a1a9d241f0c17374\","
               "\"operation\":\"getNumDoubleBondedNegativelyChargedNeighboringSi\","
               "\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"atom_index\":" << test_case.atom_index
            << "},\"output\":{\"status\":\"return\",\"count\":" << count
            << ",\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[],\"diagnostics\":[],\"exception\":null}}\n";
}

static int neighboring_si_records() {
  using T = RDKit::Bond::BondType;
  const std::vector<NeighboringSiCase> cases = {
      {"isolated", {{6, 0}}, {}, 0},
      {"single-match-reversed-endpoints",
       {{6, 0}, {14, -1}},
       {{1, 0, T::DOUBLE}},
       0},
      {"each-condition-and-two-matches",
       {{6, 0}, {14, -1}, {14, 0}, {13, -1}, {14, -1}, {14, -1}, {14, -1}},
       {{0, 1, T::DOUBLE}, {0, 2, T::DOUBLE}, {0, 3, T::DOUBLE},
        {0, 4, T::SINGLE}, {5, 0, T::DOUBLE}, {0, 6, T::TRIPLE}},
       0},
      {"all-miss",
       {{6, 0}, {14, 0}, {13, -1}, {14, -1}},
       {{0, 1, T::DOUBLE}, {0, 2, T::DOUBLE}, {0, 3, T::SINGLE}},
       0},
      {"four-matches",
       {{6, 0}, {14, -1}, {14, -1}, {14, -1}, {14, -1}},
       {{0, 1, T::DOUBLE}, {2, 0, T::DOUBLE}, {0, 3, T::DOUBLE},
        {4, 0, T::DOUBLE}},
       0},
  };
  for (const auto &test_case : cases) {
    print_neighboring_si_case(test_case);
  }
  return 0;
}

static void print_valence4n_cleanup1_case(
    const Valence4NCleanup1Case &test_case) {
  RDKit::RWMol input(test_case.atoms, test_case.bonds);
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);

  bool returned = false;
  bool result = false;
  bool threw = false;
  unsigned int error_bond_index = 0;
  int error_bond_type = 0;
  try {
    result = RDKit::_Valence4NCleanUp1(
        molecule, molecule.getAtomWithIdx(test_case.atom_index));
    returned = true;
  } catch (const RDKit::ValueErrorException &error) {
    threw = true;
    error_bond_index = error.bond_index;
    error_bond_type = error.bond_type;
  }
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"11132300ed0a447a05c38e8eed0e710762514b9180aad984f12acc697255d87d\","
               "\"operation\":\"_Valence4NCleanUp1\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(input);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(input);
  std::cout << ",\"atom_index\":" << test_case.atom_index
            << "},\"output\":{\"status\":\""
            << (returned ? "return" : "exception") << "\",\"result\":";
  if (returned) {
    std::cout << (result ? "true" : "false");
  } else {
    std::cout << "null";
  }
  std::cout << ",\"exception\":";
  if (threw) {
    std::cout << "{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << '}';
  } else {
    std::cout << "null";
  }
  std::cout << ",\"diagnostics\":";
  if (threw) {
    std::cout << "[{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "}]";
  } else {
    std::cout << "[]";
  }
  std::cout << ",\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[]}}\n";
}

static int valence4n_cleanup1_records() {
  using T = RDKit::Bond::BondType;
  const std::vector<Valence4NCleanup1Case> cases = {
      {"wrong-element-short-circuit", {{6, -1, 0}, {6, 0, 0}}, {{0, 1, T::THREECENTER}}, 0},
      {"wrong-charge-short-circuit", {{7, 0, 0}, {6, 0, 0}}, {{0, 1, T::THREECENTER}}, 0},
      {"wrong-valence", {{7, -1, 0}, {6, 0, 0}}, {{0, 1, T::TRIPLE}}, 0},
      {"unsupported-valence-bond-exception", {{7, -1, 0}, {6, 0, 0}}, {{0, 1, T::THREECENTER}}, 0},
      {"explicit-h-valence-four-no-match", {{7, -1, 1}, {6, 0, 0}}, {{0, 1, T::TRIPLE}}, 0},
      {"zero-substructure-matches", {{7, -1, 0}, {7, 0, 0}, {7, 0, 0}}, {{0, 1, T::DOUBLE}, {2, 0, T::DOUBLE}}, 0},
      {"unique-match-reversed-endpoints-extra-edge-and-charges",
       {{6, 4, 0}, {7, 2, 0}, {7, -1, 0}, {7, -3, 0}, {7, 1, 0}, {8, 0, 0}},
       {{1, 0, T::SINGLE, RDKit::Bond::BEGINDASH}, {2, 1, T::DOUBLE}, {3, 2, T::DOUBLE},
        {4, 3, T::SINGLE}, {0, 4, T::DOUBLE}, {0, 5, T::SINGLE}},
       2},
      {"multiple-unique-atom-set-matches",
       {{6, 0, 0}, {7, 0, 0}, {7, -1, 0}, {7, 0, 0}, {7, 0, 0},
        {6, 0, 0}, {7, 0, 0}},
       {{0, 1, T::SINGLE}, {1, 2, T::DOUBLE}, {2, 3, T::DOUBLE},
        {3, 4, T::SINGLE}, {4, 0, T::DOUBLE}, {1, 5, T::SINGLE},
        {5, 6, T::DOUBLE}, {6, 3, T::SINGLE}},
       2},
      {"unique-match-does-not-contain-argument-atom",
       {{7, -1, 0}, {8, 0, 0}, {8, 0, 0}, {6, 0, 0}, {7, 0, 0},
        {50, 0, 0}, {7, 0, 0}, {7, 0, 0}},
       {{0, 1, T::DOUBLE}, {2, 0, T::DOUBLE}, {3, 4, T::SINGLE},
        {4, 5, T::DOUBLE}, {5, 6, T::DOUBLE}, {6, 7, T::SINGLE},
        {7, 3, T::DOUBLE}},
       0},
  };
  for (const auto &test_case : cases) {
    print_valence4n_cleanup1_case(test_case);
  }
  return 0;
}

static void print_valence4n_cleanup2_case(
    const Valence4NCleanup2Case &test_case) {
  RDKit::RWMol input(test_case.atoms, test_case.bonds);
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);
  const bool result = RDKit::_Valence4NCleanUp2(
      molecule, molecule.getAtomWithIdx(test_case.atom_index));
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"c2fe8b4125473236e4fa1c638d93c23365c1c7df416f17537259dcc24c4d7587\","
               "\"operation\":\"_Valence4NCleanUp2\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(input);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(input);
  std::cout << ",\"atom_index\":" << test_case.atom_index
            << "},\"output\":{\"status\":\"return\",\"result\":"
            << (result ? "true" : "false")
            << ",\"exception\":null,\"diagnostics\":[],\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[]}}\n";
}

static int valence4n_cleanup2_records() {
  using D = RDKit::Bond::BondDir;
  using T = RDKit::Bond::BondType;
  const std::vector<Valence4NCleanup2Case> cases = {
      {"isolated-no-target", {{7, -1, 2}}, {}, 0},
      {"all-target-predicates-miss",
       {{7, -1, 1}, {6, 0, 2}, {7, 1, 3}, {7, 0, 4}, {7, 0, 5}},
       {{0, 1, T::DOUBLE, D::BEGINDASH},
        {0, 2, T::DOUBLE, D::ENDUPRIGHT},
        {0, 3, T::SINGLE, D::BEGINWEDGE},
        {0, 4, T::TRIPLE, D::ENDDOWNRIGHT}},
       0},
      {"unique-target-reversed-endpoints-preserves-fields",
       {{6, 5, 2}, {7, 0, 1}, {8, -2, 3}},
       {{1, 0, T::DOUBLE, D::BEGINWEDGE},
        {0, 2, T::SINGLE, D::ENDDOWNRIGHT}},
       0},
      {"multiple-targets-first-insertion",
       {{7, -1, 0}, {7, 0, 1}, {7, 0, 2}},
       {{0, 2, T::DOUBLE, D::BEGINDASH},
        {0, 1, T::DOUBLE, D::ENDUPRIGHT}},
       0},
      {"multiple-targets-reversed-insertion",
       {{7, -1, 0}, {7, 0, 1}, {7, 0, 2}},
       {{0, 1, T::DOUBLE, D::ENDUPRIGHT},
        {0, 2, T::DOUBLE, D::BEGINDASH}},
       0},
      {"miss-before-later-target",
       {{8, -4, 6}, {6, 0, 0}, {7, 0, 7}},
       {{0, 1, T::DOUBLE, D::EITHERDOUBLE},
        {2, 0, T::DOUBLE, D::UNKNOWN}},
       0},
  };
  for (const auto &test_case : cases) {
    print_valence4n_cleanup2_case(test_case);
  }
  return 0;
}

static void print_valence5n_cleanup1_case(
    const Valence5NCleanup1Case &test_case) {
  RDKit::RWMol input(test_case.atoms, test_case.bonds);
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);
  bool result = false;
  bool threw = false;
  unsigned int error_bond_index = 0;
  int error_bond_type = 0;
  try {
    result = RDKit::_Valence5NCleanUp1(
        molecule, molecule.getAtomWithIdx(test_case.atom_index));
  } catch (const RDKit::ValueErrorException &error) {
    threw = true;
    error_bond_index = error.bond_index;
    error_bond_type = error.bond_type;
  }
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"58b4273fc255b5ed398a0d81b30069e58d68478d083b499f41a77ec508250f06\","
               "\"operation\":\"_Valence5NCleanUp1\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(input);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(input);
  std::cout << ",\"atom_index\":" << test_case.atom_index
            << "},\"output\":{\"status\":\""
            << (threw ? "exception" : "return") << "\",\"result\":";
  if (threw) {
    std::cout << "null,\"exception\":{\"kind\":\"ValueErrorException\","
                 "\"message\":\"Bad bond type\",\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "},\"diagnostics\":[{\"kind\":\"ValueErrorException\","
                 "\"message\":\"Bad bond type\",\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "}]";
  } else {
    std::cout << (result ? "true" : "false")
              << ",\"exception\":null,\"diagnostics\":[]";
  }
  std::cout << ",\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[]}}\n";
}

static int valence5n_cleanup1_records() {
  using D = RDKit::Bond::BondDir;
  using T = RDKit::Bond::BondType;
  const std::vector<Valence5NCleanup1Case> cases = {
      {"no-target", {{7, 0, 0}, {7, 0, 0}, {7, 1, 0}},
       {{0, 1, T::DOUBLE}, {0, 2, T::SINGLE}}, 0},
      {"direct-reversed-preserves-fields", {{6, -3, 2}, {7, 1, 1}, {8, 2, 0}},
       {{1, 0, T::DOUBLE, D::BEGINDASH}, {1, 2, T::SINGLE, D::ENDUPRIGHT}}, 0},
      {"alternating-path", {{7, -1, 0}, {6, 0, 0}, {6, 0, 0}, {7, 1, 0}},
       {{0, 1, T::DOUBLE}, {1, 2, T::SINGLE}, {2, 3, T::DOUBLE}}, 0},
      {"over-depth-limit", {{7, -1, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0},
                            {6, 0, 0}, {6, 0, 0}, {7, 1, 0}},
       {{0, 1, T::DOUBLE}, {1, 2, T::SINGLE}, {2, 3, T::DOUBLE},
        {3, 4, T::SINGLE}, {4, 5, T::DOUBLE}, {5, 6, T::SINGLE}}, 0},
      {"valence-exception-partial-mutation", {{7, -1, 0}, {7, 1, 0}, {6, 0, 0}},
       {{0, 1, T::DOUBLE}, {1, 2, T::THREECENTER}}, 0},
  };
  for (const auto &test_case : cases) print_valence5n_cleanup1_case(test_case);
  return 0;
}

static void print_valence5n_cleanup2_case(
    const Valence5NCleanup2Case &test_case) {
  RDKit::RWMol input(test_case.atoms, test_case.bonds);
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);
  bool result = false;
  bool threw = false;
  unsigned int error_bond_index = 0;
  int error_bond_type = 0;
  try {
    result = RDKit::_Valence5NCleanUp2(
        molecule, molecule.getAtomWithIdx(test_case.atom_index));
  } catch (const RDKit::ValueErrorException &error) {
    threw = true;
    error_bond_index = error.bond_index;
    error_bond_type = error.bond_type;
  }
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"07654fddc42017f1aa76c39ec11013301f92e460b199218c6a2f3bf6bcfb4484\","
               "\"operation\":\"_Valence5NCleanUp2\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(input);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(input);
  std::cout << ",\"atom_index\":" << test_case.atom_index
            << "},\"output\":{\"status\":\""
            << (threw ? "exception" : "return") << "\",\"result\":";
  if (threw) {
    std::cout << "null,\"exception\":{\"kind\":\"ValueErrorException\","
                 "\"message\":\"Bad bond type\",\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "},\"diagnostics\":[{\"kind\":\"ValueErrorException\","
                 "\"message\":\"Bad bond type\",\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "}]";
  } else {
    std::cout << (result ? "true" : "false")
              << ",\"exception\":null,\"diagnostics\":[]";
  }
  std::cout << ",\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[]}}\n";
}

static int valence5n_cleanup2_records() {
  using D = RDKit::Bond::BondDir;
  using T = RDKit::Bond::BondType;
  const std::vector<Valence5NCleanup2Case> cases = {
      {"no-target",
       {{7, 0, 1}, {6, 4, 2}, {7, -1, 3}, {7, 0, 4}},
       {{0, 1, T::DOUBLE, D::BEGINDASH},
        {1, 2, T::SINGLE, D::ENDUPRIGHT},
        {0, 3, T::TRIPLE, D::UNKNOWN}},
       0},
      {"over-depth-limit",
       {{7, 0, 0}, {6, 0, 0}, {6, 0, 0}, {7, -1, 0}},
       {{0, 1, T::TRIPLE}, {1, 2, T::SINGLE}, {2, 3, T::DOUBLE}},
       0},
      {"begin-is-root-preserves-fields",
       {{7, 3, 1}, {6, 4, 2}, {7, -1, 3}, {8, 2, 4}},
       {{0, 1, T::TRIPLE, D::BEGINWEDGE},
        {1, 2, T::SINGLE, D::ENDDOWNRIGHT},
        {1, 3, T::SINGLE, D::BEGINDASH}},
       0},
      {"multiple-targets-first-adjacency-reversed-root",
       {{7, 0, 0}, {6, 0, 0}, {7, -1, 0}, {7, -1, 0}},
       {{1, 0, T::TRIPLE, D::ENDUPRIGHT},
        {2, 1, T::SINGLE, D::BEGINWEDGE},
        {1, 3, T::SINGLE, D::BEGINDASH}},
       0},
      {"target-valence-exception-partial-mutation",
       {{7, 0, 0}, {6, 0, 0}, {7, -1, 0}, {6, 0, 0}},
       {{0, 1, T::TRIPLE}, {1, 2, T::SINGLE}, {2, 3, T::THREECENTER}},
       0},
      {"root-valence-exception-after-target-check",
       {{7, 0, 0}, {6, 0, 0}, {7, -1, 0}, {6, 0, 0}},
       {{0, 1, T::TRIPLE}, {1, 2, T::SINGLE}, {0, 3, T::THREECENTER}},
       0},
  };
  for (const auto &test_case : cases) print_valence5n_cleanup2_case(test_case);
  return 0;
}

static void print_valence5n_cleanup3_case(
    const Valence5NCleanup3Case &test_case) {
  RDKit::RWMol input(test_case.atoms, test_case.bonds);
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);
  bool result = false;
  bool threw = false;
  unsigned int error_bond_index = 0;
  int error_bond_type = 0;
  try {
    result = RDKit::_Valence5NCleanUp3(
        molecule, molecule.getAtomWithIdx(test_case.atom_index));
  } catch (const RDKit::ValueErrorException &error) {
    threw = true;
    error_bond_index = error.bond_index;
    error_bond_type = error.bond_type;
  }
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"bf11af3caae7ff7d1dd7b381defd449e6cb9b193e9dab328a35a6ddd5d73e25d\","
               "\"operation\":\"_Valence5NCleanUp3\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(input);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(input);
  std::cout << ",\"atom_index\":" << test_case.atom_index
            << "},\"output\":{\"status\":\""
            << (threw ? "exception" : "return") << "\",\"result\":";
  if (threw) {
    std::cout << "null,\"exception\":{\"kind\":\"ValueErrorException\","
                 "\"message\":\"Bad bond type\",\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "},\"diagnostics\":[{\"kind\":\"ValueErrorException\","
                 "\"message\":\"Bad bond type\",\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "}]";
  } else {
    std::cout << (result ? "true" : "false")
              << ",\"exception\":null,\"diagnostics\":[]";
  }
  std::cout << ",\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[]}}\n";
}

static int valence5n_cleanup3_records() {
  using D = RDKit::Bond::BondDir;
  using T = RDKit::Bond::BondType;
  const std::vector<Valence5NCleanup3Case> cases = {
      {"no-neutral-nitrogen", {{7, 0, 0}, {7, 1, 0}, {8, 0, 0}},
       {{0, 1, T::DOUBLE}, {0, 2, T::DOUBLE}}, 0},
      {"oxygen-guard-true-no-mutation",
       {{7, -3, 1}, {7, 0, 2}, {8, 0, 3}, {6, 2, 4}},
       {{1, 0, T::DOUBLE, D::BEGINWEDGE},
        {0, 2, T::DOUBLE, D::ENDDOWNRIGHT},
        {0, 3, T::TRIPLE, D::BEGINDASH}},
       0},
      {"charged-oxygen-does-not-guard",
       {{7, -3, 1}, {7, 0, 2}, {8, 1, 3}, {6, 2, 0}},
       {{1, 0, T::DOUBLE, D::ENDUPRIGHT},
        {0, 2, T::DOUBLE, D::ENDDOWNRIGHT},
        {0, 3, T::SINGLE}},
       0},
      {"target-valence-exception-after-charge",
       {{7, 0, 0}, {7, 0, 0}, {6, 0, 0}},
       {{0, 1, T::DOUBLE}, {1, 2, T::THREECENTER}}, 0},
      {"root-valence-exception-after-bond-and-charges",
       {{7, 0, 0}, {7, 0, 0}, {6, 0, 0}},
       {{0, 1, T::DOUBLE, D::BEGINDASH}, {0, 2, T::THREECENTER}}, 0},
  };
  for (const auto &test_case : cases) print_valence5n_cleanup3_case(test_case);
  return 0;
}

static void print_valence5n_cleanup4_case(
    const Valence5NCleanup4Case &test_case) {
  RDKit::RWMol input(test_case.atoms, test_case.bonds);
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);
  const bool result = RDKit::_Valence5NCleanUp4(
      molecule, molecule.getAtomWithIdx(test_case.atom_index));
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"701347673195ee11dec11ba3013f6eb10746d62e0a879ce63d0de98ec0642746\","
               "\"operation\":\"_Valence5NCleanUp4\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(input);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(input);
  std::cout << ",\"atom_index\":" << test_case.atom_index
            << "},\"output\":{\"status\":\"return\",\"result\":"
            << (result ? "true" : "false")
            << ",\"exception\":null,\"diagnostics\":[],\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[]}}\n";
}

static int valence5n_cleanup4_records() {
  using D = RDKit::Bond::BondDir;
  using T = RDKit::Bond::BondType;
  const std::vector<Valence5NCleanup4Case> cases = {
      {"no-match-all-predicates",
       {{7, 3, 1}, {14, -1, 2}, {14, 0, 3}, {6, -1, 4}},
       {{0, 1, T::SINGLE, D::BEGINWEDGE},
        {0, 2, T::DOUBLE, D::ENDUPRIGHT},
        {0, 3, T::DOUBLE, D::BEGINDASH}},
       0},
      {"one-match",
       {{7, 3, 1}, {14, -1, 2}, {14, -1, 3}, {14, 0, 4}},
       {{0, 1, T::DOUBLE, D::ENDDOWNRIGHT},
        {0, 2, T::SINGLE, D::BEGINWEDGE},
        {0, 3, T::DOUBLE, D::BEGINDASH}},
       0},
      {"exactly-two-reversed-preserves-fields",
       {{7, 3, 1}, {6, 0, 0}, {14, -1, 2}, {14, 0, 3},
        {14, -1, 4}, {16, 1, 2}, {6, 0, 3}},
       {{0, 1, T::DOUBLE, D::BEGINWEDGE},
        {2, 0, T::DOUBLE, D::ENDUPRIGHT},
        {0, 3, T::DOUBLE, D::BEGINDASH},
        {4, 0, T::DOUBLE, D::ENDDOWNRIGHT},
        {5, 6, T::DOUBLE, D::UNKNOWN}},
       0},
      {"three-matches-early-return",
       {{7, 3, 1}, {14, -1, 2}, {14, -1, 3}, {14, -1, 4}},
       {{1, 0, T::DOUBLE, D::BEGINWEDGE},
        {0, 2, T::DOUBLE, D::ENDUPRIGHT},
        {3, 0, T::DOUBLE, D::BEGINDASH}},
       0},
      {"exactly-two-with-late-nonmatch",
       {{7, -4, 1}, {14, -1, 2}, {14, -1, 3}, {14, -1, 4}},
       {{0, 1, T::DOUBLE, D::BEGINDASH},
        {0, 2, T::DOUBLE, D::ENDUPRIGHT},
        {0, 3, T::TRIPLE, D::ENDDOWNRIGHT}},
       0},
  };
  for (const auto &test_case : cases) print_valence5n_cleanup4_case(test_case);
  return 0;
}

static void print_valence5n_cleanup5_case(
    const Valence5NCleanup5Case &test_case) {
  RDKit::RWMol input(test_case.atoms, test_case.bonds);
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);
  bool result = false;
  std::string exception_kind;
  std::string exception_expression;
  std::string exception_message;
  unsigned int error_bond_index = 0;
  int error_bond_type = 0;
  try {
    result = RDKit::_Valence5NCleanUp5(
        molecule, molecule.getAtomWithIdx(test_case.atom_index),
        test_case.atomic_number);
  } catch (const RDKit::InvariantException &error) {
    exception_kind = error.kind;
    exception_expression = error.expression;
    exception_message = error.message;
  } catch (const RDKit::ValueErrorException &error) {
    exception_kind = "ValueErrorException";
    exception_message = error.what();
    error_bond_index = error.bond_index;
    error_bond_type = error.bond_type;
  }
  const bool threw = !exception_kind.empty();
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"572b129f3625cbffe00811214403e660abdcf283c582ff90287011755698b73e\","
               "\"operation\":\"_Valence5NCleanUp5\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(input);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(input);
  std::cout << ",\"atom_index\":" << test_case.atom_index
            << ",\"atomic_number\":" << test_case.atomic_number
            << "},\"output\":{\"status\":\""
            << (threw ? "exception" : "return") << "\",\"result\":";
  if (!threw) {
    std::cout << (result ? "true" : "false")
              << ",\"exception\":null,\"diagnostics\":[]";
  } else if (exception_kind == "ValueErrorException") {
    std::cout << "null,\"exception\":{\"kind\":\"ValueErrorException\","
                 "\"message\":\"Bad bond type\",\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "},\"diagnostics\":[{\"kind\":\"ValueErrorException\","
                 "\"message\":\"Bad bond type\",\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "}]";
  } else {
    std::cout << "null,\"exception\":{\"kind\":\"" << exception_kind
              << "\",\"expression\":\"" << exception_expression
              << "\",\"message\":\"" << exception_message
              << "\"},\"diagnostics\":[{\"kind\":\"" << exception_kind
              << "\",\"expression\":\"" << exception_expression
              << "\",\"message\":\"" << exception_message << "\"}]";
  }
  std::cout << ",\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[]}}\n";
}

static int valence5n_cleanup5_records() {
  using D = RDKit::Bond::BondDir;
  using T = RDKit::Bond::BondType;
  const std::vector<Valence5NCleanup5Case> cases = {
      {"precondition-invalid", {{7, 0, 0}, {8, 0, 0}},
       {{0, 1, T::DOUBLE}}, 0, 7},
      {"no-target-oxygen", {{7, -2, 1}, {6, 1, 2}},
       {{0, 1, T::DOUBLE, D::BEGINWEDGE}}, 0, 8},
      {"no-target-sulfur", {{7, -2, 1}, {6, 1, 2}},
       {{0, 1, T::DOUBLE, D::BEGINDASH}}, 0, 16},
      {"no-target-fluorine", {{7, -2, 1}, {6, 1, 2}},
       {{0, 1, T::DOUBLE, D::ENDUPRIGHT}}, 0, 9},
      {"no-target-chlorine", {{7, -2, 1}, {6, 1, 2}},
       {{0, 1, T::DOUBLE, D::ENDDOWNRIGHT}}, 0, 17},
      {"only-uncharged-alternating",
       {{7, -3, 1}, {6, 2, 0}, {6, -2, 0}, {8, 0, 4}},
       {{0, 1, T::DOUBLE, D::BEGINDASH},
        {1, 2, T::SINGLE, D::ENDUPRIGHT},
        {2, 3, T::DOUBLE, D::BEGINWEDGE}},
       0, 8},
      {"only-charged-reversed",
       {{7, -3, 0}, {16, 1, 3}},
       {{1, 0, T::DOUBLE, D::ENDDOWNRIGHT}}, 0, 16},
      {"both-targets-neutral-path-only",
       {{7, -3, 0}, {9, 0, 5}, {9, 1, 2}, {6, 4, 0}},
       {{0, 1, T::DOUBLE, D::BEGINWEDGE},
        {2, 0, T::DOUBLE, D::ENDUPRIGHT},
        {0, 3, T::SINGLE, D::BEGINDASH}},
       0, 9},
      {"depth-seven-inclusive",
       {{7, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0},
        {6, 0, 0}, {6, 0, 0}, {6, 0, 0}, {17, 0, 0}},
       {{0, 1, T::DOUBLE}, {1, 2, T::SINGLE}, {2, 3, T::DOUBLE},
        {3, 4, T::SINGLE}, {4, 5, T::DOUBLE}, {5, 6, T::SINGLE},
        {6, 7, T::DOUBLE}},
       0, 17},
      {"over-depth-limit",
       {{7, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0},
        {6, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0}, {8, 0, 0}},
       {{0, 1, T::DOUBLE}, {1, 2, T::SINGLE}, {2, 3, T::DOUBLE},
        {3, 4, T::SINGLE}, {4, 5, T::DOUBLE}, {5, 6, T::SINGLE},
        {6, 7, T::DOUBLE}, {7, 8, T::SINGLE}, {8, 9, T::DOUBLE}},
       0, 8},
      {"only-uncharged-valence-exception",
       {{7, 0, 0}, {8, 0, 0}, {6, 0, 0}},
       {{0, 1, T::DOUBLE}, {1, 2, T::THREECENTER}}, 0, 8},
      {"only-charged-valence-exception",
       {{7, 0, 0}, {8, 1, 2}, {6, 0, 0}},
       {{0, 1, T::DOUBLE}, {1, 2, T::THREECENTER}}, 0, 8},
      {"both-charged-first-exception",
       {{7, 0, 0}, {8, 0, 3}, {8, 1, 4}, {6, 0, 0}, {6, 0, 0}},
       {{0, 1, T::DOUBLE}, {0, 2, T::DOUBLE},
        {1, 3, T::THREECENTER}, {2, 4, T::THREECENTER}},
       0, 8},
      {"both-uncharged-second-exception",
       {{7, 0, 0}, {8, 0, 0}, {8, 1, 0}, {6, 0, 0}},
       {{0, 1, T::DOUBLE}, {0, 2, T::DOUBLE},
        {1, 3, T::THREECENTER}},
       0, 8},
  };
  for (const auto &test_case : cases) print_valence5n_cleanup5_case(test_case);
  return 0;
}

static Valence5NCleanup6Case make_valence5n_cleanup6_case(
    std::string case_id, int root_atomic_number,
    RDKit::Bond::BondType bridge_type) {
  using D = RDKit::Bond::BondDir;
  using T = RDKit::Bond::BondType;
  return {std::move(case_id),
          {{6, -2, 1}, {6, 3, 0}, {root_atomic_number, 0, 0},
           {6, -4, 0}, {6, 5, 2}, {7, -6, 3}, {6, 7, 0}},
          {{1, 0, T::SINGLE, D::BEGINWEDGE},
           {2, 1, T::DOUBLE, D::BEGINDASH},
           {2, 3, T::DOUBLE, D::ENDDOWNRIGHT},
           {4, 3, bridge_type, D::ENDUPRIGHT},
           {4, 5, T::SINGLE, D::EITHERDOUBLE},
           {0, 5, T::DOUBLE, D::UNKNOWN},
           {6, 2, T::SINGLE, D::BEGINWEDGE}},
          2};
}

static void append_valence5n_cleanup6_component(
    Valence5NCleanup6Case &destination,
    const Valence5NCleanup6Case &component) {
  const auto offset = static_cast<unsigned int>(destination.atoms.size());
  destination.atoms.insert(destination.atoms.end(), component.atoms.begin(),
                           component.atoms.end());
  for (auto bond : component.bonds) {
    bond.begin_atom_index += offset;
    bond.end_atom_index += offset;
    destination.bonds.push_back(bond);
  }
}

static void print_valence5n_cleanup6_case(
    const Valence5NCleanup6Case &test_case) {
  RDKit::RWMol input(test_case.atoms, test_case.bonds);
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);

  bool returned = false;
  bool result = false;
  bool threw = false;
  unsigned int error_bond_index = 0;
  int error_bond_type = 0;
  try {
    result = RDKit::_Valence5NCleanUp6(
        molecule, molecule.getAtomWithIdx(test_case.atom_index));
    returned = true;
  } catch (const RDKit::ValueErrorException &error) {
    threw = true;
    error_bond_index = error.bond_index;
    error_bond_type = error.bond_type;
  }
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"33f731b27422f52bbe7b9304df66f3126a325a3458a577bcf9b29c01d4b9f47d\","
               "\"operation\":\"_Valence5NCleanUp6\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(input);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(input);
  std::cout << ",\"atom_index\":" << test_case.atom_index
            << "},\"output\":{\"status\":\""
            << (returned ? "return" : "exception") << "\",\"result\":";
  if (returned) {
    std::cout << (result ? "true" : "false");
  } else {
    std::cout << "null";
  }
  std::cout << ",\"exception\":";
  if (threw) {
    std::cout << "{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << '}';
  } else {
    std::cout << "null";
  }
  std::cout << ",\"diagnostics\":";
  if (threw) {
    std::cout << "[{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "}]";
  } else {
    std::cout << "[]";
  }
  std::cout << ",\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[]}}\n";
}

static int valence5n_cleanup6_records() {
  using T = RDKit::Bond::BondType;
  std::vector<Valence5NCleanup6Case> cases = {
      {"wrong-element-short-circuit", {{6, 0, 0}, {6, 0, 0}},
       {{0, 1, T::THREECENTER}}, 0},
      {"wrong-charge-short-circuit", {{7, 1, 0}, {6, 0, 0}},
       {{0, 1, T::THREECENTER}}, 0},
      {"wrong-valence",
       {{7, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0}},
       {{0, 1, T::SINGLE}, {0, 2, T::SINGLE}, {0, 3, T::SINGLE},
        {0, 4, T::SINGLE}},
       0},
      {"unsupported-valence-bond-exception", {{7, 0, 0}, {6, 0, 0}},
       {{0, 1, T::THREECENTER}}, 0},
      {"zero-substructure-matches",
       {{7, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0},
        {6, 0, 0}},
       {{0, 1, T::SINGLE}, {2, 0, T::SINGLE}, {0, 3, T::SINGLE},
        {4, 0, T::SINGLE}, {0, 5, T::SINGLE}},
       0},
  };

  for (unsigned int atom_index : {0U, 1U, 3U, 4U, 5U, 6U}) {
    auto test_case = make_valence5n_cleanup6_case(
        "atom-predicate-miss-" + std::to_string(atom_index), 7, T::TRIPLE);
    test_case.atoms[atom_index].atomic_number = 8;
    cases.push_back(std::move(test_case));
  }
  for (unsigned int bond_index : {0U, 4U, 5U}) {
    auto test_case = make_valence5n_cleanup6_case(
        "nonroot-bond-predicate-miss-" + std::to_string(bond_index), 7,
        T::TRIPLE);
    test_case.bonds[bond_index].bond_type =
        bond_index == 5 ? T::SINGLE : T::DOUBLE;
    cases.push_back(std::move(test_case));
  }
  auto root_bonds_one =
      make_valence5n_cleanup6_case("root-bonds-triple-single-single", 7,
                                   T::TRIPLE);
  root_bonds_one.bonds[1].bond_type = T::TRIPLE;
  root_bonds_one.bonds[2].bond_type = T::SINGLE;
  cases.push_back(std::move(root_bonds_one));
  auto root_bonds_two = make_valence5n_cleanup6_case(
      "root-bonds-double-single-double", 7, T::TRIPLE);
  root_bonds_two.bonds[2].bond_type = T::SINGLE;
  root_bonds_two.bonds[6].bond_type = T::DOUBLE;
  cases.push_back(std::move(root_bonds_two));

  for (int bond_type = static_cast<int>(T::UNSPECIFIED);
       bond_type <= static_cast<int>(T::ZERO); ++bond_type) {
    cases.push_back(make_valence5n_cleanup6_case(
        "unspecified-bridge-type-" + std::to_string(bond_type), 7,
        static_cast<T>(bond_type)));
  }

  auto multiple =
      make_valence5n_cleanup6_case("multiple-unique-atom-set-matches", 7,
                                   T::TRIPLE);
  const auto second =
      make_valence5n_cleanup6_case("unused", 50, T::AROMATIC);
  append_valence5n_cleanup6_component(multiple, second);
  cases.push_back(std::move(multiple));

  Valence5NCleanup6Case unrelated = {
      "unique-match-does-not-contain-argument-atom",
      {{7, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0},
       {6, 0, 0}},
      {{0, 1, T::SINGLE}, {2, 0, T::SINGLE}, {0, 3, T::SINGLE},
       {4, 0, T::SINGLE}, {0, 5, T::SINGLE}},
      0};
  append_valence5n_cleanup6_component(unrelated, second);
  cases.push_back(std::move(unrelated));

  for (const auto &test_case : cases) {
    print_valence5n_cleanup6_case(test_case);
  }
  return 0;
}

static Valence5NCleanup7Case make_valence5n_cleanup7_case(
    std::string case_id, int root_atomic_number,
    RDKit::Bond::BondType bridge_type) {
  using D = RDKit::Bond::BondDir;
  using T = RDKit::Bond::BondType;
  return {std::move(case_id),
          {{6, -2, 1}, {6, 3, 0}, {root_atomic_number, 0, 0},
           {7, -4, 0}, {6, 5, 0}, {8, -6, 2}, {6, 7, 0}, {8, 0, 0}},
          {{1, 0, bridge_type, D::BEGINWEDGE},
           {2, 1, T::DOUBLE, D::BEGINDASH},
           {2, 3, T::DOUBLE, D::ENDDOWNRIGHT},
           {3, 4, T::SINGLE, D::ENDUPRIGHT},
           {4, 5, T::SINGLE, D::EITHERDOUBLE},
           {5, 0, T::SINGLE, D::UNKNOWN},
           {6, 2, T::SINGLE, D::BEGINWEDGE},
           {4, 7, T::DOUBLE, D::BEGINDASH}},
          2};
}

static void append_valence5n_cleanup7_component(
    Valence5NCleanup7Case &destination,
    const Valence5NCleanup7Case &component) {
  const auto offset = static_cast<unsigned int>(destination.atoms.size());
  destination.atoms.insert(destination.atoms.end(), component.atoms.begin(),
                           component.atoms.end());
  for (auto bond : component.bonds) {
    bond.begin_atom_index += offset;
    bond.end_atom_index += offset;
    destination.bonds.push_back(bond);
  }
}

static void print_valence5n_cleanup7_case(
    const Valence5NCleanup7Case &test_case) {
  RDKit::RWMol input(test_case.atoms, test_case.bonds);
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);

  bool returned = false;
  bool result = false;
  bool threw = false;
  unsigned int error_bond_index = 0;
  int error_bond_type = 0;
  try {
    result = RDKit::_Valence5NCleanUp7(
        molecule, molecule.getAtomWithIdx(test_case.atom_index));
    returned = true;
  } catch (const RDKit::ValueErrorException &error) {
    threw = true;
    error_bond_index = error.bond_index;
    error_bond_type = error.bond_type;
  }
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"52035cacfe7ced70b44df5ac007a41f101a968fbee2a5338abf5ff81a4efa562\","
               "\"operation\":\"_Valence5NCleanUp7\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(input);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(input);
  std::cout << ",\"atom_index\":" << test_case.atom_index
            << "},\"output\":{\"status\":\""
            << (returned ? "return" : "exception") << "\",\"result\":";
  if (returned) {
    std::cout << (result ? "true" : "false");
  } else {
    std::cout << "null";
  }
  std::cout << ",\"exception\":";
  if (threw) {
    std::cout << "{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << '}';
  } else {
    std::cout << "null";
  }
  std::cout << ",\"diagnostics\":";
  if (threw) {
    std::cout << "[{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "}]";
  } else {
    std::cout << "[]";
  }
  std::cout << ",\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[]}}\n";
}

static int valence5n_cleanup7_records() {
  using T = RDKit::Bond::BondType;
  std::vector<Valence5NCleanup7Case> cases = {
      {"no-target-before-preconditions", {{7, 0, 0}, {6, 0, 0}},
       {{0, 1, T::THREECENTER}}, 0},
  };

  cases.push_back(
      make_valence5n_cleanup7_case("wrong-element-after-search", 6, T::TRIPLE));
  auto wrong_charge =
      make_valence5n_cleanup7_case("wrong-charge-after-search", 7, T::TRIPLE);
  wrong_charge.atoms[2].formal_charge = 1;
  cases.push_back(std::move(wrong_charge));
  auto wrong_valence =
      make_valence5n_cleanup7_case("wrong-valence-after-search", 7, T::TRIPLE);
  wrong_valence.atoms[2].num_explicit_hydrogens = 1;
  cases.push_back(std::move(wrong_valence));
  cases.push_back({"unsupported-valence-bond-exception-after-search",
                   {{7, 0, 0}, {6, 0, 0}, {6, 0, 0}, {8, 0, 0},
                    {6, 0, 0}},
                   {{0, 1, T::DOUBLE}, {1, 2, T::SINGLE},
                    {2, 3, T::DOUBLE}, {0, 4, T::THREECENTER}},
                   0});

  auto zero_matches =
      make_valence5n_cleanup7_case("zero-substructure-matches", 7, T::TRIPLE);
  zero_matches.atoms[5].atomic_number = 9;
  cases.push_back(std::move(zero_matches));

  for (unsigned int atom_index : {0U, 1U, 3U, 4U, 5U, 6U}) {
    auto test_case = make_valence5n_cleanup7_case(
        "atom-predicate-miss-" + std::to_string(atom_index), 7, T::TRIPLE);
    test_case.atoms[atom_index].atomic_number = 9;
    cases.push_back(std::move(test_case));
  }
  for (unsigned int bond_index : {4U, 5U}) {
    auto test_case = make_valence5n_cleanup7_case(
        "nonpath-bond-predicate-miss-" + std::to_string(bond_index), 7,
        T::TRIPLE);
    test_case.bonds[bond_index].bond_type = T::DOUBLE;
    cases.push_back(std::move(test_case));
  }

  for (int bond_type = static_cast<int>(T::UNSPECIFIED);
       bond_type <= static_cast<int>(T::ZERO); ++bond_type) {
    cases.push_back(make_valence5n_cleanup7_case(
        "unspecified-bridge-type-" + std::to_string(bond_type), 7,
        static_cast<T>(bond_type)));
  }

  auto depth_limit = make_valence5n_cleanup7_case(
      "depth-five-nontarget-stops-search", 7, T::TRIPLE);
  depth_limit.bonds.pop_back();
  depth_limit.atoms.insert(depth_limit.atoms.end(),
                           {{6, 0, 0}, {6, 0, 0}, {6, 0, 0}});
  depth_limit.bonds.insert(
      depth_limit.bonds.end(),
      {{4, 8, T::DOUBLE}, {8, 9, T::SINGLE}, {9, 10, T::DOUBLE}});
  cases.push_back(std::move(depth_limit));

  auto multiple = make_valence5n_cleanup7_case(
      "multiple-unique-atom-set-matches", 7, T::TRIPLE);
  const auto second =
      make_valence5n_cleanup7_case("unused", 50, T::AROMATIC);
  append_valence5n_cleanup7_component(multiple, second);
  cases.push_back(std::move(multiple));

  Valence5NCleanup7Case unrelated = {
      "unique-match-does-not-contain-argument-or-target",
      {{7, 0, 0}, {6, 0, 0}, {6, 0, 0}, {8, 0, 0}, {6, 0, 0},
       {6, 0, 0}},
      {{0, 1, T::DOUBLE}, {1, 2, T::SINGLE}, {2, 3, T::DOUBLE},
       {0, 4, T::DOUBLE}, {0, 5, T::SINGLE}},
      0};
  append_valence5n_cleanup7_component(unrelated, second);
  cases.push_back(std::move(unrelated));

  for (const auto &test_case : cases) {
    print_valence5n_cleanup7_case(test_case);
  }
  return 0;
}

static Valence5NCleanup8Case make_valence5n_cleanup8_case(
    std::string case_id, int root_atomic_number) {
  using D = RDKit::Bond::BondDir;
  using T = RDKit::Bond::BondType;
  return {std::move(case_id),
          {{6, -2, 0}, {7, 3, 0}, {6, -4, 0},
           {7, 5, 0},  {7, -6, 0}, {root_atomic_number, 0, 3}},
          {{0, 1, T::SINGLE, D::BEGINWEDGE},
           {1, 2, T::DOUBLE, D::BEGINDASH},
           {2, 3, T::SINGLE, D::ENDDOWNRIGHT},
           {3, 4, T::DOUBLE, D::ENDUPRIGHT},
           {4, 0, T::SINGLE, D::EITHERDOUBLE},
           {5, 0, T::DOUBLE, D::UNKNOWN}},
          5};
}

static void append_valence5n_cleanup8_component(
    Valence5NCleanup8Case &destination,
    const Valence5NCleanup8Case &component) {
  const auto offset = static_cast<unsigned int>(destination.atoms.size());
  destination.atoms.insert(destination.atoms.end(), component.atoms.begin(),
                           component.atoms.end());
  for (auto bond : component.bonds) {
    bond.begin_atom_index += offset;
    bond.end_atom_index += offset;
    destination.bonds.push_back(bond);
  }
}

static void print_valence5n_cleanup8_case(
    const Valence5NCleanup8Case &test_case) {
  RDKit::RWMol input(test_case.atoms, test_case.bonds);
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);

  bool returned = false;
  bool result = false;
  bool threw = false;
  unsigned int error_bond_index = 0;
  int error_bond_type = 0;
  try {
    result = RDKit::_Valence5NCleanUp8(
        molecule, molecule.getAtomWithIdx(test_case.atom_index));
    returned = true;
  } catch (const RDKit::ValueErrorException &error) {
    threw = true;
    error_bond_index = error.bond_index;
    error_bond_type = error.bond_type;
  }
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"9874b6974fcf056f0a7313671fcb32d818b00e9d66c7e11934a5919569440c65\","
               "\"operation\":\"_Valence5NCleanUp8\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(input);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(input);
  std::cout << ",\"atom_index\":" << test_case.atom_index
            << "},\"output\":{\"status\":\""
            << (returned ? "return" : "exception") << "\",\"result\":";
  if (returned) {
    std::cout << (result ? "true" : "false");
  } else {
    std::cout << "null";
  }
  std::cout << ",\"exception\":";
  if (threw) {
    std::cout << "{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << '}';
  } else {
    std::cout << "null";
  }
  std::cout << ",\"diagnostics\":";
  if (threw) {
    std::cout << "[{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "}]";
  } else {
    std::cout << "[]";
  }
  std::cout << ",\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[]}}\n";
}

static int valence5n_cleanup8_records() {
  using T = RDKit::Bond::BondType;
  std::vector<Valence5NCleanup8Case> cases;

  cases.push_back(make_valence5n_cleanup8_case("wrong-element", 6));
  auto wrong_charge = make_valence5n_cleanup8_case("wrong-charge", 7);
  wrong_charge.atoms[5].formal_charge = 1;
  cases.push_back(std::move(wrong_charge));
  auto wrong_valence = make_valence5n_cleanup8_case("wrong-valence", 7);
  wrong_valence.atoms[5].num_explicit_hydrogens = 2;
  cases.push_back(std::move(wrong_valence));
  cases.push_back({"unsupported-valence-bond-exception",
                   {{7, 0, 0}, {6, 0, 0}},
                   {{0, 1, T::THREECENTER}},
                   0});

  for (unsigned int atom_index = 0; atom_index < 5; ++atom_index) {
    auto test_case = make_valence5n_cleanup8_case(
        "atom-predicate-miss-" + std::to_string(atom_index), 7);
    test_case.atoms[atom_index].atomic_number = 8;
    cases.push_back(std::move(test_case));
  }
  for (unsigned int bond_index = 0; bond_index < 5; ++bond_index) {
    auto test_case = make_valence5n_cleanup8_case(
        "bond-predicate-miss-" + std::to_string(bond_index), 7);
    test_case.bonds[bond_index].bond_type =
        test_case.bonds[bond_index].bond_type == T::SINGLE ? T::DOUBLE
                                                           : T::SINGLE;
    cases.push_back(std::move(test_case));
  }
  auto root_bond_miss =
      make_valence5n_cleanup8_case("root-bond-predicate-miss", 7);
  root_bond_miss.bonds[5].bond_type = T::SINGLE;
  root_bond_miss.atoms[5].num_explicit_hydrogens = 4;
  cases.push_back(std::move(root_bond_miss));

  cases.push_back(make_valence5n_cleanup8_case("unique-match", 7));

  auto multiple =
      make_valence5n_cleanup8_case("multiple-unique-atom-set-matches", 7);
  const auto remote = make_valence5n_cleanup8_case("unused", 50);
  append_valence5n_cleanup8_component(multiple, remote);
  cases.push_back(std::move(multiple));

  Valence5NCleanup8Case unrelated = {
      "unique-match-does-not-contain-argument-atom",
      {{7, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0}},
      {{0, 1, T::DOUBLE}, {0, 2, T::SINGLE}, {0, 3, T::SINGLE},
       {0, 4, T::SINGLE}},
      0};
  append_valence5n_cleanup8_component(unrelated, remote);
  cases.push_back(std::move(unrelated));

  for (const auto &test_case : cases) {
    print_valence5n_cleanup8_case(test_case);
  }
  return 0;
}

static Valence5NCleanup9Case make_valence5n_cleanup9_case(
    std::string case_id, int root_atomic_number) {
  using D = RDKit::Bond::BondDir;
  using T = RDKit::Bond::BondType;
  return {std::move(case_id),
          {{6, -2, 0}, {7, 3, 0}, {7, -4, 0},
           {6, 5, 0},  {6, -6, 0}, {root_atomic_number, 0, 3}},
          {{1, 0, T::SINGLE, D::BEGINWEDGE},
           {1, 2, T::DOUBLE, D::BEGINDASH},
           {3, 2, T::SINGLE, D::ENDDOWNRIGHT},
           {3, 4, T::DOUBLE, D::ENDUPRIGHT},
           {0, 4, T::SINGLE, D::EITHERDOUBLE},
           {5, 0, T::DOUBLE, D::UNKNOWN}},
          5};
}

static void append_valence5n_cleanup9_component(
    Valence5NCleanup9Case &destination,
    const Valence5NCleanup9Case &component) {
  const auto offset = static_cast<unsigned int>(destination.atoms.size());
  destination.atoms.insert(destination.atoms.end(), component.atoms.begin(),
                           component.atoms.end());
  for (auto bond : component.bonds) {
    bond.begin_atom_index += offset;
    bond.end_atom_index += offset;
    destination.bonds.push_back(bond);
  }
}

static void print_valence5n_cleanup9_case(
    const Valence5NCleanup9Case &test_case) {
  RDKit::RWMol input(test_case.atoms, test_case.bonds);
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);

  bool returned = false;
  bool result = false;
  bool threw = false;
  unsigned int error_bond_index = 0;
  int error_bond_type = 0;
  try {
    result = RDKit::_Valence5NCleanUp9(
        molecule, molecule.getAtomWithIdx(test_case.atom_index));
    returned = true;
  } catch (const RDKit::ValueErrorException &error) {
    threw = true;
    error_bond_index = error.bond_index;
    error_bond_type = error.bond_type;
  }
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"c2b7f87dfdd510b485982edfc6a4f0327ca042e86b5e35eee7484e2193096ef0\","
               "\"operation\":\"_Valence5NCleanUp9\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(input);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(input);
  std::cout << ",\"atom_index\":" << test_case.atom_index
            << "},\"output\":{\"status\":\""
            << (returned ? "return" : "exception") << "\",\"result\":";
  if (returned) {
    std::cout << (result ? "true" : "false");
  } else {
    std::cout << "null";
  }
  std::cout << ",\"exception\":";
  if (threw) {
    std::cout << "{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << '}';
  } else {
    std::cout << "null";
  }
  std::cout << ",\"diagnostics\":";
  if (threw) {
    std::cout << "[{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "}]";
  } else {
    std::cout << "[]";
  }
  std::cout << ",\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[]}}\n";
}

static int valence5n_cleanup9_records() {
  using T = RDKit::Bond::BondType;
  std::vector<Valence5NCleanup9Case> cases;

  cases.push_back(make_valence5n_cleanup9_case("wrong-element", 6));
  auto wrong_charge = make_valence5n_cleanup9_case("wrong-charge", 7);
  wrong_charge.atoms[5].formal_charge = 1;
  cases.push_back(std::move(wrong_charge));
  auto wrong_valence = make_valence5n_cleanup9_case("wrong-valence", 7);
  wrong_valence.atoms[5].num_explicit_hydrogens = 2;
  cases.push_back(std::move(wrong_valence));
  cases.push_back({"unsupported-valence-bond-exception",
                   {{7, 0, 0}, {6, 0, 0}},
                   {{0, 1, T::THREECENTER}},
                   0});

  for (unsigned int atom_index = 0; atom_index < 5; ++atom_index) {
    auto test_case = make_valence5n_cleanup9_case(
        "atom-predicate-miss-" + std::to_string(atom_index), 7);
    test_case.atoms[atom_index].atomic_number = 8;
    cases.push_back(std::move(test_case));
  }
  for (unsigned int bond_index = 0; bond_index < 5; ++bond_index) {
    auto test_case = make_valence5n_cleanup9_case(
        "bond-predicate-miss-" + std::to_string(bond_index), 7);
    test_case.bonds[bond_index].bond_type =
        test_case.bonds[bond_index].bond_type == T::SINGLE ? T::DOUBLE
                                                           : T::SINGLE;
    cases.push_back(std::move(test_case));
  }
  auto root_bond_miss =
      make_valence5n_cleanup9_case("root-bond-predicate-miss", 7);
  root_bond_miss.bonds[5].bond_type = T::SINGLE;
  root_bond_miss.atoms[5].num_explicit_hydrogens = 4;
  cases.push_back(std::move(root_bond_miss));

  cases.push_back(make_valence5n_cleanup9_case("unique-match", 7));

  auto multiple =
      make_valence5n_cleanup9_case("multiple-unique-atom-set-matches", 7);
  const auto remote = make_valence5n_cleanup9_case("unused", 50);
  append_valence5n_cleanup9_component(multiple, remote);
  cases.push_back(std::move(multiple));

  Valence5NCleanup9Case unrelated = {
      "unique-match-does-not-contain-argument-atom",
      {{7, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0}},
      {{0, 1, T::DOUBLE}, {0, 2, T::SINGLE}, {0, 3, T::SINGLE},
       {0, 4, T::SINGLE}},
      0};
  append_valence5n_cleanup9_component(unrelated, remote);
  cases.push_back(std::move(unrelated));

  for (const auto &test_case : cases) {
    print_valence5n_cleanup9_case(test_case);
  }
  return 0;
}

static Valence5NCleanupACase make_valence5n_cleanupa_path_case(
    std::string case_id, unsigned int path_length) {
  using T = RDKit::Bond::BondType;
  const unsigned int target_index = path_length;
  const unsigned int partner_index = target_index + 1;
  Valence5NCleanupACase test_case{std::move(case_id), {{7, 0, 0}}, {}, 0};
  for (unsigned int index = 1; index < path_length; ++index) {
    test_case.atoms.push_back({6, 0, 0});
  }
  test_case.atoms.insert(test_case.atoms.end(),
                         {{7, 0, 0}, {7, 0, 0}, {6, 0, 0}, {6, 0, 0},
                          {6, 0, 0}});
  for (unsigned int index = 0; index < path_length; ++index) {
    test_case.bonds.push_back(
        {index, index + 1, index % 2 == 0 ? T::DOUBLE : T::SINGLE});
  }
  test_case.bonds.push_back({target_index, partner_index, T::DOUBLE});
  test_case.bonds.insert(test_case.bonds.end(),
                         {{0, partner_index + 1, T::SINGLE},
                          {0, partner_index + 2, T::SINGLE},
                          {0, partner_index + 3, T::SINGLE}});
  return test_case;
}

static void append_valence5n_cleanupa_component(
    Valence5NCleanupACase &destination,
    const Valence5NCleanupACase &component) {
  const auto offset = static_cast<unsigned int>(destination.atoms.size());
  destination.atoms.insert(destination.atoms.end(), component.atoms.begin(),
                           component.atoms.end());
  for (auto bond : component.bonds) {
    bond.begin_atom_index += offset;
    bond.end_atom_index += offset;
    destination.bonds.push_back(bond);
  }
}

static void print_valence5n_cleanupa_case(
    const Valence5NCleanupACase &test_case) {
  RDKit::RWMol input(test_case.atoms, test_case.bonds);
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);

  bool returned = false;
  bool result = false;
  bool threw = false;
  unsigned int error_bond_index = 0;
  int error_bond_type = 0;
  try {
    result = RDKit::_Valence5NCleanUpA(
        molecule, molecule.getAtomWithIdx(test_case.atom_index));
    returned = true;
  } catch (const RDKit::ValueErrorException &error) {
    threw = true;
    error_bond_index = error.bond_index;
    error_bond_type = error.bond_type;
  }
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"2d4518bea1a8abb8d5736f1c784fa7e30cb033950dcd6ba5da13846bb46d6392\","
               "\"operation\":\"_Valence5NCleanUpA\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(input);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(input);
  std::cout << ",\"atom_index\":" << test_case.atom_index
            << "},\"output\":{\"status\":\""
            << (returned ? "return" : "exception") << "\",\"result\":";
  if (returned) {
    std::cout << (result ? "true" : "false");
  } else {
    std::cout << "null";
  }
  std::cout << ",\"exception\":";
  if (threw) {
    std::cout << "{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << '}';
  } else {
    std::cout << "null";
  }
  std::cout << ",\"diagnostics\":";
  if (threw) {
    std::cout << "[{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "}]";
  } else {
    std::cout << "[]";
  }
  std::cout << ",\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[]}}\n";
}

static int valence5n_cleanupa_records() {
  using D = RDKit::Bond::BondDir;
  using T = RDKit::Bond::BondType;
  std::vector<Valence5NCleanupACase> cases;

  cases.push_back(
      {"wrong-element", {{6, 0, 0}, {6, 0, 0}}, {{0, 1, T::THREECENTER}}, 0});
  auto wrong_charge = make_valence5n_cleanupa_path_case("wrong-charge", 1);
  wrong_charge.atoms[0].formal_charge = 1;
  cases.push_back(std::move(wrong_charge));
  auto wrong_valence = make_valence5n_cleanupa_path_case("wrong-valence", 1);
  wrong_valence.atoms[0].num_explicit_hydrogens = 1;
  cases.push_back(std::move(wrong_valence));
  cases.push_back({"unsupported-valence-bond-exception",
                   {{7, 0, 0}, {6, 0, 0}},
                   {{0, 1, T::THREECENTER}},
                   0});
  cases.push_back({"zero-substructure-matches",
                   {{7, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0}},
                   {{0, 1, T::DOUBLE}, {0, 2, T::SINGLE},
                    {0, 3, T::SINGLE}, {0, 4, T::SINGLE}},
                   0});
  cases.push_back({"only-match-contains-root",
                   {{7, 0, 0}, {7, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0}},
                   {{1, 0, T::DOUBLE}, {0, 2, T::SINGLE},
                    {0, 3, T::SINGLE}, {0, 4, T::SINGLE}},
                   0});

  Valence5NCleanupACase no_path = {
      "remote-match-without-path",
      {{7, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0}},
      {{0, 1, T::DOUBLE}, {0, 2, T::SINGLE}, {0, 3, T::SINGLE},
       {0, 4, T::SINGLE}},
      0};
  append_valence5n_cleanupa_component(
      no_path, {"unused", {{7, 0, 0}, {7, 0, 0}}, {{1, 0, T::DOUBLE}}, 0});
  cases.push_back(std::move(no_path));

  auto direct = make_valence5n_cleanupa_path_case("direct-path", 1);
  direct.bonds[0].direction = D::BEGINWEDGE;
  direct.bonds[1].direction = D::BEGINDASH;
  direct.bonds[2].direction = D::ENDDOWNRIGHT;
  direct.bonds[3].direction = D::ENDUPRIGHT;
  cases.push_back(std::move(direct));
  cases.push_back(make_valence5n_cleanupa_path_case("path-length-9", 9));
  cases.push_back(make_valence5n_cleanupa_path_case("path-length-11", 11));

  cases.push_back({"later-shorter-path",
                   {{7, 0, 0}, {6, 0, 0}, {6, 0, 0}, {7, 0, 0},
                    {7, 0, 0}, {7, 0, 0}, {7, 0, 0}, {6, 0, 0}},
                   {{0, 1, T::DOUBLE}, {1, 2, T::SINGLE},
                    {2, 3, T::DOUBLE}, {3, 4, T::DOUBLE},
                    {0, 5, T::DOUBLE}, {5, 6, T::DOUBLE},
                    {0, 7, T::SINGLE}},
                   0});
  cases.push_back({"equal-path-retains-first",
                   {{7, 0, 0}, {6, 0, 0}, {6, 0, 0}, {7, 0, 0},
                    {7, 0, 0}, {6, 0, 0}, {6, 0, 0}, {7, 0, 0},
                    {7, 0, 0}, {6, 0, 0}},
                   {{0, 1, T::DOUBLE}, {1, 2, T::SINGLE},
                    {2, 3, T::DOUBLE}, {3, 4, T::DOUBLE},
                    {0, 5, T::DOUBLE}, {5, 6, T::SINGLE},
                    {6, 7, T::DOUBLE}, {7, 8, T::DOUBLE},
                    {0, 9, T::SINGLE}},
                   0});

  for (const auto &test_case : cases) {
    print_valence5n_cleanupa_case(test_case);
  }
  return 0;
}

static void print_valence5n_cleanupb_case(
    const Valence5NCleanupBCase &test_case) {
  RDKit::RWMol input(test_case.atoms, test_case.bonds);
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);

  bool returned = false;
  bool result = false;
  bool threw = false;
  unsigned int error_bond_index = 0;
  int error_bond_type = 0;
  try {
    result = RDKit::_Valence5NCleanUpB(
        molecule, molecule.getAtomWithIdx(test_case.atom_index));
    returned = true;
  } catch (const RDKit::ValueErrorException &error) {
    threw = true;
    error_bond_index = error.bond_index;
    error_bond_type = error.bond_type;
  }
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"9d097bf179047ec55e31d8f45ffe573f449e3d5828833e5430d2faf355a5defc\","
               "\"operation\":\"_Valence5NCleanUpB\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(input);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(input);
  std::cout << ",\"atom_index\":" << test_case.atom_index
            << "},\"output\":{\"status\":\""
            << (returned ? "return" : "exception") << "\",\"result\":";
  if (returned) {
    std::cout << (result ? "true" : "false");
  } else {
    std::cout << "null";
  }
  std::cout << ",\"exception\":";
  if (threw) {
    std::cout << "{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << '}';
  } else {
    std::cout << "null";
  }
  std::cout << ",\"diagnostics\":";
  if (threw) {
    std::cout << "[{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "}]";
  } else {
    std::cout << "[]";
  }
  std::cout << ",\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[]}}\n";
}

static int valence5n_cleanupb_records() {
  using D = RDKit::Bond::BondDir;
  using T = RDKit::Bond::BondType;
  const std::vector<Valence5NCleanupBCase> cases = {
      {"no-target",
       {{7, -4, 0}, {6, 0, 0}, {8, 0, 0}},
       {{0, 1, T::SINGLE}, {0, 2, T::THREECENTER}},
       0},
      {"charged-carbon", {{16, 8, 0}, {6, 1, 0}}, {{0, 1, T::DOUBLE}}, 0},
      {"over-depth",
       {{7, 0, 0}, {7, 0, 0}, {6, 0, 0}},
       {{0, 1, T::DOUBLE}, {1, 2, T::DOUBLE}},
       0},
      {"success",
       {{8, -5, 2}, {6, 0, 3}},
       {{0, 1, T::DOUBLE, D::BEGINWEDGE}},
       0},
      {"first-in-adjacency",
       {{15, 3, 0}, {6, 0, 0}, {6, 0, 0}},
       {{0, 2, T::DOUBLE, D::BEGINDASH},
        {0, 1, T::DOUBLE, D::ENDDOWNRIGHT}},
       0},
      {"target-valence-exception",
       {{7, 0, 0}, {6, 0, 0}, {8, 0, 0}},
       {{0, 1, T::DOUBLE}, {1, 2, T::THREECENTER}},
       0},
      {"root-valence-exception",
       {{7, -3, 0}, {6, 0, 0}, {8, 0, 0}},
       {{0, 1, T::DOUBLE}, {0, 2, T::THREECENTER}},
       0},
  };
  for (const auto &test_case : cases) {
    print_valence5n_cleanupb_case(test_case);
  }
  return 0;
}

static void print_valence7s_cleanup1_case(
    const Valence7SCleanup1Case &test_case) {
  RDKit::RWMol input(test_case.atoms, test_case.bonds);
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);

  bool returned = false;
  bool result = false;
  bool threw = false;
  unsigned int error_bond_index = 0;
  int error_bond_type = 0;
  try {
    result = RDKit::_Valence7SCleanUp1(
        molecule, molecule.getAtomWithIdx(test_case.atom_index));
    returned = true;
  } catch (const RDKit::ValueErrorException &error) {
    threw = true;
    error_bond_index = error.bond_index;
    error_bond_type = error.bond_type;
  }
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"a03efbad1fb09fa47513380a625fefccda34d23e3d4f5fd0b2f2b2cf3f8d076e\","
               "\"operation\":\"_Valence7SCleanUp1\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(input);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(input);
  std::cout << ",\"atom_index\":" << test_case.atom_index
            << "},\"output\":{\"status\":\""
            << (returned ? "return" : "exception") << "\",\"result\":";
  if (returned) {
    std::cout << (result ? "true" : "false");
  } else {
    std::cout << "null";
  }
  std::cout << ",\"exception\":";
  if (threw) {
    std::cout << "{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << '}';
  } else {
    std::cout << "null";
  }
  std::cout << ",\"diagnostics\":";
  if (threw) {
    std::cout << "[{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "}]";
  } else {
    std::cout << "[]";
  }
  std::cout << ",\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[]}}\n";
}

static int valence7s_cleanup1_records() {
  using D = RDKit::Bond::BondDir;
  using T = RDKit::Bond::BondType;
  const std::vector<Valence7SCleanup1Case> cases = {
      {"wrong-element", {{15, -1, 0}, {8, 0, 0}},
       {{0, 1, T::THREECENTER}}, 0},
      {"wrong-charge", {{16, 0, 0}, {8, 0, 0}},
       {{0, 1, T::THREECENTER}}, 0},
      {"wrong-valence", {{16, -1, 0}, {8, 0, 0}},
       {{0, 1, T::DOUBLE}}, 0},
      {"entry-valence-exception", {{16, -1, 0}, {8, 0, 0}},
       {{0, 1, T::THREECENTER}}, 0},
      {"oxygen-nondouble", {{16, -1, 6}, {8, 0, 0}},
       {{0, 1, T::SINGLE}}, 0},
      {"carbon-nonsingle", {{16, -1, 5}, {6, 0, 0}},
       {{0, 1, T::DOUBLE}}, 0},
      {"other-element", {{16, -1, 6}, {7, 0, 0}},
       {{0, 1, T::SINGLE}}, 0},
      {"no-oxygen", {{16, -1, 6}, {6, 0, 0}},
       {{0, 1, T::SINGLE}}, 0},
      {"criteria-false", {{16, -1, 3}, {8, 0, 0}, {8, 0, 0}},
       {{0, 1, T::DOUBLE}, {0, 2, T::DOUBLE}}, 0},
      {"one-carbon", {{16, -1, 4}, {8, 4, 2}, {6, -5, 0}},
       {{0, 1, T::DOUBLE, D::BEGINWEDGE},
        {0, 2, T::SINGLE, D::BEGINDASH}},
       0},
      {"three-oxygen",
       {{16, -1, 1}, {8, 2, 0}, {8, 3, 0}, {8, 4, 0}},
       {{0, 2, T::DOUBLE, D::ENDDOWNRIGHT},
        {0, 1, T::DOUBLE, D::ENDUPRIGHT}, {0, 3, T::DOUBLE}},
       0},
      {"oxygen-valence-exception",
       {{16, -1, 4}, {6, 0, 0}, {8, 0, 0}, {7, 0, 0}},
       {{0, 1, T::SINGLE}, {0, 2, T::DOUBLE},
        {2, 3, T::THREECENTER}},
       0},
  };
  for (const auto &test_case : cases) {
    print_valence7s_cleanup1_case(test_case);
  }
  return 0;
}

static void print_valence7s_cleanup2_case(
    const Valence7SCleanup2Case &test_case) {
  RDKit::RWMol input(test_case.atoms, test_case.bonds);
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);

  bool returned = false;
  bool result = false;
  bool threw = false;
  unsigned int error_bond_index = 0;
  int error_bond_type = 0;
  try {
    result = RDKit::_Valence7SCleanUp2(
        molecule, molecule.getAtomWithIdx(test_case.atom_index));
    returned = true;
  } catch (const RDKit::ValueErrorException &error) {
    threw = true;
    error_bond_index = error.bond_index;
    error_bond_type = error.bond_type;
  }
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"d8e8ca7577874896efc2f1df36738d0bee34a8dcabb35493be0095017a2e9618\","
               "\"operation\":\"_Valence7SCleanUp2\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(input);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(input);
  std::cout << ",\"atom_index\":" << test_case.atom_index
            << "},\"output\":{\"status\":\""
            << (returned ? "return" : "exception") << "\",\"result\":";
  if (returned) {
    std::cout << (result ? "true" : "false");
  } else {
    std::cout << "null";
  }
  std::cout << ",\"exception\":";
  if (threw) {
    std::cout << "{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << '}';
  } else {
    std::cout << "null";
  }
  std::cout << ",\"diagnostics\":";
  if (threw) {
    std::cout << "[{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "}]";
  } else {
    std::cout << "[]";
  }
  std::cout << ",\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[]}}\n";
}

static int valence7s_cleanup2_records() {
  using D = RDKit::Bond::BondDir;
  using T = RDKit::Bond::BondType;
  const std::vector<Valence7SCleanup2Case> cases = {
      {"wrong-element", {{15, -1, 0}, {7, 0, 0}},
       {{0, 1, T::THREECENTER}}, 0},
      {"wrong-charge", {{16, 0, 0}, {7, 0, 0}},
       {{0, 1, T::THREECENTER}}, 0},
      {"wrong-valence", {{16, -1, 0}, {7, 0, 0}},
       {{0, 1, T::DOUBLE}}, 0},
      {"entry-valence-exception", {{16, -1, 0}, {7, 0, 0}},
       {{0, 1, T::THREECENTER}}, 0},
      {"no-target", {{16, -1, 5}, {6, 0, 0}}, {{0, 1, T::DOUBLE}}, 0},
      {"over-depth",
       {{16, -1, 5}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0}, {7, 0, 0}},
       {{0, 1, T::DOUBLE}, {1, 2, T::SINGLE}, {2, 3, T::DOUBLE},
        {3, 4, T::TRIPLE}},
       0},
      {"direct",
       {{16, -1, 4}, {7, 0, 3}},
       {{0, 1, T::TRIPLE, D::BEGINWEDGE}},
       0},
      {"two-bonds",
       {{16, -1, 5}, {6, -3, 0}, {7, 0, 0}},
       {{0, 1, T::DOUBLE, D::BEGINDASH},
        {1, 2, T::TRIPLE, D::ENDDOWNRIGHT}},
       0},
      {"at-limit",
       {{16, -1, 5}, {6, 2, 0}, {6, -4, 0}, {7, 0, 0}},
       {{0, 1, T::DOUBLE, D::BEGINWEDGE},
        {1, 2, T::SINGLE, D::ENDUPRIGHT},
        {2, 3, T::TRIPLE, D::BEGINDASH}},
       0},
      {"first-equal-path",
       {{16, -1, 1}, {7, 0, 0}, {7, 0, 0}},
       {{0, 2, T::TRIPLE, D::ENDDOWNRIGHT},
        {0, 1, T::TRIPLE, D::ENDUPRIGHT}},
       0},
  };
  for (const auto &test_case : cases) {
    print_valence7s_cleanup2_case(test_case);
  }
  return 0;
}

static void print_valence7s_cleanup3_case(
    const Valence7SCleanup3Case &test_case) {
  RDKit::RWMol input(test_case.atoms, test_case.bonds);
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);

  bool returned = false;
  bool result = false;
  bool threw = false;
  unsigned int error_bond_index = 0;
  int error_bond_type = 0;
  try {
    result = RDKit::_Valence7SCleanUp3(
        molecule, molecule.getAtomWithIdx(test_case.atom_index));
    returned = true;
  } catch (const RDKit::ValueErrorException &error) {
    threw = true;
    error_bond_index = error.bond_index;
    error_bond_type = error.bond_type;
  }
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"b61785dfecb61f3b567183f96ec1290ef56f576436bfcd07a808907ebbf4b4ed\","
               "\"operation\":\"_Valence7SCleanUp3\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(input);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(input);
  std::cout << ",\"atom_index\":" << test_case.atom_index
            << "},\"output\":{\"status\":\""
            << (returned ? "return" : "exception") << "\",\"result\":";
  if (returned) {
    std::cout << (result ? "true" : "false");
  } else {
    std::cout << "null";
  }
  std::cout << ",\"exception\":";
  if (threw) {
    std::cout << "{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << '}';
  } else {
    std::cout << "null";
  }
  std::cout << ",\"diagnostics\":";
  if (threw) {
    std::cout << "[{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "}]";
  } else {
    std::cout << "[]";
  }
  std::cout << ",\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[]}}\n";
}

static int valence7s_cleanup3_records() {
  using D = RDKit::Bond::BondDir;
  using T = RDKit::Bond::BondType;
  const std::vector<Valence7SCleanup3Case> cases = {
      {"wrong-element", {{15, -1, 0}, {7, 0, 0}},
       {{0, 1, T::THREECENTER}}, 0},
      {"wrong-charge", {{16, 0, 0}, {7, 0, 0}},
       {{0, 1, T::THREECENTER}}, 0},
      {"wrong-valence", {{16, -1, 0}, {7, 0, 0}},
       {{0, 1, T::DOUBLE}}, 0},
      {"entry-valence-exception", {{16, -1, 0}, {7, 0, 0}},
       {{0, 1, T::THREECENTER}}, 0},
      {"no-target", {{16, -1, 5}, {6, 0, 0}}, {{0, 1, T::DOUBLE}}, 0},
      {"charged-target", {{16, -1, 5}, {7, 1, 0}},
       {{0, 1, T::DOUBLE}}, 0},
      {"over-depth",
       {{16, -1, 5}, {6, 0, 0}, {7, 0, 0}},
       {{0, 1, T::DOUBLE}, {1, 2, T::DOUBLE}}, 0},
      {"success",
       {{16, -1, 5}, {7, 0, 3}},
       {{0, 1, T::DOUBLE, D::BEGINWEDGE}}, 0},
      {"first-equal-path",
       {{16, -1, 3}, {7, 0, 0}, {7, 0, 0}},
       {{0, 2, T::DOUBLE, D::BEGINDASH},
        {0, 1, T::DOUBLE, D::ENDDOWNRIGHT}},
       0},
      {"target-valence-not-recomputed",
       {{16, -1, 5}, {7, 0, 0}, {6, -4, 0}},
       {{0, 1, T::DOUBLE, D::ENDUPRIGHT}, {1, 2, T::THREECENTER}},
       0},
  };
  for (const auto &test_case : cases) {
    print_valence7s_cleanup3_case(test_case);
  }
  return 0;
}

static void print_valence8s_cleanup1_case(
    const Valence8SCleanup1Case &test_case) {
  RDKit::RWMol input(test_case.atoms, test_case.bonds);
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);

  bool returned = false;
  bool result = false;
  bool threw = false;
  unsigned int error_bond_index = 0;
  int error_bond_type = 0;
  try {
    result = RDKit::_Valence8SCleanUp1(
        molecule, molecule.getAtomWithIdx(test_case.atom_index));
    returned = true;
  } catch (const RDKit::ValueErrorException &error) {
    threw = true;
    error_bond_index = error.bond_index;
    error_bond_type = error.bond_type;
  }
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"e8bda5f9c53e83b3f5a0cce7f26861384606e709e5ef66e7adedac5c0ebcfd78\","
               "\"operation\":\"_Valence8SCleanUp1\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(input);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(input);
  std::cout << ",\"atom_index\":" << test_case.atom_index
            << "},\"output\":{\"status\":\""
            << (returned ? "return" : "exception") << "\",\"result\":";
  if (returned) {
    std::cout << (result ? "true" : "false");
  } else {
    std::cout << "null";
  }
  std::cout << ",\"exception\":";
  if (threw) {
    std::cout << "{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << '}';
  } else {
    std::cout << "null";
  }
  std::cout << ",\"diagnostics\":";
  if (threw) {
    std::cout << "[{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "}]";
  } else {
    std::cout << "[]";
  }
  std::cout << ",\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[]}}\n";
}

static int valence8s_cleanup1_records() {
  using D = RDKit::Bond::BondDir;
  using T = RDKit::Bond::BondType;
  const std::vector<Valence8SCleanup1Case> cases = {
      {"wrong-element", {{15, -1, 0}, {7, 0, 0}},
       {{0, 1, T::THREECENTER}}, 0},
      {"wrong-charge", {{16, 0, 0}, {7, 0, 0}},
       {{0, 1, T::THREECENTER}}, 0},
      {"wrong-valence", {{16, -1, 0}, {7, 0, 0}},
       {{0, 1, T::DOUBLE}}, 0},
      {"entry-valence-exception", {{16, -1, 0}, {7, 0, 0}},
       {{0, 1, T::THREECENTER}}, 0},
      {"no-target", {{16, -1, 5}, {6, 0, 0}}, {{0, 1, T::DOUBLE}}, 0},
      {"charged-target", {{16, -1, 5}, {7, 1, 0}},
       {{0, 1, T::DOUBLE}}, 0},
      {"direct",
       {{16, -1, 5}, {7, 0, 3}},
       {{0, 1, T::DOUBLE, D::BEGINWEDGE}}, 0},
      {"alternating",
       {{16, -1, 5}, {6, 4, 0}, {6, -3, 0}, {7, 0, 2}},
       {{0, 1, T::DOUBLE, D::BEGINDASH},
        {1, 2, T::SINGLE, D::ENDDOWNRIGHT},
        {2, 3, T::DOUBLE, D::ENDUPRIGHT}},
       0},
      {"at-limit",
       {{16, -1, 5}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0},
        {6, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0}, {7, 0, 4}},
       {{0, 1, T::DOUBLE, D::ENDUPRIGHT},
        {1, 2, T::SINGLE, D::ENDDOWNRIGHT},
        {2, 3, T::DOUBLE, D::ENDUPRIGHT},
        {3, 4, T::SINGLE, D::ENDDOWNRIGHT},
        {4, 5, T::DOUBLE, D::ENDUPRIGHT},
        {5, 6, T::SINGLE, D::ENDDOWNRIGHT},
        {6, 7, T::DOUBLE, D::ENDUPRIGHT},
        {7, 8, T::SINGLE, D::ENDDOWNRIGHT},
        {8, 9, T::DOUBLE, D::ENDUPRIGHT}},
       0},
      {"over-limit",
       {{16, -1, 5}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0},
        {6, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0}, {6, 0, 0},
        {6, 0, 0}, {7, 0, 0}},
       {{0, 1, T::DOUBLE}, {1, 2, T::SINGLE}, {2, 3, T::DOUBLE},
        {3, 4, T::SINGLE}, {4, 5, T::DOUBLE}, {5, 6, T::SINGLE},
        {6, 7, T::DOUBLE}, {7, 8, T::SINGLE}, {8, 9, T::DOUBLE},
        {9, 10, T::SINGLE}, {10, 11, T::DOUBLE}},
       0},
      {"first-equal-path",
       {{16, -1, 3}, {7, 0, 0}, {7, 0, 2}},
       {{0, 2, T::DOUBLE, D::ENDUPRIGHT},
        {0, 1, T::DOUBLE, D::ENDDOWNRIGHT}},
       0},
      {"target-valence-exception",
       {{16, -1, 5}, {7, 0, 4}, {6, -4, 0}},
       {{0, 1, T::DOUBLE, D::BEGINDASH}, {1, 2, T::THREECENTER}},
       0},
  };
  for (const auto &test_case : cases) {
    print_valence8s_cleanup1_case(test_case);
  }
  return 0;
}

static void print_valence8cl_cleanup1_case(
    const Valence8ClCleanup1Case &test_case) {
  RDKit::RWMol input(test_case.atoms, test_case.bonds);
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);

  bool returned = false;
  bool result = false;
  bool threw = false;
  unsigned int error_bond_index = 0;
  int error_bond_type = 0;
  try {
    result = RDKit::_Valence8ClCleanUp1(
        molecule, molecule.getAtomWithIdx(test_case.atom_index));
    returned = true;
  } catch (const RDKit::ValueErrorException &error) {
    threw = true;
    error_bond_index = error.bond_index;
    error_bond_type = error.bond_type;
  }
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"8e97a538eab77c78c2fb02be47d47cf6be989761cc280cf4857b9500885b1f3a\","
               "\"operation\":\"_Valence8ClCleanUp1\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(input);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(input);
  std::cout << ",\"atom_index\":" << test_case.atom_index
            << "},\"output\":{\"status\":\""
            << (returned ? "return" : "exception") << "\",\"result\":";
  if (returned) {
    std::cout << (result ? "true" : "false");
  } else {
    std::cout << "null";
  }
  std::cout << ",\"exception\":";
  if (threw) {
    std::cout << "{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << '}';
  } else {
    std::cout << "null";
  }
  std::cout << ",\"diagnostics\":";
  if (threw) {
    std::cout << "[{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "}]";
  } else {
    std::cout << "[]";
  }
  std::cout << ",\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[]}}\n";
}

static int valence8cl_cleanup1_records() {
  using D = RDKit::Bond::BondDir;
  using T = RDKit::Bond::BondType;
  const std::vector<Valence8ClCleanup1Case> cases = {
      {"wrong-valence", {{17, -1, 7}}, {}, 0},
      {"wrong-charge", {{17, 0, 8}}, {}, 0},
      {"valence-before-charge-exception",
       {{17, 0, 0}, {8, 0, 0}},
       {{0, 1, T::THREECENTER}}, 0},
      {"non-oxygen", {{17, -1, 6}, {7, 0, 0}},
       {{0, 1, T::DOUBLE}}, 0},
      {"empty-neighbors-no-element-guard", {{15, -1, 8}}, {}, 0},
      {"mixed-bonds",
       {{17, -1, 5}, {8, 4, 2}, {8, 5, 3}},
       {{0, 1, T::DOUBLE, D::BEGINWEDGE},
        {0, 2, T::SINGLE, D::BEGINDASH}},
       0},
      {"ordered-partial-error",
       {{17, -1, 4}, {8, 2, 0}, {8, 3, 0}, {6, 0, 0}},
       {{0, 1, T::DOUBLE, D::ENDDOWNRIGHT},
        {0, 2, T::DOUBLE, D::BEGINDASH},
        {2, 3, T::THREECENTER}},
       0},
  };
  for (const auto &test_case : cases) {
    print_valence8cl_cleanup1_case(test_case);
  }
  return 0;
}

static void print_valence5cl_cleanup1_case(
    const Valence5ClCleanup1Case &test_case) {
  RDKit::RWMol input(test_case.atoms, test_case.bonds);
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);

  bool returned = false;
  bool result = false;
  bool threw = false;
  unsigned int error_bond_index = 0;
  int error_bond_type = 0;
  try {
    result = RDKit::_Valence5ClCleanUp1(
        molecule, molecule.getAtomWithIdx(test_case.atom_index));
    returned = true;
  } catch (const RDKit::ValueErrorException &error) {
    threw = true;
    error_bond_index = error.bond_index;
    error_bond_type = error.bond_type;
  }
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"102c9d4a05c61706c9b0de30d61fc725fc394f80d532810d0b94b0103d6ff250\","
               "\"operation\":\"_Valence5ClCleanUp1\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(input);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(input);
  std::cout << ",\"atom_index\":" << test_case.atom_index
            << "},\"output\":{\"status\":\""
            << (returned ? "return" : "exception") << "\",\"result\":";
  if (returned) {
    std::cout << (result ? "true" : "false");
  } else {
    std::cout << "null";
  }
  std::cout << ",\"exception\":";
  if (threw) {
    std::cout << "{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << '}';
  } else {
    std::cout << "null";
  }
  std::cout << ",\"diagnostics\":";
  if (threw) {
    std::cout << "[{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "}]";
  } else {
    std::cout << "[]";
  }
  std::cout << ",\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[]}}\n";
}

static int valence5cl_cleanup1_records() {
  using D = RDKit::Bond::BondDir;
  using T = RDKit::Bond::BondType;
  const std::vector<Valence5ClCleanup1Case> cases = {
      {"entry-valence-exception",
       {{17, 0, 0}, {8, -1, 0}},
       {{0, 1, T::THREECENTER}}, 0},
      {"wrong-valence", {{17, 1, 4}, {8, -1, 0}},
       {{0, 1, T::SINGLE}}, 0},
      {"wrong-charge", {{17, 0, 5}, {8, -1, 0}},
       {{0, 1, T::SINGLE}}, 0},
      {"wrong-element", {{17, 1, 5}, {7, -1, 0}},
       {{0, 1, T::SINGLE}}, 0},
      {"wrong-target-charge", {{17, 1, 5}, {8, 0, 0}},
       {{0, 1, T::SINGLE}}, 0},
      {"wrong-bond", {{17, 1, 4}, {8, -1, 0}},
       {{0, 1, T::DOUBLE}}, 0},
      {"success-no-element-guard-target-not-recomputed",
       {{15, 1, 5}, {8, -1, 3}, {6, 4, 0}},
       {{0, 1, T::SINGLE, D::BEGINWEDGE},
        {1, 2, T::THREECENTER, D::BEGINDASH}},
       0},
      {"first-equal-target",
       {{17, 1, 4}, {8, -1, 2}, {8, -1, 3}},
       {{0, 1, T::SINGLE, D::ENDDOWNRIGHT},
        {0, 2, T::SINGLE, D::BEGINDASH}},
       0},
  };
  for (const auto &test_case : cases) {
    print_valence5cl_cleanup1_case(test_case);
  }
  return 0;
}

static void print_valence3cl_cleanup1_case(
    const Valence3ClCleanup1Case &test_case) {
  RDKit::RWMol input(test_case.atoms, test_case.bonds);
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);

  bool returned = false;
  bool result = false;
  bool threw = false;
  unsigned int error_bond_index = 0;
  int error_bond_type = 0;
  try {
    result = RDKit::_Valence3ClCleanUp1(
        molecule, molecule.getAtomWithIdx(test_case.atom_index));
    returned = true;
  } catch (const RDKit::ValueErrorException &error) {
    threw = true;
    error_bond_index = error.bond_index;
    error_bond_type = error.bond_type;
  }
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"078efac5f5acece4f6323b171f6aab2b1795873b55a505e56fb44aa54fa7d982\","
               "\"operation\":\"_Valence3ClCleanUp1\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(input);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(input);
  std::cout << ",\"atom_index\":" << test_case.atom_index
            << "},\"output\":{\"status\":\""
            << (returned ? "return" : "exception") << "\",\"result\":";
  if (returned) {
    std::cout << (result ? "true" : "false");
  } else {
    std::cout << "null";
  }
  std::cout << ",\"exception\":";
  if (threw) {
    std::cout << "{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << '}';
  } else {
    std::cout << "null";
  }
  std::cout << ",\"diagnostics\":";
  if (threw) {
    std::cout << "[{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "}]";
  } else {
    std::cout << "[]";
  }
  std::cout << ",\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[]}}\n";
}

static int valence3cl_cleanup1_records() {
  using D = RDKit::Bond::BondDir;
  using T = RDKit::Bond::BondType;
  const std::vector<Valence3ClCleanup1Case> cases = {
      {"entry-valence-exception",
       {{17, 1, 0}, {16, 0, 0}},
       {{0, 1, T::THREECENTER}}, 0},
      {"wrong-valence", {{17, 0, 0}, {16, 0, 0}},
       {{0, 1, T::DOUBLE}}, 0},
      {"wrong-charge", {{17, 1, 0}, {16, 0, 0}},
       {{0, 1, T::TRIPLE}}, 0},
      {"wrong-element", {{17, 0, 0}, {8, 0, 0}},
       {{0, 1, T::TRIPLE}}, 0},
      {"wrong-target-charge", {{17, 0, 0}, {16, 1, 0}},
       {{0, 1, T::TRIPLE}}, 0},
      {"wrong-bond", {{17, 0, 1}, {16, 0, 0}},
       {{0, 1, T::DOUBLE}}, 0},
      {"success-no-element-guard-target-not-recomputed",
       {{15, 0, 0}, {16, 0, 3}, {6, 4, 0}},
       {{0, 1, T::TRIPLE, D::BEGINWEDGE},
        {1, 2, T::THREECENTER, D::BEGINDASH}},
       0},
  };
  for (const auto &test_case : cases) {
    print_valence3cl_cleanup1_case(test_case);
  }
  return 0;
}

static void print_clean_up_case(const CleanUpCase &test_case) {
  RDKit::RWMol input(test_case.atoms, test_case.bonds);
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  const auto atoms_before = atom_snapshot(molecule);
  const auto bonds_before = bond_snapshot(molecule);

  bool returned = false;
  bool threw = false;
  unsigned int error_bond_index = 0;
  int error_bond_type = 0;
  try {
    RDKit::cleanUp(molecule);
    returned = true;
  } catch (const RDKit::ValueErrorException &error) {
    threw = true;
    error_bond_index = error.bond_index;
    error_bond_type = error.bond_type;
  }
  const bool graph_unchanged = atoms_before == atom_snapshot(molecule) &&
                               bonds_before == bond_snapshot(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"25458d16f3e2888aed0a60ba3b1c8f0a4e435a912215af7da93fbc491b34b371\","
               "\"operation\":\"cleanUp\",\"case_id\":\""
            << test_case.case_id << "\",\"input\":{\"atom_fields\":";
  print_atom_fields(input);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(input);
  std::cout << "},\"output\":{\"status\":\""
            << (returned ? "return" : "exception")
            << "\",\"result\":null,\"exception\":";
  if (threw) {
    std::cout << "{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << '}';
  } else {
    std::cout << "null";
  }
  std::cout << ",\"diagnostics\":";
  if (threw) {
    std::cout << "[{\"kind\":\"ValueErrorException\",\"message\":\"Bad bond type\","
                 "\"bond_index\":"
              << error_bond_index << ",\"bond_type\":" << error_bond_type
              << "}]";
  } else {
    std::cout << "[]";
  }
  std::cout << ",\"graph_unchanged\":"
            << (graph_unchanged ? "true" : "false")
            << ",\"atom_count\":" << molecule.atoms().size()
            << ",\"bond_count\":" << molecule.bonds().size()
            << ",\"atom_fields\":";
  print_atom_fields(molecule);
  std::cout << ",\"bond_fields\":";
  print_graph_bond_fields(molecule);
  std::cout << ",\"stereo_fields\":{\"bond_directions\":";
  print_directions(bond_directions(molecule));
  std::cout << "},\"properties\":[]}}\n";
}

static int clean_up_records() {
  using D = RDKit::Bond::BondDir;
  using T = RDKit::Bond::BondType;
  std::vector<CleanUpCase> cases = {
      {"empty", {}, {}},
      {"ignored-aromatic-carbon", {{6, 2, 0, true}}, {}},
      {"entry-valence-exception",
       {{7, 0, 0, true}, {6, 0, 0}},
       {{0, 1, T::THREECENTER}}},
      {"valence4-cleanup1",
       {{6, 4, 0}, {7, 2, 0}, {7, -1, 0, true}, {7, -3, 0},
        {7, 1, 0}, {8, 0, 0}},
       {{1, 0, T::SINGLE, D::BEGINDASH}, {2, 1, T::DOUBLE},
        {3, 2, T::DOUBLE}, {4, 3, T::SINGLE}, {0, 4, T::DOUBLE},
        {0, 5, T::SINGLE}}},
      {"valence4-cleanup2",
       {{7, -1, 1, true}, {7, 0, 0}, {8, -2, 0}},
       {{1, 0, T::DOUBLE}, {0, 2, T::SINGLE}}},
      {"valence4-no-match", {{7, 0, 4, true}}, {}},
      {"charged-nitrogen", {{7, 1, 0, true}}, {}},
      {"valence5-all-miss", {{7, 0, 5, true}}, {}},
      {"valence5-cleanup1",
       {{7, 0, 3, true}, {7, 1, 0}},
       {{0, 1, T::DOUBLE}}},
      {"valence5-cleanup2",
       {{7, 0, 2, true}, {6, 0, 0}, {7, -1, 0}},
       {{0, 1, T::TRIPLE}, {1, 2, T::SINGLE}}},
      {"valence5-cleanup3",
       {{7, 0, 3, true}, {7, 0, 0}},
       {{0, 1, T::DOUBLE}}},
      {"valence5-cleanup4",
       {{7, 0, 1, true}, {14, -1, 0}, {14, -1, 0}},
       {{0, 1, T::DOUBLE}, {0, 2, T::DOUBLE}}},
      {"valence5-cleanup5-oxygen",
       {{7, 0, 3, true}, {8, 0, 0}},
       {{0, 1, T::DOUBLE}}},
      {"valence5-cleanup5-sulfur",
       {{7, 0, 3, true}, {16, 0, 0}},
       {{0, 1, T::DOUBLE}}},
      {"valence5-cleanup5-fluorine",
       {{7, 0, 3, true}, {9, 0, 0}},
       {{0, 1, T::DOUBLE}}},
      {"valence5-cleanup5-chlorine",
       {{7, 0, 3, true}, {17, 0, 0}},
       {{0, 1, T::DOUBLE}}},
      {"valence5-cleanupb",
       {{7, 0, 3, true}, {6, 0, 0}},
       {{0, 1, T::DOUBLE}}},
      {"ordered-partial-error",
       {{7, 0, 3, true}, {8, 0, 0}, {6, 0, 0}},
       {{0, 1, T::DOUBLE}, {1, 2, T::THREECENTER}}},
      {"chlorine8-empty-neighbors", {{17, -1, 8}}, {}},
      {"chlorine8-all-oxygen",
       {{17, -1, 0}, {8, 0, 0}, {8, 0, 0}, {8, 0, 0}, {8, 0, 0}},
       {{0, 1, T::DOUBLE}, {0, 2, T::DOUBLE}, {0, 3, T::DOUBLE},
        {0, 4, T::DOUBLE}}},
      {"chlorine5-dispatch-helper-mismatch", {{17, 1, 5}}, {}},
      {"chlorine3-cleanup",
       {{17, 0, 0}, {16, 0, 0}},
       {{0, 1, T::TRIPLE, D::BEGINWEDGE}}},
      {"sulfur7-cleanup1",
       {{16, -1, 4}, {8, 4, 2}, {6, -5, 0}},
       {{0, 1, T::DOUBLE, D::BEGINWEDGE}, {0, 2, T::SINGLE}}},
      {"sulfur7-cleanup2",
       {{16, -1, 4}, {7, 0, 0}},
       {{0, 1, T::TRIPLE}}},
      {"sulfur7-cleanup3",
       {{16, -1, 5}, {7, 0, 3}},
       {{0, 1, T::DOUBLE}}},
      {"sulfur7-all-miss", {{16, -1, 7}}, {}},
      {"sulfur7-long-path-cleanup",
       {{16, -1, 5}, {6, 0, 0}, {6, 0, 0}, {7, 0, 2}},
       {{0, 1, T::DOUBLE}, {1, 2, T::SINGLE}, {2, 3, T::DOUBLE}}},
      {"sulfur8-source-defined-no-op",
       {{16, -1, 6}, {7, 0, 3}},
       {{0, 1, T::DOUBLE}}},
      {"bromine-selenium",
       {{35, 0, 0}, {34, 0, 0}},
       {{0, 1, T::TRIPLE, D::BEGINDASH}}},
      {"bromine-charged-miss",
       {{35, 1, 0}, {34, 0, 0}},
       {{0, 1, T::TRIPLE}}},
      {"bromine-degree-miss",
       {{35, 0, 0}, {34, 0, 0}, {6, 0, 0}},
       {{0, 1, T::DOUBLE}, {0, 2, T::SINGLE}}},
      {"bromine-element-miss",
       {{35, 0, 0}, {6, 0, 0}},
       {{0, 1, T::TRIPLE}}},
  };

  auto cleanup6 = make_valence5n_cleanup6_case(
      "valence5-cleanup6", 7, T::TRIPLE);
  cleanup6.atoms[2].is_aromatic = true;
  cases.push_back({cleanup6.case_id, cleanup6.atoms, cleanup6.bonds});

  auto cleanup7 = make_valence5n_cleanup7_case(
      "valence5-cleanup7", 7, T::TRIPLE);
  cleanup7.atoms[2].is_aromatic = true;
  cases.push_back({cleanup7.case_id, cleanup7.atoms, cleanup7.bonds});

  auto cleanup8 = make_valence5n_cleanup8_case("valence5-cleanup8", 7);
  cleanup8.atoms[5].is_aromatic = true;
  cases.push_back({cleanup8.case_id, cleanup8.atoms, cleanup8.bonds});

  auto cleanup9 = make_valence5n_cleanup9_case("valence5-cleanup9", 7);
  cleanup9.atoms[5].is_aromatic = true;
  cases.push_back({cleanup9.case_id, cleanup9.atoms, cleanup9.bonds});

  auto cleanupa =
      make_valence5n_cleanupa_path_case("valence5-cleanupa", 1);
  cleanupa.atoms[0].is_aromatic = true;
  cases.push_back({cleanupa.case_id, cleanupa.atoms, cleanupa.bonds});

  for (const auto &test_case : cases) {
    print_clean_up_case(test_case);
  }
  return 0;
}

struct InchiToMolCase {
  explicit InchiToMolCase(std::string id) : case_id(std::move(id)) {}

  std::string case_id;
  std::string inchi = "InChI=1S/test";
  ScriptedInchiOutput source_output;
  bool sanitize = false;
  bool remove_hs = false;
  std::vector<unsigned int> ranks;
  std::string fail_on;
};

static inchi_Atom make_inchi_atom(const char *element) {
  inchi_Atom atom{};
  std::strncpy(atom.elname, element, ATOM_EL_LEN - 1);
  return atom;
}

static void connect_inchi_atoms(std::vector<inchi_Atom> &atoms,
                                unsigned int first, unsigned int second,
                                signed char bond_type,
                                signed char first_stereo = 0,
                                signed char second_stereo = 0) {
  auto &left = atoms.at(first);
  auto &right = atoms.at(second);
  const auto left_offset = static_cast<unsigned int>(left.num_bonds++);
  const auto right_offset = static_cast<unsigned int>(right.num_bonds++);
  left.neighbor[left_offset] = static_cast<AT_NUM>(second);
  left.bond_type[left_offset] = bond_type;
  left.bond_stereo[left_offset] = first_stereo;
  right.neighbor[right_offset] = static_cast<AT_NUM>(first);
  right.bond_type[right_offset] = bond_type;
  right.bond_stereo[right_offset] = second_stereo;
}

static inchi_Stereo0D make_stereo(std::array<AT_NUM, 4> neighbors,
                                  AT_NUM central_atom, signed char type,
                                  signed char parity) {
  inchi_Stereo0D stereo{};
  std::copy(neighbors.begin(), neighbors.end(), stereo.neighbor);
  stereo.central_atom = central_atom;
  stereo.type = type;
  stereo.parity = parity;
  return stereo;
}

static void print_json_string(const std::string &value) {
  std::cout << '"';
  for (unsigned char byte : value) {
    switch (byte) {
      case '"':
        std::cout << "\\\"";
        break;
      case '\\':
        std::cout << "\\\\";
        break;
      case '\b':
        std::cout << "\\b";
        break;
      case '\f':
        std::cout << "\\f";
        break;
      case '\n':
        std::cout << "\\n";
        break;
      case '\r':
        std::cout << "\\r";
        break;
      case '\t':
        std::cout << "\\t";
        break;
      default:
        if (byte < 0x20) {
          const char *hex = "0123456789abcdef";
          std::cout << "\\u00" << hex[byte >> 4] << hex[byte & 0x0f];
        } else {
          std::cout << static_cast<char>(byte);
        }
    }
  }
  std::cout << '"';
}

static void print_string_array(const std::vector<std::string> &values) {
  std::cout << '[';
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (index != 0) std::cout << ',';
    print_json_string(values[index]);
  }
  std::cout << ']';
}

static void print_byte_array(const std::vector<unsigned char> &values) {
  std::cout << '[';
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (index != 0) std::cout << ',';
    std::cout << static_cast<unsigned int>(values[index]);
  }
  std::cout << ']';
}

static void print_source_atoms(const std::vector<inchi_Atom> &atoms) {
  std::cout << '[';
  for (std::size_t index = 0; index < atoms.size(); ++index) {
    if (index != 0) std::cout << ',';
    const auto &atom = atoms[index];
    std::cout << "{\"element\":";
    print_json_string(atom.elname);
    std::cout << ",\"num_bonds\":" << atom.num_bonds
              << ",\"neighbors\":[";
    for (int offset = 0; offset < atom.num_bonds; ++offset) {
      if (offset != 0) std::cout << ',';
      std::cout << atom.neighbor[offset];
    }
    std::cout << "],\"bond_types\":[";
    for (int offset = 0; offset < atom.num_bonds; ++offset) {
      if (offset != 0) std::cout << ',';
      std::cout << static_cast<int>(atom.bond_type[offset]);
    }
    std::cout << "],\"bond_stereo\":[";
    for (int offset = 0; offset < atom.num_bonds; ++offset) {
      if (offset != 0) std::cout << ',';
      std::cout << static_cast<int>(atom.bond_stereo[offset]);
    }
    std::cout << "],\"num_iso_h\":["
              << static_cast<int>(atom.num_iso_H[0]) << ','
              << static_cast<int>(atom.num_iso_H[1]) << ','
              << static_cast<int>(atom.num_iso_H[2]) << ','
              << static_cast<int>(atom.num_iso_H[3])
              << "],\"isotopic_mass\":" << atom.isotopic_mass
              << ",\"radical\":" << static_cast<int>(atom.radical)
              << ",\"charge\":" << static_cast<int>(atom.charge) << '}';
  }
  std::cout << ']';
}

static void print_source_stereo(
    const std::vector<inchi_Stereo0D> &stereo_entries) {
  std::cout << '[';
  for (std::size_t index = 0; index < stereo_entries.size(); ++index) {
    if (index != 0) std::cout << ',';
    const auto &stereo = stereo_entries[index];
    std::cout << "{\"neighbors\":[" << stereo.neighbor[0] << ','
              << stereo.neighbor[1] << ',' << stereo.neighbor[2] << ','
              << stereo.neighbor[3]
              << "],\"central_atom\":" << stereo.central_atom
              << ",\"type\":" << static_cast<int>(stereo.type)
              << ",\"parity\":" << static_cast<int>(stereo.parity) << '}';
  }
  std::cout << ']';
}

static void print_inchi_to_mol_atoms(const RDKit::RWMol &molecule) {
  std::cout << '[';
  for (std::size_t index = 0; index < molecule.atoms().size(); ++index) {
    if (index != 0) std::cout << ',';
    const auto &atom = molecule.atoms()[index];
    std::cout << "{\"index\":" << atom.getIdx()
              << ",\"atomic_number\":" << atom.getAtomicNum()
              << ",\"formal_charge\":" << atom.getFormalCharge()
              << ",\"num_explicit_hydrogens\":"
              << atom.getNumExplicitHs()
              << ",\"is_aromatic\":"
              << (atom.getIsAromatic() ? "true" : "false")
              << ",\"isotope\":" << atom.getIsotope()
              << ",\"num_radical_electrons\":"
              << atom.getNumRadicalElectrons()
              << ",\"no_implicit\":"
              << (atom.getNoImplicit() ? "true" : "false")
              << ",\"chiral_tag\":"
              << static_cast<int>(atom.getChiralTag()) << ",\"cip_rank\":";
    if (atom.hasCIPRank()) {
      std::cout << atom.getCIPRank();
    } else {
      std::cout << "null";
    }
    std::cout << '}';
  }
  std::cout << ']';
}

static void print_inchi_to_mol_bonds(const RDKit::RWMol &molecule) {
  std::cout << '[';
  for (std::size_t index = 0; index < molecule.bonds().size(); ++index) {
    if (index != 0) std::cout << ',';
    const auto &bond = molecule.bonds()[index];
    std::cout << "{\"index\":" << bond.getIdx()
              << ",\"begin_atom_index\":" << bond.getBeginAtomIdx()
              << ",\"end_atom_index\":" << bond.getEndAtomIdx()
              << ",\"bond_type\":" << static_cast<int>(bond.getBondType())
              << ",\"direction\":" << static_cast<int>(bond.getBondDir())
              << ",\"is_aromatic\":"
              << (bond.getIsAromatic() ? "true" : "false")
              << ",\"stereo\":" << static_cast<int>(bond.getStereo())
              << ",\"stereo_atoms\":";
    print_numbers(bond.getStereoAtoms());
    std::cout << '}';
  }
  std::cout << ']';
}

static void print_inchi_to_mol_case(const InchiToMolCase &test_case) {
  scripted_inchi_output = test_case.source_output;
  scripted_seen_inchi.clear();
  scripted_seen_options.clear();
  scripted_get_count = 0;
  scripted_free_count = 0;
  scripted_outstanding_outputs = 0;
  RDKit::inchi_to_mol_control = {};
  RDKit::inchi_to_mol_control.ranks = test_case.ranks;
  RDKit::inchi_to_mol_control.fail_on = test_case.fail_on;
  RDKit::rdWarningLog.str("");
  RDKit::rdWarningLog.clear();
  RDKit::rdErrorLog.str("");
  RDKit::rdErrorLog.clear();

  RDKit::ExtraInchiReturnValues return_values;
  return_values.returnCode = -77;
  return_values.messagePtr = "old-message";
  return_values.logPtr = "old-log";
  return_values.auxInfoPtr = "old-aux";
  RDKit::RWMol *molecule = nullptr;
  std::string status = "return";
  std::string exception_kind;
  std::string exception_message;
  try {
    molecule = RDKit::InchiToMol(test_case.inchi, return_values,
                                 test_case.sanitize, test_case.remove_hs);
  } catch (const RDKit::MolSanitizeException &error) {
    status = "exception";
    exception_kind = "MolSanitizeException";
    exception_message = error.what();
  } catch (const RDKit::InvariantException &error) {
    status = "exception";
    exception_kind = error.kind;
    exception_message = error.what();
  } catch (const RDKit::ValueErrorException &error) {
    status = "exception";
    exception_kind = "ValueErrorException";
    exception_message = error.what();
  } catch (const std::exception &error) {
    status = "exception";
    exception_kind = "std::exception";
    exception_message = error.what();
  }

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"3765efd5f4f8a855a2c2231c1d62a81e7d5880da76804f47f400ed8f5a6464a1\","
               "\"operation\":\"InchiToMol\",\"case_id\":";
  print_json_string(test_case.case_id);
  std::cout << ",\"input\":{\"inchi\":";
  print_json_string(test_case.inchi);
  std::cout << ",\"return_code\":" << test_case.source_output.return_code
            << ",\"atoms\":";
  print_source_atoms(test_case.source_output.atoms);
  std::cout << ",\"stereo0d\":";
  print_source_stereo(test_case.source_output.stereo0d);
  std::cout << ",\"message\":";
  if (test_case.source_output.has_message) {
    print_json_string(test_case.source_output.message);
  } else {
    std::cout << "null";
  }
  std::cout << ",\"log\":";
  if (test_case.source_output.has_log) {
    print_json_string(test_case.source_output.log);
  } else {
    std::cout << "null";
  }
  std::cout << ",\"sanitize\":"
            << (test_case.sanitize ? "true" : "false")
            << ",\"remove_hs\":"
            << (test_case.remove_hs ? "true" : "false")
            << ",\"ranks\":";
  print_numbers(test_case.ranks);
  std::cout << ",\"fail_on\":";
  if (test_case.fail_on.empty()) {
    std::cout << "null";
  } else {
    print_json_string(test_case.fail_on);
  }
  std::cout << ",\"throw_on_get\":"
            << (test_case.source_output.throw_on_get ? "true" : "false")
            << ",\"throw_on_free\":"
            << (test_case.source_output.throw_on_free ? "true" : "false")
            << "},\"output\":{\"status\":";
  print_json_string(status);
  std::cout << ",\"exception\":";
  if (status == "exception") {
    std::cout << "{\"kind\":";
    print_json_string(exception_kind);
    std::cout << ",\"message\":";
    print_json_string(exception_message);
    std::cout << '}';
  } else {
    std::cout << "null";
  }
  std::cout << ",\"return_values\":{\"return_code\":"
            << return_values.returnCode << ",\"message\":";
  print_json_string(return_values.messagePtr);
  std::cout << ",\"log\":";
  print_json_string(return_values.logPtr);
  std::cout << ",\"aux_info\":";
  print_json_string(return_values.auxInfoPtr);
  std::cout << "},\"molecule\":";
  if (molecule == nullptr) {
    std::cout << "null";
  } else {
    std::cout << "{\"atom_count\":" << molecule->getNumAtoms()
              << ",\"bond_count\":" << molecule->getNumBonds()
              << ",\"atom_fields\":";
    print_inchi_to_mol_atoms(*molecule);
    std::cout << ",\"bond_fields\":";
    print_inchi_to_mol_bonds(*molecule);
    std::cout << '}';
  }
  std::cout << ",\"warning_log\":";
  print_json_string(RDKit::rdWarningLog.str());
  std::cout << ",\"error_log\":";
  print_json_string(RDKit::rdErrorLog.str());
  std::cout << ",\"calls\":";
  print_string_array(RDKit::inchi_to_mol_control.calls);
  std::cout << ",\"get_count\":" << scripted_get_count
            << ",\"free_count\":" << scripted_free_count
            << ",\"outstanding_outputs\":" << scripted_outstanding_outputs
            << ",\"seen_inchi\":";
  print_json_string(scripted_seen_inchi);
  std::cout << ",\"seen_options\":";
  print_json_string(scripted_seen_options);
  std::cout << "}}\n";
  delete molecule;
}

static std::vector<inchi_Atom> make_double_stereo_graph() {
  std::vector<inchi_Atom> atoms(6);
  for (auto &atom : atoms) atom = make_inchi_atom("C");
  connect_inchi_atoms(atoms, 0, 1, INCHI_BOND_TYPE_DOUBLE);
  connect_inchi_atoms(atoms, 0, 2, INCHI_BOND_TYPE_SINGLE);
  connect_inchi_atoms(atoms, 0, 3, INCHI_BOND_TYPE_SINGLE);
  connect_inchi_atoms(atoms, 1, 4, INCHI_BOND_TYPE_SINGLE);
  connect_inchi_atoms(atoms, 1, 5, INCHI_BOND_TYPE_SINGLE);
  return atoms;
}

static int inchi_to_mol_records() {
  std::vector<InchiToMolCase> cases;

  InchiToMolCase failure{"failure-return"};
  failure.source_output.return_code = 2;
  failure.source_output.has_log = true;
  failure.source_output.log = "new-log";
  failure.sanitize = true;
  failure.remove_hs = true;
  cases.push_back(failure);

  InchiToMolCase mapping{"atom-bond-isotope-radical-mapping"};
  mapping.source_output.atoms.resize(10);
  for (auto &atom : mapping.source_output.atoms) atom = make_inchi_atom("C");
  mapping.source_output.atoms[0].isotopic_mass = ISOTOPIC_SHIFT_FLAG + 1;
  mapping.source_output.atoms[0].charge = -2;
  mapping.source_output.atoms[0].radical = 2;
  mapping.source_output.atoms[0].num_iso_H[0] = 2;
  mapping.source_output.atoms[0].num_iso_H[1] = 2;
  mapping.source_output.atoms[0].num_iso_H[2] = 3;
  mapping.source_output.atoms[0].num_iso_H[3] = 4;
  mapping.source_output.atoms[1].radical = 3;
  mapping.source_output.atoms[2].radical = 4;
  mapping.source_output.atoms[2].num_iso_H[2] = 1;
  mapping.source_output.atoms[3].num_iso_H[3] = 1;
  const std::vector<std::pair<int, int>> mapping_bonds = {
      {0, 0}, {1, 1}, {2, -6}, {3, 6}, {4, -1},
      {1, 4}, {2, 3}, {1, -4}, {1, 99}};
  for (std::size_t index = 0; index < mapping_bonds.size(); ++index) {
    connect_inchi_atoms(mapping.source_output.atoms, 0,
                        static_cast<unsigned int>(index + 1),
                        static_cast<signed char>(mapping_bonds[index].first),
                        static_cast<signed char>(mapping_bonds[index].second));
  }
  mapping.source_output.has_message = true;
  mapping.source_output.message = "message";
  mapping.source_output.has_log = true;
  mapping.source_output.log = "log";
  mapping.remove_hs = true;
  cases.push_back(mapping);

  InchiToMolCase illegal{"illegal-bond-type"};
  illegal.source_output.atoms = {make_inchi_atom("C"),
                                 make_inchi_atom("C")};
  connect_inchi_atoms(illegal.source_output.atoms, 0, 1, 5);
  illegal.sanitize = true;
  illegal.remove_hs = true;
  cases.push_back(illegal);

  for (const auto &[case_id, sanitize, remove_hs] :
       std::vector<std::tuple<std::string, bool, bool>>{
           {"simple-no-sanitize", false, false},
           {"simple-sanitize", true, false},
           {"simple-remove-hs", true, true}}) {
    InchiToMolCase simple{case_id};
    simple.source_output.atoms = {make_inchi_atom("C")};
    simple.sanitize = sanitize;
    simple.remove_hs = remove_hs;
    cases.push_back(simple);
  }

  InchiToMolCase warning{"warning-return"};
  warning.source_output.return_code = inchi_Ret_WARNING;
  warning.source_output.atoms = {make_inchi_atom("O")};
  warning.source_output.has_message = true;
  warning.source_output.message = "source warning";
  cases.push_back(warning);

  InchiToMolCase empty{"empty-success"};
  cases.push_back(empty);

  for (const auto &[case_id, parity, original_left, ranks] :
       std::vector<std::tuple<std::string, signed char, AT_NUM,
                              std::vector<unsigned int>>>{
           {"double-z", INCHI_PARITY_ODD, 2, {0, 0, 9, 3, 8, 2}},
           {"double-e", INCHI_PARITY_EVEN, 2, {0, 0, 9, 3, 8, 2}},
           {"double-switch-ez", INCHI_PARITY_ODD, 3, {0, 0, 9, 3, 8, 2}},
           {"double-unknown", INCHI_PARITY_UNKNOWN, 2, {0, 0, 9, 3, 8, 2}},
           {"double-gcc-rank-boundary", INCHI_PARITY_ODD, 3,
            {0, 0, std::numeric_limits<unsigned int>::max(), 7, 8, 2}}}) {
    InchiToMolCase stereo{case_id};
    stereo.source_output.atoms = make_double_stereo_graph();
    stereo.source_output.stereo0d = {make_stereo(
        {original_left, 0, 1, 4}, NO_ATOM, INCHI_StereoType_DoubleBond,
        parity)};
    stereo.ranks = ranks;
    cases.push_back(stereo);
  }

  for (const auto &[case_id, parity] :
       std::vector<std::pair<std::string, signed char>>{
           {"parity-none", INCHI_PARITY_NONE},
           {"parity-undefined", INCHI_PARITY_UNDEFINED}}) {
    InchiToMolCase skipped{case_id};
    skipped.source_output.atoms = make_double_stereo_graph();
    skipped.source_output.stereo0d = {make_stereo(
        {2, 0, 1, 4}, NO_ATOM, INCHI_StereoType_DoubleBond, parity)};
    cases.push_back(skipped);
  }

  InchiToMolCase absent{"extended-double-missing-bond"};
  absent.source_output.atoms = {make_inchi_atom("C"),
                                make_inchi_atom("C")};
  absent.source_output.stereo0d = {make_stereo(
      {0, 0, 1, 1}, NO_ATOM, INCHI_StereoType_DoubleBond,
      INCHI_PARITY_ODD)};
  cases.push_back(absent);

  InchiToMolCase insufficient{"double-missing-neighbor"};
  insufficient.source_output.atoms = {
      make_inchi_atom("C"), make_inchi_atom("C"), make_inchi_atom("C")};
  connect_inchi_atoms(insufficient.source_output.atoms, 0, 1,
                      INCHI_BOND_TYPE_DOUBLE);
  connect_inchi_atoms(insufficient.source_output.atoms, 0, 2,
                      INCHI_BOND_TYPE_SINGLE);
  insufficient.source_output.stereo0d = {make_stereo(
      {2, 0, 1, 2}, NO_ATOM, INCHI_StereoType_DoubleBond,
      INCHI_PARITY_ODD)};
  cases.push_back(insufficient);

  InchiToMolCase tetra{"tetrahedral-four-neighbor"};
  tetra.source_output.atoms.resize(5);
  for (auto &atom : tetra.source_output.atoms) atom = make_inchi_atom("C");
  for (unsigned int neighbor = 1; neighbor < 5; ++neighbor) {
    connect_inchi_atoms(tetra.source_output.atoms, 0, neighbor,
                        INCHI_BOND_TYPE_SINGLE);
  }
  tetra.source_output.stereo0d = {make_stereo(
      {1, 2, 3, 4}, 0, INCHI_StereoType_Tetrahedral, INCHI_PARITY_ODD)};
  cases.push_back(tetra);

  InchiToMolCase tetra3{"tetrahedral-three-neighbor"};
  tetra3.source_output.atoms.resize(4);
  for (auto &atom : tetra3.source_output.atoms) atom = make_inchi_atom("S");
  for (unsigned int neighbor = 1; neighbor < 4; ++neighbor) {
    connect_inchi_atoms(tetra3.source_output.atoms, 0, neighbor,
                        INCHI_BOND_TYPE_SINGLE);
  }
  tetra3.source_output.stereo0d = {make_stereo(
      {0, 1, 2, 3}, 0, INCHI_StereoType_Tetrahedral, INCHI_PARITY_ODD)};
  cases.push_back(tetra3);

  for (const auto &[case_id, type] :
       std::vector<std::pair<std::string, signed char>>{
           {"stereo-none", INCHI_StereoType_None},
           {"stereo-allene", INCHI_StereoType_Allene},
           {"stereo-unrecognized", 99}}) {
    InchiToMolCase stereo_type{case_id};
    stereo_type.source_output.atoms = {make_inchi_atom("C")};
    stereo_type.source_output.stereo0d = {
        make_stereo({0, 0, 0, 0}, 0, type, INCHI_PARITY_ODD)};
    cases.push_back(stereo_type);
  }

  InchiToMolCase direction_conflict{"assign-bond-dirs-conflict"};
  direction_conflict.source_output.atoms = make_double_stereo_graph();
  direction_conflict.source_output.atoms[0].bond_stereo[1] =
      INCHI_BOND_STEREO_SINGLE_1UP;
  direction_conflict.source_output.stereo0d = {make_stereo(
      {2, 0, 1, 4}, NO_ATOM, INCHI_StereoType_DoubleBond,
      INCHI_PARITY_ODD)};
  direction_conflict.ranks = {0, 0, 9, 3, 8, 2};
  cases.push_back(direction_conflict);

  for (const auto &[case_id, fail_on, sanitize, remove_hs] :
       std::vector<std::tuple<std::string, std::string, bool, bool>>{
           {"atomic-number-exception", "atomic_number", false, false},
           {"average-weight-exception", "average_atomic_weight", false, false},
           {"property-cache-exception", "update_property_cache", false, false},
           {"cip-rank-exception", "assign_atom_cip_ranks", false, false},
           {"sanitize-exception", "sanitize_molecule", true, false},
           {"remove-hs-exception", "remove_hydrogens", true, true},
           {"assign-stereo-exception", "assign_stereochemistry", false,
            false}}) {
    InchiToMolCase exception_case{case_id};
    exception_case.source_output.atoms = {make_inchi_atom("C")};
    if (fail_on == "assign_atom_cip_ranks") {
      exception_case.source_output.stereo0d = {make_stereo(
          {0, 0, 0, 0}, 0, INCHI_StereoType_None, INCHI_PARITY_ODD)};
    }
    exception_case.fail_on = fail_on;
    exception_case.sanitize = sanitize;
    exception_case.remove_hs = remove_hs;
    cases.push_back(exception_case);
  }

  InchiToMolCase invalid_element{"invalid-element-exception"};
  invalid_element.source_output.atoms = {make_inchi_atom("X")};
  cases.push_back(invalid_element);

  InchiToMolCase get_exception{"get-struct-exception"};
  get_exception.source_output.throw_on_get = true;
  cases.push_back(get_exception);

  InchiToMolCase free_exception{"free-struct-exception"};
  free_exception.source_output.atoms = {make_inchi_atom("C")};
  free_exception.source_output.throw_on_free = true;
  cases.push_back(free_exception);

  for (const auto &test_case : cases) print_inchi_to_mol_case(test_case);
  return 0;
}

struct FixOptionSymbolCase {
  std::string case_id;
  std::vector<unsigned char> input;
  std::size_t output_size;
};

static void print_fix_option_symbol_case(const FixOptionSymbolCase &test_case) {
  auto input = test_case.input;
  std::vector<unsigned char> output(test_case.output_size, 0xa5);
  const auto output_before = output;
  RDKit::fixOptionSymbol(reinterpret_cast<const char *>(input.data()),
                         reinterpret_cast<char *>(output.data()));

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"5770d6ab655210bd67e5351ef3ebe39aadc0c8b8e96c1b6e2ec0935acb8844d9\","
               "\"operation\":\"fixOptionSymbol\",\"case_id\":";
  print_json_string(test_case.case_id);
  std::cout << ",\"input\":{\"input_bytes\":";
  print_byte_array(test_case.input);
  std::cout << ",\"output_bytes\":";
  print_byte_array(output_before);
  std::cout << "},\"output\":{\"input_bytes\":";
  print_byte_array(input);
  std::cout << ",\"output_bytes\":";
  print_byte_array(output);
  std::cout << "}}\n";
}

static int fix_option_symbol_records() {
  std::vector<FixOptionSymbolCase> cases{
      {"empty", {0}, 3},
      {"slash-and-ordinary", {'/', 'A', '-', '/', 0}, 7},
      {"first-nul-stops", {'/', 'a', 0, '/', 'b', 0}, 8},
      {"exact-capacity", {'/', '-', 'x', 0}, 4},
  };
  std::vector<unsigned char> all_nonzero_bytes;
  for (unsigned int byte = 1; byte <= 255; ++byte) {
    all_nonzero_bytes.push_back(static_cast<unsigned char>(byte));
  }
  all_nonzero_bytes.push_back(0);
  cases.push_back(
      {"all-nonzero-bytes", std::move(all_nonzero_bytes), 260});

  for (const auto &test_case : cases) {
    print_fix_option_symbol_case(test_case);
  }
  return 0;
}

static void print_r_cleanup_atoms(const RDKit::RWMol &molecule) {
  std::cout << '[';
  for (std::size_t index = 0; index < molecule.atoms().size(); ++index) {
    if (index != 0) std::cout << ',';
    const auto &atom = molecule.atoms()[index];
    std::cout << "{\"index\":" << atom.getIdx()
              << ",\"atomic_number\":" << atom.getAtomicNum()
              << ",\"formal_charge\":" << atom.getFormalCharge()
              << ",\"num_explicit_hydrogens\":"
              << atom.getNumExplicitHs() << ",\"is_aromatic\":"
              << (atom.getIsAromatic() ? "true" : "false")
              << ",\"isotope\":" << atom.getIsotope()
              << ",\"num_radical_electrons\":"
              << atom.getNumRadicalElectrons() << ",\"no_implicit\":"
              << (atom.getNoImplicit() ? "true" : "false")
              << ",\"chiral_tag\":" << static_cast<int>(atom.getChiralTag())
              << ",\"cip_rank\":";
    if (atom.hasCIPRank()) {
      std::cout << atom.getCIPRank();
    } else {
      std::cout << "null";
    }
    std::cout << '}';
  }
  std::cout << ']';
}

static void print_r_cleanup_bonds(const RDKit::RWMol &molecule) {
  std::cout << '[';
  for (std::size_t index = 0; index < molecule.bonds().size(); ++index) {
    if (index != 0) std::cout << ',';
    const auto &bond = molecule.bonds()[index];
    std::cout << "{\"index\":" << bond.getIdx()
              << ",\"begin_atom_index\":" << bond.getBeginAtomIdx()
              << ",\"end_atom_index\":" << bond.getEndAtomIdx()
              << ",\"bond_type\":" << static_cast<int>(bond.getBondType())
              << ",\"direction\":" << static_cast<int>(bond.getBondDir())
              << ",\"is_aromatic\":"
              << (bond.getIsAromatic() ? "true" : "false")
              << ",\"stereo\":" << static_cast<int>(bond.getStereo())
              << ",\"stereo_atoms\":";
    std::cout << '[';
    for (std::size_t stereo_index = 0;
         stereo_index < bond.getStereoAtoms().size(); ++stereo_index) {
      if (stereo_index != 0) std::cout << ',';
      std::cout << bond.getStereoAtoms()[stereo_index];
    }
    std::cout << ']';
    std::cout << '}';
  }
  std::cout << ']';
}

struct RCleanUpCase {
  std::string case_id;
  std::vector<RDKit::GraphAtom> atoms;
  std::vector<RDKit::GraphBond> bonds;
  bool decorate = false;
};

static RCleanUpCase make_r_cleanup_star(std::string case_id,
                                        int center_charge,
                                        const std::vector<int> &charges,
                                        bool reverse_bonds = false) {
  RCleanUpCase result;
  result.case_id = std::move(case_id);
  result.atoms.push_back({8, charges.at(0), 0, false});
  result.atoms.push_back({17, center_charge, 0, false});
  for (std::size_t index = 1; index < charges.size(); ++index) {
    result.atoms.push_back({8, charges[index], 0, false});
  }
  for (std::size_t oxygen = 0; oxygen < charges.size(); ++oxygen) {
    const auto oxygen_index =
        oxygen == 0 ? 0U : static_cast<unsigned int>(oxygen + 1);
    const bool reverse = reverse_bonds && oxygen % 2 == 0;
    result.bonds.push_back(
        {reverse ? oxygen_index : 1U, reverse ? 1U : oxygen_index,
         RDKit::Bond::SINGLE, RDKit::Bond::NONE});
  }
  return result;
}

static void decorate_r_cleanup_molecule(RDKit::RWMol &molecule) {
  for (unsigned int index = 0; index < molecule.getNumAtoms(); ++index) {
    auto *atom = molecule.getAtomWithIdx(index);
    atom->setNumExplicitHs(index + 1);
    atom->setIsAromatic(index % 2 != 0);
    atom->setIsotope(18 + index);
    atom->setNumRadicalElectrons(index + 2);
    atom->setNoImplicit(index % 2 == 0);
    atom->setChiralTag(index % 2 == 0 ? RDKit::Atom::CHI_TETRAHEDRAL_CW
                                      : RDKit::Atom::CHI_TETRAHEDRAL_CCW);
    atom->setCIPRank(91 + index);
  }
  for (unsigned int index = 0; index < molecule.getNumBonds(); ++index) {
    auto *bond = molecule.getBondWithIdx(index);
    bond->setBondDir(static_cast<RDKit::Bond::BondDir>((index + 1) % 7));
    bond->setIsAromatic(index % 2 == 0);
    bond->setStereo(static_cast<RDKit::Bond::BondStereo>((index + 1) % 8));
    bond->getStereoAtoms() = {index % molecule.getNumAtoms(),
                              (index + 3) % molecule.getNumAtoms()};
  }
}

static void print_r_cleanup_case(const RCleanUpCase &test_case) {
  RDKit::RWMol molecule(test_case.atoms, test_case.bonds);
  if (test_case.decorate) decorate_r_cleanup_molecule(molecule);

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"81165fe80fd36e56b94f44fc7f062b0971bdd86d96fec1ebbdd947c954c64858\","
               "\"dependency_sha256\":{"
               "\"atom_cpp\":\"bfa51d29a6b18c4bfa373b2803ee79499b7ffbfa40d5838efb273bba2078961e\","
               "\"smiles_parse_cpp\":\"270b0457fbee294a9a1a1d7ada856538d4faec96fd782ffff9d2c5bb9ec0ca59\","
               "\"substruct_match_cpp\":\"a9c2829c484c98df285fa503870351d1ac8d947cfe46ad38f3cb6664f1c9ede2\","
               "\"substruct_match_h\":\"9c0a6d4c56f8fcbc8b4c3144ffd2893547780ebba5602858b2c5415beaaca1cc\","
               "\"substruct_utils_cpp\":\"76dbf84efc3b628f206de39b213373d55031e2762b116026c2fd05581af95c5f\"},"
               "\"operation\":\"rCleanUp\",\"case_id\":";
  print_json_string(test_case.case_id);
  std::cout << ",\"input\":{\"atoms\":";
  print_r_cleanup_atoms(molecule);
  std::cout << ",\"bonds\":";
  print_r_cleanup_bonds(molecule);
  std::cout << "},\"output\":{\"status\":\"returned\",\"atoms\":";
  RDKit::rCleanUp(molecule);
  print_r_cleanup_atoms(molecule);
  std::cout << ",\"bonds\":";
  print_r_cleanup_bonds(molecule);
  std::cout << "}}\n";
}

static int r_clean_up_records() {
  std::vector<RCleanUpCase> cases;
  cases.push_back(make_r_cleanup_star("wrong-center-charge", 2,
                                      {-1, -1, -1, 0}));
  cases.push_back(make_r_cleanup_star("insufficient-negative-oxygens", 3,
                                      {-1, -1, 0, 0}));
  cases.push_back(make_r_cleanup_star("arbitrary-nonmatching-nonzero", 3,
                                      {-2, -3, 1, 2}));
  cases.push_back(make_r_cleanup_star("arbitrary-nonmatching-one-neutral", 3,
                                      {-2, -3, 1, 0}));
  cases.push_back(
      make_r_cleanup_star("all-negative", 3, {-1, -1, -1, -1}));
  cases.push_back(
      make_r_cleanup_star("nonzero-wildcard", 3, {-1, -1, -1, 2}));
  for (std::size_t neutral = 0; neutral < 4; ++neutral) {
    std::vector<int> charges(4, -1);
    charges[neutral] = 0;
    cases.push_back(make_r_cleanup_star(
        "neutral-target-" + std::to_string(neutral), 3, charges));
  }
  cases.push_back(make_r_cleanup_star("reversed-bond-endpoints", 3,
                                      {-1, -1, -1, 0}, true));
  auto preserved =
      make_r_cleanup_star("preserved-complete-fields", 3, {-1, -1, -1, 0});
  preserved.decorate = true;
  cases.push_back(std::move(preserved));

  auto topology_miss =
      make_r_cleanup_star("topology-bond-type-miss", 3, {-1, -1, -1, 0});
  topology_miss.bonds[2].bond_type = RDKit::Bond::DOUBLE;
  cases.push_back(std::move(topology_miss));

  cases.push_back(make_r_cleanup_star("overlapping-precollected-matches", 3,
                                      {-1, -1, -1, 0, 0}));

  auto first = make_r_cleanup_star("two-disconnected-matches", 3,
                                   {-1, -1, -1, 0});
  const auto second =
      make_r_cleanup_star("unused", 3, {-1, -1, -1, 0}, true);
  const auto atom_offset = static_cast<unsigned int>(first.atoms.size());
  first.atoms.insert(first.atoms.end(), second.atoms.begin(),
                     second.atoms.end());
  for (auto bond : second.bonds) {
    bond.begin_atom_index += atom_offset;
    bond.end_atom_index += atom_offset;
    first.bonds.push_back(bond);
  }
  cases.push_back(std::move(first));

  for (const auto &test_case : cases) print_r_cleanup_case(test_case);
  return 0;
}

struct MolToInchiAtomSpec {
  MolToInchiAtomSpec() = default;
  MolToInchiAtomSpec(int atomic_number, int formal_charge = 0)
      : atomic_number(atomic_number), formal_charge(formal_charge) {}

  int atomic_number = 6;
  int formal_charge = 0;
  unsigned int explicit_hydrogens = 0;
  bool aromatic = false;
  unsigned int isotope = 0;
  unsigned int radical_electrons = 0;
  bool no_implicit = false;
  RDKit::Atom::ChiralType chiral_tag = RDKit::Atom::CHI_UNSPECIFIED;
  unsigned int total_hydrogens = 0;
  std::optional<unsigned int> total_degree;
};

struct MolToInchiBondSpec {
  MolToInchiBondSpec() = default;
  MolToInchiBondSpec(
      unsigned int begin_atom_index, unsigned int end_atom_index,
      RDKit::Bond::BondType bond_type,
      RDKit::Bond::BondDir direction = RDKit::Bond::NONE,
      bool aromatic = false,
      RDKit::Bond::BondStereo stereo = RDKit::Bond::STEREONONE,
      std::vector<unsigned int> stereo_atoms = {})
      : begin_atom_index(begin_atom_index),
        end_atom_index(end_atom_index),
        bond_type(bond_type),
        direction(direction),
        aromatic(aromatic),
        stereo(stereo),
        stereo_atoms(std::move(stereo_atoms)) {}

  unsigned int begin_atom_index = 0;
  unsigned int end_atom_index = 0;
  RDKit::Bond::BondType bond_type = RDKit::Bond::UNSPECIFIED;
  RDKit::Bond::BondDir direction = RDKit::Bond::NONE;
  bool aromatic = false;
  RDKit::Bond::BondStereo stereo = RDKit::Bond::STEREONONE;
  std::vector<unsigned int> stereo_atoms;
};

struct MolToInchiCase {
  std::string case_id;
  std::vector<MolToInchiAtomSpec> atoms;
  std::vector<MolToInchiBondSpec> bonds;
  std::vector<std::vector<RDKit::RDGeom::Point3D>> conformers;
  bool needs_update_property_cache = false;
  std::optional<std::vector<unsigned char>> options;
  ScriptedMolToInchiOutput scripted_output;
  std::string fail_on;
};

static MolToInchiCase make_mol_to_inchi_case(std::string case_id) {
  MolToInchiCase result;
  result.case_id = std::move(case_id);
  result.scripted_output.return_code = 17;
  result.scripted_output.has_inchi = true;
  result.scripted_output.has_message = true;
  result.scripted_output.has_log = true;
  result.scripted_output.has_aux_info = true;
  result.scripted_output.inchi = "InChI=1S/scripted";
  result.scripted_output.message = "new-message";
  result.scripted_output.log = "new-log";
  result.scripted_output.aux_info = "new-aux";
  return result;
}

static RDKit::RWMol make_mol_to_inchi_molecule(
    const MolToInchiCase &test_case) {
  std::vector<RDKit::GraphAtom> atoms;
  for (const auto &atom : test_case.atoms) {
    atoms.push_back({atom.atomic_number, atom.formal_charge,
                     atom.explicit_hydrogens, atom.aromatic});
  }
  std::vector<RDKit::GraphBond> bonds;
  for (const auto &bond : test_case.bonds) {
    bonds.push_back({bond.begin_atom_index, bond.end_atom_index, bond.bond_type,
                     bond.direction});
  }
  RDKit::RWMol molecule(atoms, bonds);
  for (unsigned int index = 0; index < test_case.atoms.size(); ++index) {
    const auto &spec = test_case.atoms[index];
    auto *atom = molecule.getAtomWithIdx(index);
    atom->setIsotope(spec.isotope);
    atom->setNumRadicalElectrons(spec.radical_electrons);
    atom->setNoImplicit(spec.no_implicit);
    atom->setChiralTag(spec.chiral_tag);
    atom->setTotalNumHs(spec.total_hydrogens);
    if (spec.total_degree) atom->setTotalDegree(*spec.total_degree);
  }
  for (unsigned int index = 0; index < test_case.bonds.size(); ++index) {
    const auto &spec = test_case.bonds[index];
    auto *bond = molecule.getBondWithIdx(index);
    bond->setIsAromatic(spec.aromatic);
    bond->setStereo(spec.stereo);
    bond->getStereoAtoms() = spec.stereo_atoms;
  }
  molecule.setConformers(test_case.conformers);
  molecule.setNeedsUpdatePropertyCache(test_case.needs_update_property_cache);
  return molecule;
}

static void print_optional_bytes(const std::optional<std::string> &value) {
  if (!value) {
    std::cout << "null";
    return;
  }
  print_byte_array(std::vector<unsigned char>(value->begin(), value->end()));
}

static void print_scripted_output(const ScriptedMolToInchiOutput &output) {
  std::cout << "{\"return_code\":" << output.return_code << ",\"inchi\":";
  print_optional_bytes(output.has_inchi
                           ? std::optional<std::string>(output.inchi)
                           : std::nullopt);
  std::cout << ",\"message\":";
  print_optional_bytes(output.has_message
                           ? std::optional<std::string>(output.message)
                           : std::nullopt);
  std::cout << ",\"log\":";
  print_optional_bytes(output.has_log ? std::optional<std::string>(output.log)
                                      : std::nullopt);
  std::cout << ",\"aux_info\":";
  print_optional_bytes(output.has_aux_info
                           ? std::optional<std::string>(output.aux_info)
                           : std::nullopt);
  std::cout << ",\"throw_on_get\":"
            << (output.throw_on_get ? "true" : "false")
            << ",\"throw_on_free\":"
            << (output.throw_on_free ? "true" : "false") << '}';
}

static void print_mol_to_inchi_atom_specs(
    const std::vector<MolToInchiAtomSpec> &atoms) {
  std::cout << '[';
  for (std::size_t index = 0; index < atoms.size(); ++index) {
    if (index != 0) std::cout << ',';
    const auto &atom = atoms[index];
    std::cout << "{\"atomic_number\":" << atom.atomic_number
              << ",\"formal_charge\":" << atom.formal_charge
              << ",\"explicit_hydrogens\":" << atom.explicit_hydrogens
              << ",\"aromatic\":" << (atom.aromatic ? "true" : "false")
              << ",\"isotope\":" << atom.isotope
              << ",\"radical_electrons\":" << atom.radical_electrons
              << ",\"no_implicit\":"
              << (atom.no_implicit ? "true" : "false")
              << ",\"chiral_tag\":" << static_cast<int>(atom.chiral_tag)
              << ",\"total_hydrogens\":" << atom.total_hydrogens
              << ",\"total_degree\":";
    if (atom.total_degree) {
      std::cout << *atom.total_degree;
    } else {
      std::cout << "null";
    }
    std::cout << '}';
  }
  std::cout << ']';
}

static void print_mol_to_inchi_bond_specs(
    const std::vector<MolToInchiBondSpec> &bonds) {
  std::cout << '[';
  for (std::size_t index = 0; index < bonds.size(); ++index) {
    if (index != 0) std::cout << ',';
    const auto &bond = bonds[index];
    std::cout << "{\"begin_atom_index\":" << bond.begin_atom_index
              << ",\"end_atom_index\":" << bond.end_atom_index
              << ",\"bond_type\":" << static_cast<int>(bond.bond_type)
              << ",\"direction\":" << static_cast<int>(bond.direction)
              << ",\"aromatic\":" << (bond.aromatic ? "true" : "false")
              << ",\"stereo\":" << static_cast<int>(bond.stereo)
              << ",\"stereo_atoms\":";
    print_numbers(bond.stereo_atoms);
    std::cout << '}';
  }
  std::cout << ']';
}

static void print_json_double(double value) {
  const auto flags = std::cout.flags();
  const auto precision = std::cout.precision();
  std::cout << std::showpoint
            << std::setprecision(std::numeric_limits<double>::max_digits10)
            << value;
  std::cout.flags(flags);
  std::cout.precision(precision);
}

static void print_conformers(
    const std::vector<std::vector<RDKit::RDGeom::Point3D>> &conformers) {
  std::cout << '[';
  for (std::size_t conformer_index = 0; conformer_index < conformers.size();
       ++conformer_index) {
    if (conformer_index != 0) std::cout << ',';
    std::cout << '[';
    for (std::size_t atom_index = 0;
         atom_index < conformers[conformer_index].size(); ++atom_index) {
      if (atom_index != 0) std::cout << ',';
      const auto &point = conformers[conformer_index][atom_index];
      std::cout << '[';
      print_json_double(point[0]);
      std::cout << ',';
      print_json_double(point[1]);
      std::cout << ',';
      print_json_double(point[2]);
      std::cout << ']';
    }
    std::cout << ']';
  }
  std::cout << ']';
}

static void print_generation_atoms(const std::vector<inchi_Atom> &atoms) {
  std::cout << '[';
  for (std::size_t index = 0; index < atoms.size(); ++index) {
    if (index != 0) std::cout << ',';
    const auto &atom = atoms[index];
    std::cout << "{\"x\":";
    print_json_double(atom.x);
    std::cout << ",\"y\":";
    print_json_double(atom.y);
    std::cout << ",\"z\":";
    print_json_double(atom.z);
    std::cout << ",\"element\":";
    print_json_string(atom.elname);
    std::cout << ",\"isotopic_mass\":" << atom.isotopic_mass
              << ",\"charge\":" << static_cast<int>(atom.charge)
              << ",\"radical\":" << static_cast<int>(atom.radical)
              << ",\"num_iso_h\":["
              << static_cast<int>(atom.num_iso_H[0]) << ','
              << static_cast<int>(atom.num_iso_H[1]) << ','
              << static_cast<int>(atom.num_iso_H[2]) << ','
              << static_cast<int>(atom.num_iso_H[3])
              << "],\"num_bonds\":" << atom.num_bonds
              << ",\"neighbors\":[";
    for (int offset = 0; offset < atom.num_bonds; ++offset) {
      if (offset != 0) std::cout << ',';
      std::cout << atom.neighbor[offset];
    }
    std::cout << "],\"bond_types\":[";
    for (int offset = 0; offset < atom.num_bonds; ++offset) {
      if (offset != 0) std::cout << ',';
      std::cout << static_cast<int>(atom.bond_type[offset]);
    }
    std::cout << "],\"bond_stereo\":[";
    for (int offset = 0; offset < atom.num_bonds; ++offset) {
      if (offset != 0) std::cout << ',';
      std::cout << static_cast<int>(atom.bond_stereo[offset]);
    }
    std::cout << "]}";
  }
  std::cout << ']';
}

static void print_return_values(const RDKit::ExtraInchiReturnValues &values) {
  std::cout << "{\"return_code\":" << values.returnCode
            << ",\"message\":";
  print_byte_array(
      std::vector<unsigned char>(values.messagePtr.begin(), values.messagePtr.end()));
  std::cout << ",\"log\":";
  print_byte_array(
      std::vector<unsigned char>(values.logPtr.begin(), values.logPtr.end()));
  std::cout << ",\"aux_info\":";
  print_byte_array(
      std::vector<unsigned char>(values.auxInfoPtr.begin(), values.auxInfoPtr.end()));
  std::cout << '}';
}

static void print_mol_to_inchi_case(const MolToInchiCase &test_case) {
  auto molecule = make_mol_to_inchi_molecule(test_case);
  scripted_mol_to_inchi_output = test_case.scripted_output;
  scripted_seen_atoms.clear();
  scripted_seen_stereo0d.clear();
  scripted_seen_null_options = true;
  scripted_seen_generation_options.clear();
  scripted_generation_get_count = 0;
  scripted_generation_free_count = 0;
  scripted_generation_outstanding_outputs = 0;
  scripted_generation_calls.clear();
  RDKit::mol_to_inchi_control = {};
  RDKit::mol_to_inchi_control.active = true;
  RDKit::mol_to_inchi_control.needs_update_property_cache =
      test_case.needs_update_property_cache;
  RDKit::mol_to_inchi_control.fail_on = test_case.fail_on;
  RDKit::rdWarningLog.str("");
  RDKit::rdWarningLog.clear();
  RDKit::rdErrorLog.str("");
  RDKit::rdErrorLog.clear();

  RDKit::ExtraInchiReturnValues values;
  values.returnCode = -77;
  values.messagePtr = "old-message";
  values.logPtr = "old-log";
  values.auxInfoPtr = "old-aux";
  std::string inchi;
  std::string status = "returned";
  std::string exception_kind;
  std::string exception_message;
  try {
    const char *options =
        test_case.options
            ? reinterpret_cast<const char *>(test_case.options->data())
            : nullptr;
    inchi = RDKit::MolToInchi(molecule, values, options);
  } catch (const RDKit::MolSanitizeException &error) {
    status = "exception";
    exception_kind = "MolSanitizeException";
    exception_message = error.what();
  } catch (const RDKit::InvariantException &error) {
    status = "exception";
    exception_kind = error.kind;
    exception_message = error.what();
  } catch (const std::runtime_error &error) {
    status = "exception";
    exception_kind = "runtime_error";
    exception_message = error.what();
  }

  std::cout << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
               "\"rdkit_version\":\"2026.03.1\","
               "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\","
               "\"source_fragment_sha256\":\"000d4eff4353e06aaa802057674fa24f91726872b0c0f9bb593d7a1937721ac2\","
               "\"operation\":\"MolToInchi\",\"case_id\":";
  print_json_string(test_case.case_id);
  std::cout << ",\"input\":{\"atoms\":";
  print_mol_to_inchi_atom_specs(test_case.atoms);
  std::cout << ",\"bonds\":";
  print_mol_to_inchi_bond_specs(test_case.bonds);
  std::cout << ",\"conformers\":";
  print_conformers(test_case.conformers);
  std::cout << ",\"needs_update_property_cache\":"
            << (test_case.needs_update_property_cache ? "true" : "false")
            << ",\"options\":";
  if (test_case.options) {
    print_byte_array(*test_case.options);
  } else {
    std::cout << "null";
  }
  std::cout << ",\"scripted_output\":";
  print_scripted_output(test_case.scripted_output);
  std::cout << ",\"fail_on\":";
  if (test_case.fail_on.empty()) {
    std::cout << "null";
  } else {
    print_json_string(test_case.fail_on);
  }
  std::cout << "},\"output\":{\"status\":";
  print_json_string(status);
  std::cout << ",\"exception_kind\":";
  if (exception_kind.empty()) {
    std::cout << "null";
  } else {
    print_json_string(exception_kind);
  }
  std::cout << ",\"exception_message\":";
  if (exception_message.empty()) {
    std::cout << "null";
  } else {
    print_json_string(exception_message);
  }
  std::cout << ",\"inchi\":";
  print_byte_array(std::vector<unsigned char>(inchi.begin(), inchi.end()));
  std::cout << ",\"return_values\":";
  print_return_values(values);
  std::cout << ",\"warning_text\":";
  print_json_string(RDKit::rdWarningLog.str());
  std::cout << ",\"error_text\":";
  print_json_string(RDKit::rdErrorLog.str());
  std::cout << ",\"toolkit_calls\":";
  print_string_array(RDKit::mol_to_inchi_control.calls);
  std::cout << ",\"generation_calls\":";
  print_string_array(scripted_generation_calls);
  std::cout << ",\"captured_atoms\":";
  print_generation_atoms(scripted_seen_atoms);
  std::cout << ",\"captured_stereo0d\":";
  print_source_stereo(scripted_seen_stereo0d);
  std::cout << ",\"captured_options\":";
  if (scripted_seen_null_options) {
    std::cout << "null";
  } else {
    print_byte_array(std::vector<unsigned char>(
        scripted_seen_generation_options.begin(),
        scripted_seen_generation_options.end()));
  }
  std::cout << ",\"get_count\":" << scripted_generation_get_count
            << ",\"free_count\":" << scripted_generation_free_count
            << ",\"outstanding_outputs\":"
            << scripted_generation_outstanding_outputs
            << ",\"original_atoms\":";
  print_inchi_to_mol_atoms(molecule);
  std::cout << ",\"original_bonds\":";
  print_inchi_to_mol_bonds(molecule);
  std::cout << ",\"original_conformers\":";
  print_conformers(test_case.conformers);
  std::cout << "}}\n";
}

static int mol_to_inchi_records() {
  std::vector<MolToInchiCase> cases;

  auto atom_fields = make_mol_to_inchi_case("atom-fields-options-and-cache");
  for (int atomic_number : {6, 7, 8, 9, 17, 35, 53, 15}) {
    MolToInchiAtomSpec atom;
    atom.atomic_number = atomic_number;
    atom.total_hydrogens = atom_fields.atoms.size() + 2;
    atom_fields.atoms.push_back(atom);
  }
  atom_fields.atoms[0].isotope = 13;
  atom_fields.atoms[0].formal_charge = 130;
  atom_fields.atoms[0].radical_electrons = 1;
  atom_fields.atoms[7].isotope = 31;
  atom_fields.atoms[7].formal_charge = -130;
  atom_fields.needs_update_property_cache = true;
  atom_fields.options = std::vector<unsigned char>{
      '/', 'A', 'u', 'x', 'N', 'o', 'n', 'e', ' ', '-', 'F', 'i', 'x',
      'e', 'd', 'H', 0, 'i', 'g', 'n', 'o', 'r', 'e', 'd', 0};
  atom_fields.scripted_output.has_message = false;
  atom_fields.scripted_output.has_aux_info = false;
  cases.push_back(std::move(atom_fields));

  auto cleanup = make_mol_to_inchi_case("r-clean-up-on-clone");
  cleanup.atoms = {{8, -1}, {17, 3}, {8, 0}, {8, -1}, {8, -1}};
  cleanup.bonds = {{0, 1, RDKit::Bond::SINGLE},
                   {1, 2, RDKit::Bond::SINGLE},
                   {1, 3, RDKit::Bond::SINGLE},
                   {1, 4, RDKit::Bond::SINGLE}};
  cases.push_back(std::move(cleanup));

  auto conformers = make_mol_to_inchi_case("first-conformer-null-output");
  conformers.atoms = {{1}};
  conformers.atoms[0].total_hydrogens = 3;
  conformers.conformers = {{{1.25, -2.5, 3.75}}, {{9.0, 9.0, 9.0}}};
  conformers.scripted_output = {};
  cases.push_back(std::move(conformers));

  for (unsigned int degree : {3U, 4U}) {
    for (auto tag : {RDKit::Atom::CHI_TETRAHEDRAL_CW,
                     RDKit::Atom::CHI_TETRAHEDRAL_CCW}) {
      auto tetra = make_mol_to_inchi_case(
          "tetra-degree-" + std::to_string(degree) + "-tag-" +
          std::to_string(static_cast<int>(tag)));
      tetra.atoms.resize(degree + 1);
      tetra.atoms[0].chiral_tag = tag;
      tetra.atoms[0].total_degree = degree;
      for (unsigned int neighbor = 1; neighbor <= degree; ++neighbor) {
        tetra.bonds.push_back({0, neighbor, RDKit::Bond::SINGLE});
      }
      cases.push_back(std::move(tetra));
    }
  }
  for (unsigned int total_degree : {2U, 5U}) {
    auto tetra = make_mol_to_inchi_case(
        "tetra-invalid-total-degree-" + std::to_string(total_degree));
    const unsigned int graph_degree = total_degree == 2 ? 2 : 4;
    tetra.atoms.resize(graph_degree + 1);
    tetra.atoms[0].chiral_tag = total_degree == 2
                                    ? RDKit::Atom::CHI_TETRAHEDRAL_CW
                                    : RDKit::Atom::CHI_TETRAHEDRAL_CCW;
    tetra.atoms[0].total_degree = total_degree;
    for (unsigned int neighbor = 1; neighbor <= graph_degree; ++neighbor) {
      tetra.bonds.push_back({0, neighbor, RDKit::Bond::SINGLE});
    }
    cases.push_back(std::move(tetra));
  }

  for (int bond_type = RDKit::Bond::UNSPECIFIED;
       bond_type <= RDKit::Bond::ZERO; ++bond_type) {
    auto bond_case = make_mol_to_inchi_case(
        "bond-type-" + std::to_string(bond_type));
    bond_case.atoms.resize(2);
    bond_case.bonds.push_back(
        {0, 1, static_cast<RDKit::Bond::BondType>(bond_type)});
    cases.push_back(std::move(bond_case));
  }

  for (int direction = RDKit::Bond::NONE; direction <= RDKit::Bond::UNKNOWN;
       ++direction) {
    for (bool reversed : {false, true}) {
      auto direction_case = make_mol_to_inchi_case(
          "bond-direction-" + std::to_string(direction) +
          (reversed ? "-reversed" : "-forward"));
      direction_case.atoms.resize(2);
      direction_case.bonds.push_back(
          {reversed ? 1U : 0U, reversed ? 0U : 1U, RDKit::Bond::SINGLE,
           static_cast<RDKit::Bond::BondDir>(direction)});
      cases.push_back(std::move(direction_case));
    }
  }

  for (int stereo = RDKit::Bond::STEREOZ;
       stereo <= RDKit::Bond::STEREOATROPCCW; ++stereo) {
    auto stereo_case = make_mol_to_inchi_case(
        "double-bond-stereo-" + std::to_string(stereo));
    stereo_case.atoms.resize(4);
    stereo_case.bonds = {
        {0, 1, RDKit::Bond::SINGLE},
        {1, 2, RDKit::Bond::DOUBLE, RDKit::Bond::NONE, false,
         static_cast<RDKit::Bond::BondStereo>(stereo), {0, 3}},
        {2, 3, RDKit::Bond::SINGLE}};
    cases.push_back(std::move(stereo_case));
  }
  auto short_stereo = make_mol_to_inchi_case("short-stereo-atom-list");
  short_stereo.atoms.resize(2);
  short_stereo.bonds = {{0, 1, RDKit::Bond::DOUBLE, RDKit::Bond::NONE, false,
                         RDKit::Bond::STEREOE, {0}}};
  cases.push_back(std::move(short_stereo));

  auto stereo_any = make_mol_to_inchi_case("stereo-any-coordinate-collapse");
  stereo_any.atoms.resize(2);
  stereo_any.bonds = {{1, 0, RDKit::Bond::DOUBLE, RDKit::Bond::NONE, false,
                       RDKit::Bond::STEREOANY, {0, 1}}};
  stereo_any.conformers = {{{1.0, 2.0, 3.0}, {7.0, 8.0, 9.0}}};
  cases.push_back(std::move(stereo_any));

  auto stereo_swap = make_mol_to_inchi_case("double-bond-stereo-swap");
  stereo_swap.atoms.resize(4);
  stereo_swap.bonds = {
      {1, 2, RDKit::Bond::DOUBLE, RDKit::Bond::NONE, false,
       RDKit::Bond::STEREOE, {3, 0}},
      {0, 1, RDKit::Bond::SINGLE}};
  cases.push_back(std::move(stereo_swap));

  auto maxval = make_mol_to_inchi_case("maxval-early-return");
  maxval.atoms.resize(MAXVAL + 2);
  for (unsigned int neighbor = 1; neighbor <= MAXVAL + 1; ++neighbor) {
    maxval.bonds.push_back({0, neighbor, RDKit::Bond::SINGLE});
  }
  cases.push_back(std::move(maxval));

  for (const std::string failure :
       {"needs_update_property_cache", "update_property_cache", "kekulize",
        "element_symbol", "atomic_weight", "total_num_hydrogens",
        "calc_implicit_valence", "total_degree"}) {
    auto error_case = make_mol_to_inchi_case("toolkit-exception-" + failure);
    error_case.atoms.resize(1);
    error_case.fail_on = failure;
    if (failure == "update_property_cache") {
      error_case.needs_update_property_cache = true;
    } else if (failure == "atomic_weight") {
      error_case.atoms[0].isotope = 13;
    } else if (failure == "total_num_hydrogens") {
      error_case.atoms[0].atomic_number = 15;
    } else if (failure == "calc_implicit_valence" ||
               failure == "total_degree") {
      error_case.atoms[0].chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CW;
      error_case.atoms[0].total_degree = 3;
    }
    cases.push_back(std::move(error_case));
  }

  auto get_exception = make_mol_to_inchi_case("get-inchi-exception");
  get_exception.scripted_output.throw_on_get = true;
  cases.push_back(std::move(get_exception));
  auto free_exception = make_mol_to_inchi_case("free-inchi-exception");
  free_exception.scripted_output.throw_on_free = true;
  cases.push_back(std::move(free_exception));

  for (const auto &test_case : cases) print_mol_to_inchi_case(test_case);
  return 0;
}

struct InchiToInchiKeyCase {
  std::string case_id;
  std::vector<unsigned char> inchi;
  int scripted_status;
  std::vector<unsigned char> scripted_key_buffer;
};

static std::vector<unsigned char> bytes_with_nul(const std::string &value) {
  std::vector<unsigned char> result(value.begin(), value.end());
  result.push_back(0);
  return result;
}

static void print_inchi_to_inchi_key_case(
    const InchiToInchiKeyCase &test_case) {
  scripted_inchi_key_status = test_case.scripted_status;
  scripted_inchi_key_buffer = test_case.scripted_key_buffer;
  scripted_seen_inchi_key_input.clear();
  scripted_seen_inchi_key_xtra1 = -1;
  scripted_seen_inchi_key_xtra2 = -1;
  scripted_inchi_key_call_count = 0;
  RDKit::rdErrorLog.str("");
  RDKit::rdErrorLog.clear();

  const std::string input(test_case.inchi.begin(), test_case.inchi.end());
  const std::string key = RDKit::InchiToInchiKey(input);

  std::cout
      << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\"," 
         "\"rdkit_version\":\"2026.03.1\"," 
         "\"source_sha256\":\"104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f\"," 
         "\"source_fragment_sha256\":\"b2d19205873e8444cef0ce660bbcf36bbf35cc528e77802b40b9b28ec15e67ae\"," 
         "\"operation\":\"InchiToInchiKey\",\"case_id\":";
  print_json_string(test_case.case_id);
  std::cout << ",\"input\":{\"inchi\":";
  print_byte_array(test_case.inchi);
  std::cout << ",\"scripted_status\":" << test_case.scripted_status
            << ",\"scripted_key_buffer\":";
  print_byte_array(test_case.scripted_key_buffer);
  std::cout << "},\"output\":{\"key\":";
  print_byte_array(std::vector<unsigned char>(key.begin(), key.end()));
  std::cout << ",\"error_text\":";
  print_json_string(RDKit::rdErrorLog.str());
  std::cout << ",\"captured_inchi\":";
  print_byte_array(std::vector<unsigned char>(
      scripted_seen_inchi_key_input.begin(),
      scripted_seen_inchi_key_input.end()));
  std::cout << ",\"captured_xtra1\":" << scripted_seen_inchi_key_xtra1
            << ",\"captured_xtra2\":" << scripted_seen_inchi_key_xtra2
            << ",\"call_count\":" << scripted_inchi_key_call_count
            << "}}\n";
}

static int inchi_to_inchi_key_records() {
  const std::string methane = "InChI=1S/CH4/h1H4";
  std::vector<InchiToInchiKeyCase> cases;
  cases.push_back({"success-full-key",
                   std::vector<unsigned char>(methane.begin(), methane.end()),
                   INCHIKEY_OK,
                   bytes_with_nul("VNWKTOKETHGBQD-UHFFFAOYSA-N")});
  cases.push_back({"success-first-nul",
                   {'I', 'n', 'C', 'h', 'I', '=', '1', 'S', '/', 'C', '2',
                    'H', '6'},
                   INCHIKEY_OK,
                   {'K', 'E', 'Y', 0, 'i', 'g', 'n', 'o', 'r', 'e', 'd'}});
  cases.push_back({"unknown-error", {'x'}, INCHIKEY_UNKNOWN_ERROR, {}});
  cases.push_back({"empty-input", {}, INCHIKEY_EMPTY_INPUT, {}});
  cases.push_back(
      {"invalid-prefix", {'n', 'o', 't', '-', 'i', 'n', 'c', 'h', 'i'},
       INCHIKEY_INVALID_INCHI_PREFIX, {}});
  cases.push_back({"not-enough-memory", {'x'}, INCHIKEY_NOT_ENOUGH_MEMORY,
                   {}});
  cases.push_back(
      {"invalid-inchi-embedded-nul",
       {'I', 'n', 'C', 'h', 'I', '=', '1', 'S', '/', 'C', 'H', '4', 0,
        'i', 'g', 'n', 'o', 'r', 'e', 'd'},
       INCHIKEY_INVALID_INCHI, {}});
  cases.push_back({"invalid-standard-inchi", {'I', 'n', 'C', 'h', 'I'},
                   INCHIKEY_INVALID_STD_INCHI, {}});
  cases.push_back({"unknown-status", {'x'}, 77, {}});

  for (const auto &test_case : cases) {
    print_inchi_to_inchi_key_case(test_case);
  }
  return 0;
}

struct MolToInchiKeyCase {
  MolToInchiCase generation;
  int scripted_key_status = INCHIKEY_OK;
  std::vector<unsigned char> scripted_key_buffer;
  bool throw_on_key = false;
};

static MolToInchiKeyCase make_mol_to_inchi_key_case(std::string case_id) {
  MolToInchiKeyCase result;
  result.generation = make_mol_to_inchi_case(std::move(case_id));
  result.scripted_key_buffer = bytes_with_nul("KEY");
  return result;
}

static void print_mol_to_inchi_key_case(
    const MolToInchiKeyCase &test_case) {
  auto molecule = make_mol_to_inchi_molecule(test_case.generation);
  scripted_mol_to_inchi_output = test_case.generation.scripted_output;
  scripted_seen_atoms.clear();
  scripted_seen_stereo0d.clear();
  scripted_seen_null_options = true;
  scripted_seen_generation_options.clear();
  scripted_generation_get_count = 0;
  scripted_generation_free_count = 0;
  scripted_generation_outstanding_outputs = 0;
  scripted_generation_calls.clear();
  scripted_inchi_key_status = test_case.scripted_key_status;
  scripted_inchi_key_buffer = test_case.scripted_key_buffer;
  scripted_inchi_key_throw = test_case.throw_on_key;
  scripted_seen_inchi_key_input.clear();
  scripted_seen_inchi_key_xtra1 = -1;
  scripted_seen_inchi_key_xtra2 = -1;
  scripted_inchi_key_call_count = 0;
  scripted_adapter_calls.clear();
  RDKit::mol_to_inchi_control = {};
  RDKit::mol_to_inchi_control.active = true;
  RDKit::mol_to_inchi_control.needs_update_property_cache =
      test_case.generation.needs_update_property_cache;
  RDKit::mol_to_inchi_control.fail_on = test_case.generation.fail_on;
  RDKit::rdWarningLog.str("");
  RDKit::rdWarningLog.clear();
  RDKit::rdErrorLog.str("");
  RDKit::rdErrorLog.clear();

  std::string key;
  std::string status = "returned";
  std::string exception_kind;
  std::string exception_message;
  try {
    const char *options =
        test_case.generation.options
            ? reinterpret_cast<const char *>(
                  test_case.generation.options->data())
            : nullptr;
    key = RDKit::MolToInchiKey(molecule, options);
  } catch (const RDKit::MolSanitizeException &error) {
    status = "exception";
    exception_kind = "MolSanitizeException";
    exception_message = error.what();
  } catch (const RDKit::InvariantException &error) {
    status = "exception";
    exception_kind = error.kind;
    exception_message = error.what();
  } catch (const std::runtime_error &error) {
    status = "exception";
    exception_kind = "runtime_error";
    exception_message = error.what();
  }

  std::cout
      << "{\"schema_version\":\"cosmolkit-inchi-rdkit-cpp-v1\","
         "\"rdkit_version\":\"2026.03.1\","
         "\"source_sha256\":\"27b4eaef1714869c42dbf2998807018a03389b4f9ce40438e843248ebfc3614e\","
         "\"source_fragment_sha256\":\"f6cc64a50170bb62e6c7a97325815f1399a6d3c74a9bc1b4a8bf69da90ccae51\","
         "\"dependency_fragment_sha256\":{"
         "\"MolToInchi\":\"000d4eff4353e06aaa802057674fa24f91726872b0c0f9bb593d7a1937721ac2\","
         "\"InchiToInchiKey\":\"b2d19205873e8444cef0ce660bbcf36bbf35cc528e77802b40b9b28ec15e67ae\"},"
         "\"operation\":\"MolToInchiKey\",\"case_id\":";
  print_json_string(test_case.generation.case_id);
  std::cout << ",\"input\":{\"atoms\":";
  print_mol_to_inchi_atom_specs(test_case.generation.atoms);
  std::cout << ",\"bonds\":";
  print_mol_to_inchi_bond_specs(test_case.generation.bonds);
  std::cout << ",\"conformers\":";
  print_conformers(test_case.generation.conformers);
  std::cout << ",\"needs_update_property_cache\":"
            << (test_case.generation.needs_update_property_cache ? "true"
                                                                 : "false")
            << ",\"options\":";
  if (test_case.generation.options) {
    print_byte_array(*test_case.generation.options);
  } else {
    std::cout << "null";
  }
  std::cout << ",\"scripted_generation_output\":";
  print_scripted_output(test_case.generation.scripted_output);
  std::cout << ",\"generation_fail_on\":";
  if (test_case.generation.fail_on.empty()) {
    std::cout << "null";
  } else {
    print_json_string(test_case.generation.fail_on);
  }
  std::cout << ",\"scripted_key_status\":"
            << test_case.scripted_key_status
            << ",\"scripted_key_buffer\":";
  print_byte_array(test_case.scripted_key_buffer);
  std::cout << ",\"throw_on_key\":"
            << (test_case.throw_on_key ? "true" : "false")
            << "},\"output\":{\"status\":";
  print_json_string(status);
  std::cout << ",\"exception_kind\":";
  if (exception_kind.empty()) {
    std::cout << "null";
  } else {
    print_json_string(exception_kind);
  }
  std::cout << ",\"exception_message\":";
  if (exception_message.empty()) {
    std::cout << "null";
  } else {
    print_json_string(exception_message);
  }
  std::cout << ",\"key\":";
  print_byte_array(std::vector<unsigned char>(key.begin(), key.end()));
  std::cout << ",\"warning_text\":";
  print_json_string(RDKit::rdWarningLog.str());
  std::cout << ",\"error_text\":";
  print_json_string(RDKit::rdErrorLog.str());
  std::cout << ",\"toolkit_calls\":";
  print_string_array(RDKit::mol_to_inchi_control.calls);
  std::cout << ",\"adapter_calls\":";
  print_string_array(scripted_adapter_calls);
  std::cout << ",\"captured_atoms\":";
  print_generation_atoms(scripted_seen_atoms);
  std::cout << ",\"captured_stereo0d\":";
  print_source_stereo(scripted_seen_stereo0d);
  std::cout << ",\"captured_options\":";
  if (scripted_seen_null_options) {
    std::cout << "null";
  } else {
    print_byte_array(std::vector<unsigned char>(
        scripted_seen_generation_options.begin(),
        scripted_seen_generation_options.end()));
  }
  std::cout << ",\"captured_key_inchi\":";
  print_byte_array(std::vector<unsigned char>(scripted_seen_inchi_key_input.begin(),
                                               scripted_seen_inchi_key_input.end()));
  std::cout << ",\"captured_key_xtra1\":" << scripted_seen_inchi_key_xtra1
            << ",\"captured_key_xtra2\":" << scripted_seen_inchi_key_xtra2
            << ",\"generation_get_count\":"
            << scripted_generation_get_count
            << ",\"generation_free_count\":"
            << scripted_generation_free_count
            << ",\"generation_outstanding_outputs\":"
            << scripted_generation_outstanding_outputs
            << ",\"key_call_count\":" << scripted_inchi_key_call_count
            << ",\"original_atoms\":";
  print_inchi_to_mol_atoms(molecule);
  std::cout << ",\"original_bonds\":";
  print_inchi_to_mol_bonds(molecule);
  std::cout << ",\"original_conformers\":";
  print_conformers(test_case.generation.conformers);
  std::cout << "}}\n";
}

static int mol_to_inchi_key_records() {
  std::vector<MolToInchiKeyCase> cases;

  auto success = make_mol_to_inchi_key_case("success-complete-forwarding");
  success.generation.atoms = {{6}, {8, -1}};
  success.generation.atoms[0].isotope = 13;
  success.generation.atoms[0].radical_electrons = 1;
  success.generation.atoms[0].total_hydrogens = 3;
  success.generation.atoms[0].chiral_tag = RDKit::Atom::CHI_TETRAHEDRAL_CW;
  success.generation.atoms[0].total_degree = 3;
  success.generation.bonds = {
      {0, 1, RDKit::Bond::DOUBLE, RDKit::Bond::BEGINWEDGE, false,
       RDKit::Bond::STEREOZ, {0, 1}}};
  success.generation.conformers = {{{1.25, -2.5, 3.75}, {4.0, 5.0, 6.0}}};
  success.generation.needs_update_property_cache = true;
  success.generation.options =
      std::vector<unsigned char>{'/', 'F', 'i', 'x', 'e', 'd', 'H', 0,
                                 'i', 'g', 'n', 'o', 'r', 'e', 'd', 0};
  success.generation.scripted_output.inchi = "InChI=1S/forwarded";
  success.scripted_key_buffer = {'K', 'E', 'Y', 0, 'i', 'g', 'n', 'o', 'r', 'e', 'd'};
  cases.push_back(std::move(success));

  auto ordered_diagnostics =
      make_mol_to_inchi_key_case("generation-warning-before-key-error");
  ordered_diagnostics.generation.atoms.resize(2);
  ordered_diagnostics.generation.bonds = {
      {0, 1, RDKit::Bond::QUADRUPLE}};
  ordered_diagnostics.scripted_key_status = INCHIKEY_INVALID_INCHI;
  ordered_diagnostics.scripted_key_buffer.clear();
  cases.push_back(std::move(ordered_diagnostics));

  auto empty = make_mol_to_inchi_key_case("empty-generation-output");
  empty.generation.scripted_output.has_inchi = false;
  empty.scripted_key_status = INCHIKEY_EMPTY_INPUT;
  empty.scripted_key_buffer.clear();
  cases.push_back(std::move(empty));

  auto return_code =
      make_mol_to_inchi_key_case("generation-return-code-does-not-short-circuit");
  return_code.generation.scripted_output.return_code = inchi_Ret_FATAL;
  return_code.generation.scripted_output.inchi = "InChI=1S/still-forwarded";
  cases.push_back(std::move(return_code));

  auto maxval = make_mol_to_inchi_key_case("maxval-empty-return-still-keys");
  maxval.generation.atoms.resize(MAXVAL + 2);
  for (unsigned int neighbor = 1; neighbor <= MAXVAL + 1; ++neighbor) {
    maxval.generation.bonds.push_back({0, neighbor, RDKit::Bond::SINGLE});
  }
  maxval.scripted_key_status = INCHIKEY_EMPTY_INPUT;
  maxval.scripted_key_buffer.clear();
  cases.push_back(std::move(maxval));

  auto toolkit = make_mol_to_inchi_key_case("toolkit-exception-short-circuit");
  toolkit.generation.atoms.resize(1);
  toolkit.generation.fail_on = "kekulize";
  cases.push_back(std::move(toolkit));

  auto get_exception = make_mol_to_inchi_key_case("get-exception-short-circuit");
  get_exception.generation.scripted_output.throw_on_get = true;
  cases.push_back(std::move(get_exception));

  auto free_exception = make_mol_to_inchi_key_case("free-exception-short-circuit");
  free_exception.generation.scripted_output.throw_on_free = true;
  cases.push_back(std::move(free_exception));

  auto key_exception = make_mol_to_inchi_key_case("key-exception-after-cleanup");
  key_exception.throw_on_key = true;
  cases.push_back(std::move(key_exception));

  for (const auto &test_case : cases) {
    print_mol_to_inchi_key_case(test_case);
  }
  return 0;
}

int main(int argc, char **argv) {
  if (argc == 2 && std::string(argv[1]) == "--mol-to-inchi-key-records") {
    return mol_to_inchi_key_records();
  }
  if (argc == 2 &&
      std::string(argv[1]) == "--inchi-to-inchi-key-records") {
    return inchi_to_inchi_key_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--mol-to-inchi-records") {
    return mol_to_inchi_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--r-clean-up-records") {
    return r_clean_up_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--fix-option-symbol-records") {
    return fix_option_symbol_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--inchi-to-mol-records") {
    return inchi_to_mol_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--assign-bond-dirs-records") {
    return assign_bond_dirs_records();
  }
  if (argc == 2 &&
      std::string(argv[1]) == "--find-alternating-bonds-records") {
    return find_alternating_bonds_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--neighboring-si-records") {
    return neighboring_si_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--valence4n-cleanup1-records") {
    return valence4n_cleanup1_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--valence4n-cleanup2-records") {
    return valence4n_cleanup2_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--valence5n-cleanup1-records") {
    return valence5n_cleanup1_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--valence5n-cleanup2-records") {
    return valence5n_cleanup2_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--valence5n-cleanup3-records") {
    return valence5n_cleanup3_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--valence5n-cleanup4-records") {
    return valence5n_cleanup4_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--valence5n-cleanup5-records") {
    return valence5n_cleanup5_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--valence5n-cleanup6-records") {
    return valence5n_cleanup6_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--valence5n-cleanup7-records") {
    return valence5n_cleanup7_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--valence5n-cleanup8-records") {
    return valence5n_cleanup8_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--valence5n-cleanup9-records") {
    return valence5n_cleanup9_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--valence5n-cleanupa-records") {
    return valence5n_cleanupa_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--valence5n-cleanupb-records") {
    return valence5n_cleanupb_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--valence7s-cleanup1-records") {
    return valence7s_cleanup1_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--valence7s-cleanup2-records") {
    return valence7s_cleanup2_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--valence7s-cleanup3-records") {
    return valence7s_cleanup3_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--valence8s-cleanup1-records") {
    return valence8s_cleanup1_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--valence8cl-cleanup1-records") {
    return valence8cl_cleanup1_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--valence5cl-cleanup1-records") {
    return valence5cl_cleanup1_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--valence3cl-cleanup1-records") {
    return valence3cl_cleanup1_records();
  }
  if (argc == 2 && std::string(argv[1]) == "--clean-up-records") {
    return clean_up_records();
  }
  std::cerr << "usage: rdkit-inchi-cpp-oracle "
               "--mol-to-inchi-key-records|--inchi-to-inchi-key-records|"
               "--mol-to-inchi-records|"
               "--fix-option-symbol-records|"
               "--assign-bond-dirs-records|"
               "--find-alternating-bonds-records|"
               "--neighboring-si-records|--valence4n-cleanup1-records|"
               "--valence4n-cleanup2-records|--valence5n-cleanup1-records|"
               "--valence5n-cleanup2-records|--valence5n-cleanup3-records|"
               "--valence5n-cleanup4-records|--valence5n-cleanup5-records|"
               "--valence5n-cleanup6-records|--valence5n-cleanup7-records|"
               "--valence5n-cleanup8-records|--valence5n-cleanup9-records|"
               "--valence5n-cleanupa-records|--valence5n-cleanupb-records|"
               "--valence7s-cleanup1-records|--valence7s-cleanup2-records|"
               "--valence7s-cleanup3-records|--valence8s-cleanup1-records|"
               "--valence8cl-cleanup1-records|--valence5cl-cleanup1-records|"
               "--valence3cl-cleanup1-records|--clean-up-records|"
               "--inchi-to-mol-records\n";
  return 64;
}
