// Development-only probe TU for BIO-PDB-WRITE. Includes the pinned
// third_party/gemmi/src/to_pdb.cpp UNCHANGED in its own translation unit
// so its anonymous-namespace helpers (encode_serial_in_hybrid36,
// write_seq_id, write_chain_atoms) are reachable without modifying
// upstream source. libgemmi_cpp.a's own to_pdb.o is not pulled because
// this TU already satisfies the referenced gemmi::write_pdb symbols. No
// algorithm is copied.
#include "../../../third_party/gemmi/src/to_pdb.cpp"

namespace gemmi {

// Bridges to the anonymous-namespace helpers of the included source.
std::string bridge_serial_hybrid36(int serial) {
  std::array<char, 8> str = encode_serial_in_hybrid36(serial);
  return std::string(str.data());
}

std::string bridge_seq_id(int num, bool has_value, char icode) {
  SeqId seqid;
  seqid.num = has_value ? SeqId::OptionalNum(num) : SeqId::OptionalNum();
  seqid.icode = icode;
  std::array<char, 8> str = write_seq_id(seqid);
  return std::string(str.data());
}

std::string bridge_padded_name(const std::string& name, El el) {
  Atom atom;
  atom.name = name;
  atom.element = el;
  return atom.padded_name();
}

int bridge_use_hetatm(const std::string& res_name, char het_flag,
                     EntityType entity_type) {
  Residue res;
  res.name = res_name;
  res.het_flag = het_flag;
  res.entity_type = entity_type;
  return use_hetatm(res) ? 1 : 0;
}

struct ChainWriteResult {
  std::string bytes;
  int final_serial;
};

// ONE source write_chain_atoms invocation returning BOTH the bytes and the
// final serial (ROOT-REFERENCE-REVIEW fix 2: two separate invocations
// doubled the actual source call census per row).
static ChainWriteResult run_chain(const Chain& chain, bool ter_records,
                                  bool numbered_ter, bool ter_ignores_type,
                                  bool preserve_serial) {
  PdbWriteOptions opt;
  opt.ter_records = ter_records;
  opt.numbered_ter = numbered_ter;
  opt.ter_ignores_type = ter_ignores_type;
  opt.preserve_serial = preserve_serial;
  std::ostringstream os;
  int serial = 0;
  write_chain_atoms(chain, os, serial, opt);
  return {os.str(), serial};
}

ChainWriteResult bridge_write_chain_result(const Chain& chain,
                                            bool ter_records,
                                            bool numbered_ter,
                                            bool ter_ignores_type,
                                            bool preserve_serial) {
  return run_chain(chain, ter_records, numbered_ter, ter_ignores_type,
                   preserve_serial);
}

} // namespace gemmi
