// Development-only pinned-Gemmi PDB coordinate oracle (BIO-PDB-WRITE).
// Modes (argv[1]): n0 n1 n2 r2 r3 pdb cif. argv[2]=input argv[3]=output.
// All outputs are byte-exact escaped lines: "len\t<escaped-bytes>".
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "gemmi/mmread.hpp"
#include "gemmi/sprintf.hpp"
#include "gemmi/to_pdb.hpp"
#include "gemmi/util.hpp"

namespace gemmi {
std::string bridge_serial_hybrid36(int serial);
std::string bridge_seq_id(int num, bool has_value, char icode);
std::string bridge_padded_name(const std::string& name, gemmi::El el);
int bridge_use_hetatm(const std::string& res_name, char het_flag,
                     gemmi::EntityType entity_type);
struct ChainWriteResult {
  std::string bytes;
  int final_serial;
};
ChainWriteResult bridge_write_chain_result(const gemmi::Chain& chain,
                                            bool ter_records,
                                            bool numbered_ter,
                                            bool ter_ignores_type,
                                            bool preserve_serial);
} // namespace gemmi

namespace {

std::string escape_bytes(const std::string& raw) {
  std::string out;
  for (unsigned char c : raw) {
    if (c == '\n') {
      out += "\\n";
    } else if (c == '\t') {
      out += "\\t";
    } else if (c == '\\') {
      out += "\\\\";
    } else if (c < 0x20 || c >= 0x7f) {
      char buf[5];
      std::snprintf(buf, sizeof buf, "\\x%02x", c);
      out += buf;
    } else {
      out += static_cast<char>(c);
    }
  }
  return out;
}

// Seven-record coordinate whitelist projection: keeps only ATOM, HETATM,
// ANISOU, TER, MODEL, ENDMDL, END lines, byte-for-byte.
std::string project_seven_records(const std::string& text) {
  std::string out;
  size_t pos = 0;
  while (pos < text.size()) {
    size_t eol = text.find('\n', pos);
    if (eol == std::string::npos)
      eol = text.size();
    std::string line = text.substr(pos, eol - pos);
    pos = eol + 1;
    std::string rec = line.substr(0, 6);
    // ROOT-REFERENCE-REVIEW fix 1: the pinned source writes the END record
    // as WRITE("%-80s", "END") — a PADDED 80-byte line. Compare the exact
    // six-byte record prefix so ALL original 80 bytes + LF are preserved.
    bool keep = rec == "ATOM  " || rec == "HETATM" || rec == "ANISOU" ||
                rec == "TER   " || rec == "MODEL " || rec == "ENDMDL" ||
                rec == "END   ";
    if (keep) {
      out += line;
      out += '\n';
    }
  }
  return out;
}

// ---- Literal detached fixtures C00-C11 (frozen in the packet annex) ----

gemmi::Atom atom(const std::string& name, gemmi::El el, double x, double y,
                 double z, char altloc = '\0', float occ = 1.0f,
                 float b = 20.0f, int charge = 0, int serial = 0) {
  gemmi::Atom a;
  a.name = name;
  a.element = el;
  a.pos = gemmi::Position(x, y, z);
  a.altloc = altloc;
  a.occ = occ;
  a.b_iso = b;
  a.charge = charge;
  a.serial = serial;
  return a;
}

gemmi::Residue residue(const std::string& name, int seq, char icode = ' ',
                       char het_flag = '\0',
                       gemmi::EntityType et = gemmi::EntityType::Unknown,
                       const std::string& segment = "") {
  gemmi::Residue r;
  r.name = name;
  r.seqid = gemmi::SeqId(gemmi::SeqId::OptionalNum(seq), icode);
  r.het_flag = het_flag;
  r.entity_type = et;
  r.segment = segment;
  return r;
}

gemmi::Chain chain_of(const std::string& name,
                      const std::vector<gemmi::Residue>& residues) {
  gemmi::Chain c;
  c.name = name;
  c.residues = residues;
  return c;
}

gemmi::Structure fixture(const std::string& id) {
  gemmi::Structure st;
  st.name = id;
  if (id == "C01" || id == "C02" || id == "C03" || id == "C04" ||
      id == "C05" || id == "C06" || id == "C07" || id == "C08" ||
      id == "C09" || id == "C10" || id == "C11") {
    gemmi::Model m1(1);
    if (id == "C02") {
      gemmi::Residue r = residue("ALA", 1, ' ', 'A', gemmi::EntityType::Polymer);
      r.atoms.push_back(atom("CA", gemmi::El::C, 1.0, 2.0, 3.0));
      m1.chains.push_back(chain_of("A", {r}));
    } else if (id == "C03") {
      gemmi::Residue poly = residue("ALA", 1, ' ', 'A', gemmi::EntityType::Polymer);
      poly.atoms.push_back(atom("CA", gemmi::El::C, 1.0, 2.0, 3.0));
      gemmi::Residue lig = residue("LIG", 2, ' ', 'H', gemmi::EntityType::NonPolymer);
      lig.atoms.push_back(atom("C1", gemmi::El::C, 4.0, 5.0, 6.0));
      gemmi::Residue wat = residue("HOH", 3, ' ', '\0', gemmi::EntityType::Water);
      wat.atoms.push_back(atom("O", gemmi::El::O, 7.0, 8.0, 9.0));
      m1.chains.push_back(chain_of("A", {poly}));
      m1.chains.push_back(chain_of("B", {lig}));
      m1.chains.push_back(chain_of("C", {wat}));
    } else if (id == "C04") {
      gemmi::Model m2(7);
      for (gemmi::Model* m : {&m1, &m2}) {
        gemmi::Residue r1 = residue("ALA", 1, ' ', 'A', gemmi::EntityType::Polymer);
        r1.atoms.push_back(atom("N", gemmi::El::N, 0.5, 0.5, 0.5, '\0', 1.0f, 20.0f, 0, 5));
        r1.atoms.push_back(atom("CA", gemmi::El::C, 1.5, 1.5, 1.5, '\0', 1.0f, 20.0f, 0, 99999));
        gemmi::Residue r2 = residue("GLY", 2, ' ', 'A', gemmi::EntityType::Polymer);
        r2.atoms.push_back(atom("N", gemmi::El::N, 2.5, 2.5, 2.5, '\0', 1.0f, 20.0f, 0, 3));
        m->chains.push_back(chain_of("A", {r1, r2}));
      }
      // C04 is TWO models with noncontiguous numbers (1 and 7): push both.
      // (The ordinal RF4 validation caught the prior single-model bug where
      // only m2 was pushed; native output was faithful to that wrong input.)
      st.models.push_back(m1);
      st.models.push_back(m2);
    } else if (id == "C05") {
      gemmi::Residue r = residue("ALA", 9999, 'A', 'A', gemmi::EntityType::Polymer);
      r.atoms.push_back(atom("CA", gemmi::El::C, 1.0, 2.0, 3.0, 'a'));
      r.atoms.push_back(atom("CA", gemmi::El::C, 1.1, 2.1, 3.1, ' '));
      gemmi::Residue neg = residue("GLY", -1, ' ', 'A', gemmi::EntityType::Polymer);
      neg.atoms.push_back(atom("N", gemmi::El::N, 4.0, 5.0, 6.0));
      gemmi::Residue big = residue("GLY", 10000, ' ', 'A', gemmi::EntityType::Polymer);
      big.atoms.push_back(atom("N", gemmi::El::N, 7.0, 8.0, 9.0));
      m1.chains.push_back(chain_of("A", {r, neg, big}));
    } else if (id == "C06") {
      gemmi::Residue ion = residue("SO4", 1, ' ', 'H', gemmi::EntityType::NonPolymer, "SEG1");
      ion.atoms.push_back(atom("S", gemmi::El::S, 1.0, 2.0, 3.0, '\0', 1.0f, 20.0f, -2));
      ion.atoms.push_back(atom("O1", gemmi::El::O, 1.5, 2.5, 3.5, '\0', 1.0f, 20.0f, -1, 100001));
      gemmi::Residue hyd = residue("HOH", 2, ' ', '\0', gemmi::EntityType::Water, "SEG2");
      hyd.atoms.push_back(atom("D", gemmi::El::D, 4.0, 5.0, 6.0));
      hyd.atoms.push_back(atom("H", gemmi::El::H, 4.1, 5.1, 6.1));
      m1.chains.push_back(chain_of("A", {ion, hyd}));
    } else if (id == "C07") {
      gemmi::Residue r = residue("ALA", 1, ' ', 'A', gemmi::EntityType::Polymer);
      gemmi::Atom a = atom("CA", gemmi::El::C, 1.0, 2.0, 3.0);
      a.aniso.u11 = 0.01;
      r.atoms.push_back(a);
      gemmi::Atom sz = atom("N", gemmi::El::N, 4.0, 5.0, 6.0);
      sz.aniso.u11 = -0.0;
      sz.aniso.u22 = 0.0;
      sz.aniso.u33 = -0.0;
      r.atoms.push_back(sz);
      m1.chains.push_back(chain_of("A", {r}));
    } else if (id == "C08") {
      gemmi::Residue r = residue("ALA", 1, ' ', 'A', gemmi::EntityType::Polymer);
      r.atoms.push_back(atom("CA", gemmi::El::C, -0.0004, -0.0006, 99.99949));
      r.atoms.push_back(atom("N", gemmi::El::N, 0.0, 0.0, 0.0, '\0', 0.999999f, 999.999f));
      r.atoms.push_back(atom("O", gemmi::El::O, 1.0, 2.0, 3.0, '\0', 1.5f, 1000.0f));
      m1.chains.push_back(chain_of("A", {r}));
    } else if (id == "C09") {
      gemmi::Residue r = residue("ALA", 1, ' ', 'A', gemmi::EntityType::Polymer);
      r.atoms.push_back(atom("CA", gemmi::El::C, 123456.7, -123456.7, 99999.4));
      m1.chains.push_back(chain_of("A", {r}));
    } else if (id == "C10") {
      gemmi::Residue empty1 = residue("ALA", 1, ' ', 'A', gemmi::EntityType::Polymer);
      gemmi::Residue poly = residue("GLY", 2, ' ', 'A', gemmi::EntityType::Polymer);
      poly.atoms.push_back(atom("CA", gemmi::El::C, 1.0, 2.0, 3.0));
      gemmi::Residue empty2 = residue("VAL", 3, ' ', 'A', gemmi::EntityType::Polymer);
      gemmi::Residue wat = residue("HOH", 4, ' ', '\0', gemmi::EntityType::Water);
      wat.atoms.push_back(atom("O", gemmi::El::O, 4.0, 5.0, 6.0));
      m1.chains.push_back(chain_of("A", {empty1, poly, empty2, wat}));
      gemmi::Chain empty_chain = chain_of("Z", {});
      m1.chains.push_back(empty_chain);
    } else if (id == "C11") {
      st.cell = gemmi::UnitCell(10.0, 20.0, 30.0, 90.0, 90.0, 90.0);
      st.spacegroup_hm = "P 21 21 21";
      st.info["_entry.id"] = "C11X";
      st.info["_struct.title"] = "metadata rich";
      gemmi::Residue r = residue("ALA", 1, ' ', 'A', gemmi::EntityType::Polymer);
      r.atoms.push_back(atom("CA", gemmi::El::C, 1.0, 2.0, 3.0));
      m1.chains.push_back(chain_of("A", {r}));
    } else { // C01
      m1.chains.push_back(chain_of("A", {}));
    }
    if (id != "C04")
      st.models.push_back(m1);
  }
  return st;
}

// Chain fixtures G00-G11 (dedicated native chains).
gemmi::Chain chain_fixture(const std::string& id) {
  if (id == "G00")
    return chain_of("", {});
  if (id == "G01")
    return chain_of("AB", {});
  gemmi::Structure st = fixture("C" + id.substr(1));
  const gemmi::Model& model = st.models.at(0);
  if (id == "G10")
    return model.chains.at(0);
  return model.chains.at(0);
}

} // namespace

int main(int argc, char** argv) {
  if (argc != 4) {
    std::cerr << "usage: pdb_coordinate_oracle MODE IN OUT\n";
    return 2;
  }
  std::string mode = argv[1];
  std::ifstream in(argv[2]);
  std::ofstream out(argv[3], std::ios::binary);
  if (!in || !out) {
    std::cerr << "failed to open files\n";
    return 2;
  }
  try {
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty())
        continue;
      std::istringstream fields(line);
      std::string kind;
      fields >> kind;
      if (mode == "n0") {
        // line: BITS PRECISION
        unsigned long long bits = std::stoull(kind, nullptr, 16);
        int prec;
        fields >> prec;
        double value;
        std::memcpy(&value, &bits, 8);
        char buf[512];
        int len = gemmi::snprintf_z(buf, sizeof buf, "%.*f", prec, value);
        out << len << '\t' << escape_bytes(std::string(buf, len > 0 ? len : 0))
            << '\n';
      } else if (mode == "n1") {
        if (kind == "s") {
          int serial;
          fields >> serial;
          out << gemmi::bridge_serial_hybrid36(serial) << '\n';
        } else if (kind == "q") {
          std::string num;
          std::string icode;
          fields >> num >> icode;
          bool has = num != "_";
          int value = has ? std::stoi(num) : 0;
          char code = icode.empty() ? ' ' : icode[0];
          out << gemmi::bridge_seq_id(value, has, code) << '\n';
        }
      } else if (mode == "n2") {
        int index;
        fields >> index;
        static const char* names[] = {"C",  "CA", "FE", "FE", "H",
                                      "HD11", "D",  "H",  "1HB", "ABCD",
                                      "",   " CA "};
        static const gemmi::El els[] = {
            gemmi::El::C, gemmi::El::C, gemmi::El::Fe, gemmi::El::Fe,
            gemmi::El::H, gemmi::El::H, gemmi::El::D, gemmi::El::H,
            gemmi::El::C, gemmi::El::C, gemmi::El::C, gemmi::El::C};
        if (kind == "p") {
          out << gemmi::bridge_padded_name(names[index], els[index])
              << '\n';
        } else if (kind == "h") {
          static const char* res_names[] = {"ALA", "MSE", "HOH", "DA",
                                            "XXX", "GLY", "LIG", "ALA"};
          static const char het_flags[] = {'\0', 'A', 'H'};
          static const gemmi::EntityType entity_types[] = {
              gemmi::EntityType::Polymer, gemmi::EntityType::Polymer,
              gemmi::EntityType::Water,   gemmi::EntityType::Polymer,
              gemmi::EntityType::Unknown, gemmi::EntityType::Branched,
              gemmi::EntityType::NonPolymer, gemmi::EntityType::Unknown};
          int res_index = index / 3;
          int het_index = index % 3;
          char het = het_flags[het_index];
          out << gemmi::bridge_use_hetatm(res_names[res_index], het,
                                               entity_types[res_index])
              << '\n';
        }
      } else if (mode == "r2") {
        // line: GID TER NUMBERED IGNORES PRESERVE
        std::string ter_s, numbered_s, ignores_s, preserve_s;
        fields >> ter_s >> numbered_s >> ignores_s >> preserve_s;
        gemmi::Chain c = chain_fixture(kind);
        gemmi::ChainWriteResult result = gemmi::bridge_write_chain_result(
            c, ter_s == "1", numbered_s == "1", ignores_s == "1",
            preserve_s == "1");
        out << result.bytes.size() << '\t' << escape_bytes(result.bytes)
            << '\t' << result.final_serial << '\n';
      } else if (mode == "r1") {
        // line: CID VARIANT PRESERVE  (variant 0=as-is,1=preserved-serial
        // sense encoded via PRESERVE, 2=frozen tensor u11=0.02 when absent)
        std::string variant_s, preserve_s;
        fields >> variant_s >> preserve_s;
        int variant = std::stoi(variant_s);
        gemmi::Structure st = fixture(kind);
        const gemmi::Model& model = st.models.at(0);
        const gemmi::Chain& src = model.chains.at(0);
        gemmi::Chain single;
        single.name = src.name;
        gemmi::Residue r = src.residues.at(0);
        gemmi::Atom selected = r.atoms.at(0);
        if (variant == 2 && !selected.aniso.nonzero())
          selected.aniso.u11 = 0.02;
        r.atoms.clear();
        r.atoms.push_back(selected);
        single.residues.push_back(r);
        gemmi::ChainWriteResult result = gemmi::bridge_write_chain_result(
            single, true, true, false, preserve_s == "1");
        out << result.bytes.size() << '\t' << escape_bytes(result.bytes)
            << '\t' << result.final_serial << '\n';
      } else if (mode == "r3") {
        // line: CID TER NUMBERED IGNORES PRESERVE END
        std::string ter_s, numbered_s, ignores_s, preserve_s, end_s;
        fields >> ter_s >> numbered_s >> ignores_s >> preserve_s >> end_s;
        gemmi::Structure st = fixture(kind);
        gemmi::PdbWriteOptions opt;
        opt.ter_records = ter_s == "1";
        opt.numbered_ter = numbered_s == "1";
        opt.ter_ignores_type = ignores_s == "1";
        opt.preserve_serial = preserve_s == "1";
        opt.end_record = end_s == "1";
        // Coordinate-only reference: minimal_file disables the non-coordinate
        // records the seven-record whitelist excludes; cryst1 stays source
        // default TRUE and is removed by projection, matching the declared
        // upfront whitelist policy (NOT a metadata-free claim).
        opt.minimal_file = true;
        opt.seqres_records = false;
        opt.ssbond_records = false;
        opt.link_records = false;
        opt.cispep_records = false;
        std::ostringstream os;
        gemmi::write_pdb(st, os, opt);
        std::string projected = project_seven_records(os.str());
        out << projected.size() << '\t' << escape_bytes(projected) << '\n';
      } else if (mode == "pdb" || mode == "cif") {
        // line: INPUTPATH TER NUMBERED IGNORES PRESERVE END
        std::string path, ter_s, numbered_s, ignores_s, preserve_s, end_s;
        fields >> path >> ter_s >> numbered_s >> ignores_s >> preserve_s >> end_s;
        gemmi::Structure st = gemmi::read_structure_file(
            path.c_str(), mode == "pdb" ? gemmi::CoorFormat::Pdb
                                        : gemmi::CoorFormat::Mmcif);
        gemmi::PdbWriteOptions opt;
        opt.ter_records = ter_s == "1";
        opt.numbered_ter = numbered_s == "1";
        opt.ter_ignores_type = ignores_s == "1";
        opt.preserve_serial = preserve_s == "1";
        opt.end_record = end_s == "1";
        opt.minimal_file = true;
        opt.seqres_records = false;
        opt.ssbond_records = false;
        opt.link_records = false;
        opt.cispep_records = false;
        std::ostringstream os;
        gemmi::write_pdb(st, os, opt);
        std::string projected = project_seven_records(os.str());
        out << projected.size() << '\t' << escape_bytes(projected) << '\n';
      }
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
