#include <fstream>
#include <iostream>
#include "gemmi/mmread.hpp"
#include "gemmi/to_cif.hpp"
#include "gemmi/to_mmcif.hpp"
int main(int argc, char** argv) {
 if (argc != 6) return 2;
 try {
  auto structure=gemmi::read_structure_file(argv[1]);
  gemmi::MmcifOutputGroups groups(std::string(argv[3]) == "1");
  std::string flag=argv[4]; bool value=std::string(argv[5]) == "1";
  if (flag == "atoms") groups.atoms=value;
  else if (flag == "block_name") groups.block_name=value;
  else if (flag == "entry") groups.entry=value;
  else if (flag == "database_status") groups.database_status=value;
  else if (flag == "author") groups.author=value;
  else if (flag == "cell") groups.cell=value;
  else if (flag == "symmetry") groups.symmetry=value;
  else if (flag == "entity") groups.entity=value;
  else if (flag == "entity_poly") groups.entity_poly=value;
  else if (flag == "struct_ref") groups.struct_ref=value;
  else if (flag == "chem_comp") groups.chem_comp=value;
  else if (flag == "exptl") groups.exptl=value;
  else if (flag == "diffrn") groups.diffrn=value;
  else if (flag == "reflns") groups.reflns=value;
  else if (flag == "refine") groups.refine=value;
  else if (flag == "title_keywords") groups.title_keywords=value;
  else if (flag == "ncs") groups.ncs=value;
  else if (flag == "struct_asym") groups.struct_asym=value;
  else if (flag == "origx") groups.origx=value;
  else if (flag == "struct_conf") groups.struct_conf=value;
  else if (flag == "struct_sheet") groups.struct_sheet=value;
  else if (flag == "struct_biol") groups.struct_biol=value;
  else if (flag == "assembly") groups.assembly=value;
  else if (flag == "conn") groups.conn=value;
  else if (flag == "cis") groups.cis=value;
  else if (flag == "modres") groups.modres=value;
  else if (flag == "scale") groups.scale=value;
  else if (flag == "atom_type") groups.atom_type=value;
  else if (flag == "entity_poly_seq") groups.entity_poly_seq=value;
  else if (flag == "tls") groups.tls=value;
  else if (flag == "software") groups.software=value;
  else if (flag == "group_pdb") groups.group_pdb=value;
  else if (flag == "auth_all") groups.auth_all=value;
  else return 2;
  std::ofstream output(argv[2],std::ios::binary);
  gemmi::cif::write_cif_to_stream(output,gemmi::make_mmcif_document(structure,groups));
  return output ? 0 : 2;
 } catch(const std::exception& e) { std::cerr << e.what(); return 1; }
}
