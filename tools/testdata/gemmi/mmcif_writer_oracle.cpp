#include <fstream>
#include <iostream>
#include <string>

#include "gemmi/mmread.hpp"
#include "gemmi/to_cif.hpp"
#include "gemmi/to_mmcif.hpp"

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: mmcif_writer_oracle INPUT OUTPUT\n";
    return 2;
  }
  try {
    gemmi::Structure structure = gemmi::read_structure_file(argv[1]);
    gemmi::cif::Document document = gemmi::make_mmcif_document(structure);
    std::ofstream output(argv[2], std::ios::binary);
    if (!output) {
      std::cerr << "failed to open output: " << argv[2] << '\n';
      return 2;
    }
    gemmi::cif::write_cif_to_stream(output, document);
    if (!output) {
      std::cerr << "failed to write output: " << argv[2] << '\n';
      return 2;
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
