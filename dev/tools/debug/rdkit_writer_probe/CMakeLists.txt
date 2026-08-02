cmake_minimum_required(VERSION 3.20)
project(rdkit_writer_probe LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(RDKit REQUIRED CONFIG PATHS ../../../../rdkit/build NO_DEFAULT_PATH)

add_executable(rdkit_writer_probe main.cpp)
target_link_libraries(
  rdkit_writer_probe
  PRIVATE
    RDKit::SmilesParse
    RDKit::GraphMol
    RDKit::RDGeneral
)
