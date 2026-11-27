# Loaded into the upstream build so the oracle uses RDKit's own headers,
# compile definitions and library targets, rather than an installed wheel.
include_guard(GLOBAL)
set(ref_dir "${CMAKE_CURRENT_LIST_DIR}")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${ref_dir}/reference.cpp" "${ref_dir}/../../../testdata/reference/descriptor_foundations.json")
file(READ "${ref_dir}/../../../testdata/reference/descriptor_foundations.json" pins)
string(JSON wanted GET "${pins}" "${DESCRIPTOR_REFERENCE}")
execute_process(COMMAND git rev-parse HEAD WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE revision OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
execute_process(COMMAND git status --porcelain WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE dirty OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
if(NOT revision STREQUAL wanted OR NOT dirty STREQUAL "")
  message(FATAL_ERROR "The ${DESCRIPTOR_REFERENCE} oracle requires clean RDKit ${wanted}")
endif()

# RDKit embeds __FILE__ in diagnostics. Random private checkout paths must not
# give the same executable a different digest on each build.
add_compile_options("-ffile-prefix-map=${CMAKE_SOURCE_DIR}=rdkit" "-ffile-prefix-map=${CMAKE_BINARY_DIR}=rdkit-build")

function(add_descriptor_reference)
  file(SHA256 "${ref_dir}/reference.cpp" adapter_sha)
  file(SHA256 "${ref_dir}/rdkit.cmake" cmake_sha)
  add_executable(descriptor_oracle "${ref_dir}/reference.cpp")
  target_link_libraries(descriptor_oracle PRIVATE Descriptors_static SmilesParse_static Boost::headers Eigen3::Eigen)
  target_compile_definitions(descriptor_oracle PRIVATE RDK_HAS_EIGEN3
    CK_REFERENCE_KIND="${DESCRIPTOR_REFERENCE}" CK_RDKIT_REVISION="${revision}"
    CK_ADAPTER_SHA="${adapter_sha}" CK_CMAKE_SHA="${cmake_sha}"
    CK_CXX="${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}"
    CK_BUILD="${CMAKE_BUILD_TYPE}" CK_CXX_FLAGS="${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
  if(DESCRIPTOR_REFERENCE STREQUAL "exact-ipc")
    target_compile_definitions(descriptor_oracle PRIVATE CK_EXACT_IPC)
    set_property(SOURCE "${CMAKE_SOURCE_DIR}/Code/GraphMol/Descriptors/GraphDescriptors.cpp"
      TARGET_DIRECTORY Descriptors_static APPEND PROPERTY COMPILE_OPTIONS -ffp-contract=off)
  endif()
  set_target_properties(descriptor_oracle PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
endfunction()
cmake_language(DEFER CALL add_descriptor_reference)
