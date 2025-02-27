###############################################################################
# Testing
###############################################################################
option(ENABLE_TESTING "Enable testing" OFF)
if(ENABLE_TESTING)
  find_program(VALGRIND_BINARY NAMES valgrind)
  # generate custom suppression file
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/lammps.supp "\n")
  file(GLOB VALGRIND_SUPPRESSION_FILES ${LAMMPS_TOOLS_DIR}/valgrind/[^.]*.supp)
  foreach(SUPP ${VALGRIND_SUPPRESSION_FILES})
    file(READ ${SUPP} SUPPRESSIONS)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/lammps.supp "${SUPPRESSIONS}")
  endforeach()
  set(VALGRIND_DEFAULT_OPTIONS "--leak-check=full --show-leak-kinds=all --track-origins=yes --suppressions=${CMAKE_BINARY_DIR}/lammps.supp")

  set(MEMORYCHECK_COMMAND "${VALGRIND_BINARY}" CACHE FILEPATH "Memory Check Command")
  set(MEMORYCHECK_COMMAND_OPTIONS "${VALGRIND_DEFAULT_OPTIONS}" CACHE STRING "Memory Check Command Options")

  # check if a faster linker is available.
  # only verified with GNU and Clang compilers and new CMake
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    if((${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        OR (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
      include(CheckCXXCompilerFlag)
      set(CMAKE_CUSTOM_LINKER_DEFAULT default)
      check_cxx_compiler_flag(-fuse-ld=lld HAVE_LLD_LINKER_FLAG)
      check_cxx_compiler_flag(-fuse-ld=gold HAVE_GOLD_LINKER_FLAG)
      check_cxx_compiler_flag(-fuse-ld=bfd HAVE_BFD_LINKER_FLAG)
      find_program(HAVE_LLD_LINKER_BIN lld ld.lld)
      find_program(HAVE_GOLD_LINKER_BIN ld.gold)
      find_program(HAVE_BFD_LINKER_BIN ld.bfd)
      if(HAVE_LLD_LINKER_FLAG AND HAVE_LLD_LINKER_BIN)
        set(CMAKE_CUSTOM_LINKER_DEFAULT lld)
      elseif(HAVE_BFD_LINKER_FLAG AND HAVE_BFD_LINKER_BIN)
        set(CMAKE_CUSTOM_LINKER_DEFAULT bfd)
      elseif(HAVE_GOLD_LINKER_FLAG AND HAVE_GOLD_LINKER_BIN)
        set(CMAKE_CUSTOM_LINKER_DEFAULT gold)
      endif()
      set(CMAKE_CUSTOM_LINKER_VALUES lld gold bfd default)
      set(CMAKE_CUSTOM_LINKER ${CMAKE_CUSTOM_LINKER_DEFAULT} CACHE STRING "Choose a custom linker for faster linking (lld, gold, bfd, default)")
      validate_option(CMAKE_CUSTOM_LINKER CMAKE_CUSTOM_LINKER_VALUES)
      mark_as_advanced(CMAKE_CUSTOM_LINKER)
      if(NOT "${CMAKE_CUSTOM_LINKER}" STREQUAL "default")
        target_link_options(lammps PUBLIC -fuse-ld=${CMAKE_CUSTOM_LINKER})
      endif()
    endif()
  endif()

  include(CTest)

  enable_testing()
  get_filename_component(LAMMPS_UNITTEST_DIR ${LAMMPS_SOURCE_DIR}/../unittest ABSOLUTE)
  get_filename_component(LAMMPS_UNITTEST_BIN ${CMAKE_BINARY_DIR}/unittest ABSOLUTE)
  add_subdirectory(${LAMMPS_UNITTEST_DIR} ${LAMMPS_UNITTEST_BIN})
endif()
