if (NOT DEBUG)
  string(APPEND CFLAGS " -O1")
endif()
if (NOT DEBUG)
  string(APPEND FFLAGS " -O1")
endif()
if (COMP_NAME STREQUAL eam)
  string(APPEND FFLAGS " -vector0")
endif()
set(CXX_LINKER "CXX")
set(KOKKOS_OPTIONS "--gcc-toolchain=/opt/gcc/9.3.0/snos --with-serial")
set(CXXFLAGS "-std=c++14")
if (compile_threaded)
  string(APPEND CXXFLAGS " -fopenmp")
endif()
if (DEBUG)
  string(APPEND CXXFLAGS " -g -Wall")
endif()
if (NOT DEBUG)
  string(APPEND CXXFLAGS " -O1")
endif()