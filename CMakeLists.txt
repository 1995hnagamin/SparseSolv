# CMake version set
cmake_minimum_required(VERSION 3.21)
# determine compiler
#set(CMAKE_CXX_COMPILER icpc)
#set(CMAKE_CXX_COMPILER icpx)


#
if(CMAKE_CXX_COMPILER MATCHES "icpc")
 message("Intel classical Compiler use !")
elseif(CMAKE_CXX_COMPILER MATCHES "icpx")
 message("Intel new Compiler use !")
 set(CMAKE_AR "/opt/intel/oneapi/compiler/latest/bin/compiler/llvm-ar")
 set(CMAKE_NM "/opt/intel/oneapi/compiler/latest/bin/compiler/llvm-ar")
#set(CMAKE_RANLIB "/opt/intel/oneapi/compiler/latest/bin/compiler/llvm-ranlib")
endif()
#
# set true if you use pybind11
set(PYBIND_EXPORT 1)

# FetchContent setup
include(FetchContent)
#set(FETCHCONTENT_BASE_DIR ${SparseSolvBase_BINARY_DIR}/000_thirdparty)

# set base dir
project(SparseSolvBase CXX)

# set dir set
add_subdirectory(SparseSolv)
if(PYBIND_EXPORT)
 add_subdirectory(SparseSolvPy)
endif()

# sample 
add_subdirectory(examples/VoxelFEM)

# test setup
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/52eb8108c5bdec04579160ae17225d66034bd723.zip
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(googletest)
enable_testing()
add_subdirectory(tests)
