project(SparseSolvTest CXX)

add_executable(
  SRLfemTest
  SRLfemTest.cpp
)

target_include_directories(
    SRLfemTest PUBLIC ${SparseSolvBase_SOURCE_DIR})

target_link_libraries(
  SRLfemTest
  GTest::gtest_main
  SparseSolv
)

include(GoogleTest)
gtest_discover_tests(SRLfemTest)
