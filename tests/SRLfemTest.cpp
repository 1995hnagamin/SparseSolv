#include <gtest/gtest.h>
#include "SparseSolv/SparseMat.hpp"

TEST(SRLfemTest, IsFixedAssertions) {
  constexpr int size = 3;
  double matA[size][size] ={
    { 1.1, 1.2, 0.0 },
    { 1.2, 2.2, 0.0 },
    { 0.0, 0.0, 3.3 },
  };

  SRLfem::SparseMat matAs(size);
  for(int i = 0 ; i < size ; i++){
    for(int j = 0 ; j < size ; j++){
      if( fabs(matA[i][j]) > 1.0e-12 ){
        matAs.add(i, j, matA[i][j]);
      }
    }
  }
  EXPECT_FALSE(matAs.isFixed());
  matAs.fix();
  EXPECT_TRUE(matAs.isFixed());
}
