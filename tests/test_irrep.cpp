// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <gtest/gtest.h>
#include <vector>
#include <stdexcept>

#include "sun/irrep/irrep.h"
#include "testutils.h"


TEST(IrrepConstruction, DefaultConstructor)
{
    sun::Irrep alpha;
    EXPECT_EQ(alpha.n(), 0u);
    EXPECT_EQ(alpha.nrows(), 0u);
    EXPECT_EQ(alpha.ncols(), 0u);
    EXPECT_EQ(alpha.size(), 0u);
}


TEST(IrrepVectorConstruction, ValidConstruction)
{
    std::vector<unsigned int> v = {4, 3, 2};
    sun::Irrep alpha(v);
    EXPECT_EQ(alpha.n(), 9u);
    EXPECT_EQ(alpha.nrows(), 3u);
    EXPECT_EQ(alpha.ncols(), 4u);
    EXPECT_EQ(alpha.size(), 3u);
    expect_equal_vec(alpha.get_vector(), v);
}


TEST(IrrepVectorConstruction, NegativeValueThrowsError)
{
	std::vector<int> v = {4, 2, -1};
	
	EXPECT_THROW(sun::Irrep alpha(v), std::invalid_argument);
}


TEST(IrrepVectorConstruction, NonDescendingThrowsError)
{
	std::vector<int> v({4, 4, 3, 4});
	
	EXPECT_THROW(sun::Irrep alpha(v), std::invalid_argument);
}


TEST(IrrepInitializerListConstruction, ValidConstruction)
{
    sun::Irrep alpha({4, 3, 2});
    EXPECT_EQ(alpha.n(), 9u);
    EXPECT_EQ(alpha.nrows(), 3u);
    EXPECT_EQ(alpha.ncols(), 4u);
    EXPECT_EQ(alpha.size(), 3u);
}


TEST(IrrepInitializerListConstruction, NegativeValueThrowsError)
{
	EXPECT_THROW(sun::Irrep alpha({5, 4, -1, 3}), std::invalid_argument);
}


TEST(IrrepInitializerListConstruction, NonDescendingThrowsError)
{
	EXPECT_THROW(sun::Irrep alpha({5, 4, 2, 3, 3}), std::invalid_argument);
}


TEST(IrrepTranspose, TransposeValues)
{
    sun::Irrep alpha({5, 4, 4, 2});
    std::vector<unsigned int> alphaT = alpha.transpose();
    
    expect_equal_vec(alphaT, std::vector<unsigned int>({4, 4, 3, 3, 1}));
}


TEST(IrrepTranspose, TransposeInvolution)
{
    sun::Irrep alpha({5, 4, 4, 2});
    std::vector<unsigned int> vecT = alpha.transpose();
    sun::Irrep alphaT(vecT);
    std::vector<unsigned int> vec = alphaT.transpose();
    
    expect_equal_vec(alpha.get_vector(), vec);
}
