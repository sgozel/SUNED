// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <gtest/gtest.h>
#include <vector>
#include <utility>

#include "sun/irrep/irrep.h"
#include "sun/syt/bsyt.h"
#include "sun/syt_usage/bsyt_usage.h"


TEST(SYTUsage, GetSYT)
{
	sun::Irrep alpha({3, 2});
	std::vector<sun::tbSYT> Y = sun::get_SYT(alpha);
	
    EXPECT_EQ(Y[0].value(), 68ULL);
    EXPECT_EQ(Y[1].value(), 80ULL);
    EXPECT_EQ(Y[2].value(), 260ULL);
    EXPECT_EQ(Y[3].value(), 272ULL);
    EXPECT_EQ(Y[4].value(), 320ULL);
}


TEST(SYTUsage, OperatorLessThan)
{
    sun::Irrep alpha({3, 2});
	std::vector<sun::tbSYT> Y = sun::get_SYT(alpha);
	
    EXPECT_EQ(Y[0]<Y[1], true);
    EXPECT_EQ(Y[1]<Y[2], true);
    EXPECT_EQ(Y[2]<Y[3], true);
    EXPECT_EQ(Y[3]<Y[4], true);
    
    EXPECT_EQ(Y[3]<Y[3], false);
    EXPECT_EQ(Y[4]<Y[3], false);
}


TEST(SYTUsage, OperatorGreaterThan)
{
    sun::Irrep alpha({3, 2});
	std::vector<sun::tbSYT> Y = sun::get_SYT(alpha);
	
    EXPECT_EQ(Y[1]>Y[0], true);
    EXPECT_EQ(Y[2]>Y[1], true);
    EXPECT_EQ(Y[3]>Y[2], true);
    EXPECT_EQ(Y[4]>Y[3], true);
    
    EXPECT_EQ(Y[3]>Y[0], true);
    EXPECT_EQ(Y[2]>Y[4], false);
    EXPECT_EQ(Y[1]>Y[1], false);
}


TEST(SYTUsage, OperatorLessOrEqualThan)
{
    sun::Irrep alpha({3, 2});
	std::vector<sun::tbSYT> Y = sun::get_SYT(alpha);
	
    EXPECT_EQ(Y[0]<=Y[1], true);
    EXPECT_EQ(Y[1]<=Y[2], true);
    EXPECT_EQ(Y[2]<=Y[3], true);
    EXPECT_EQ(Y[3]<=Y[4], true);
    
    EXPECT_EQ(Y[3]<=Y[0], false);
    EXPECT_EQ(Y[2]<=Y[4], true);
    EXPECT_EQ(Y[1]<=Y[1], true);
}


TEST(SYTUsage, OperatorGreaterOrEqualThan)
{
    sun::Irrep alpha({3, 2});
	std::vector<sun::tbSYT> Y = sun::get_SYT(alpha);
	
    EXPECT_EQ(Y[1]>=Y[0], true);
    EXPECT_EQ(Y[2]>=Y[1], true);
    EXPECT_EQ(Y[3]>=Y[2], true);
    EXPECT_EQ(Y[4]>=Y[3], true);
    
    EXPECT_EQ(Y[3]>=Y[0], true);
    EXPECT_EQ(Y[2]>=Y[4], false);
    EXPECT_EQ(Y[1]>=Y[1], true);
}


TEST(SYTUsage, OperatorEqual)
{
    sun::Irrep alpha({3, 2});
	std::vector<sun::tbSYT> Y = sun::get_SYT(alpha);
	
    EXPECT_EQ(Y[0]==Y[0], true);
    EXPECT_EQ(Y[0]==Y[1], false);
    EXPECT_EQ(Y[1]==Y[1], true);
    EXPECT_EQ(Y[3]==Y[4], false);
}


TEST(SYTUsage, GetSYTParts)
{
	sun::Irrep alpha({3, 2});
	const std::vector<UINT64> expected({68ULL, 80ULL, 260ULL, 272ULL, 320ULL});
	const UINT64 falpha = 5;
	
	std::vector<sun::tbSYT> Y;
	
	UINT64 dim = 1;
	for (UINT64 from = 0; from < falpha; ++from) {
		Y = sun::get_SYT(alpha, from, dim);
		EXPECT_EQ(Y.size(), static_cast<size_t>(dim));
		EXPECT_EQ(Y[0].value(), expected[from]);
	}
	
	dim = 2;
	for (UINT64 from = 0; from < falpha-1; ++from) {
		
		UINT64 dimeff = (from + dim > falpha ? dim-1 : dim);
		
		Y = sun::get_SYT(alpha, from, dimeff);
		
		EXPECT_EQ(Y.size(), static_cast<size_t>(dimeff));
		EXPECT_EQ(Y[0].value(), expected[from]);
		if (dimeff == 2) {
			EXPECT_EQ(Y[1].value(), expected[from+1]);
		}
	}
	
	Y = sun::get_SYT(alpha, 1, 3);
	EXPECT_EQ(Y.size(), 3);
	EXPECT_EQ(Y[0].value(), expected[1]);
	EXPECT_EQ(Y[1].value(), expected[2]);
	EXPECT_EQ(Y[2].value(), expected[3]);
	
	Y = sun::get_SYT(alpha, 2, 3);
	EXPECT_EQ(Y.size(), 3);
	EXPECT_EQ(Y[0].value(), expected[2]);
	EXPECT_EQ(Y[1].value(), expected[3]);
	EXPECT_EQ(Y[2].value(), expected[4]);
}


TEST(SYTUsage, GetSYTPartsDimThrows)
{
	sun::Irrep alpha({3, 2});
	
	EXPECT_THROW(sun::get_SYT(alpha, 0, 6), std::runtime_error);
	EXPECT_THROW(sun::get_SYT(alpha, 1, 5), std::runtime_error);
	EXPECT_THROW(sun::get_SYT(alpha, 2, 4), std::runtime_error);
	EXPECT_THROW(sun::get_SYT(alpha, 3, 3), std::runtime_error);
	EXPECT_THROW(sun::get_SYT(alpha, 4, 2), std::runtime_error);
}


TEST(SYTUsage, GetSYTPartsFromThrows)
{
	sun::Irrep alpha({3, 2});
	
	EXPECT_THROW(sun::get_SYT(alpha, 5, 1), std::runtime_error);
}


TEST(SYTUsage, GetColumn)
{
	sun::Irrep alpha({3, 2});
	unsigned int Ns = alpha.n();
	
	std::vector<std::vector<int>> expected(5, std::vector<int>(Ns));
	expected[0] = {0, 0, 1, 1, 2};
	expected[1] = {0, 1, 0, 1, 2};
	expected[2] = {0, 0, 1, 2, 1};
	expected[3] = {0, 1, 0, 2, 1};
	expected[4] = {0, 1, 2, 0, 1};
	
	std::vector<sun::tbSYT> Y = sun::get_SYT(alpha);
	std::vector<std::vector<int>> CY(Y.size(), std::vector<int>(Ns));
	for (unsigned int i = 0; i < Y.size(); ++i) {
		CY[i] = sun::get_column(Y[i], alpha);
	}
	
	for (unsigned int i = 0; i < Y.size(); ++i) {
		for (unsigned int k = 0; k < Ns; ++k) {
			EXPECT_EQ(CY[i][k], expected[i][k]);
		}
	}
}


TEST(SYTUsage, GetColumnkkplusone)
{
	sun::Irrep alpha({3, 2});
	unsigned int Ns = alpha.n();
	
	std::vector<std::vector<int>> expected(5, std::vector<int>(Ns));
	expected[0] = {0, 0, 1, 1, 2};
	expected[1] = {0, 1, 0, 1, 2};
	expected[2] = {0, 0, 1, 2, 1};
	expected[3] = {0, 1, 0, 2, 1};
	expected[4] = {0, 1, 2, 0, 1};
	
	std::vector<sun::tbSYT> Y = sun::get_SYT(alpha);
	
	for (unsigned int i = 0; i < Y.size(); ++i) {
		for (unsigned int k = 0; k < Ns-1; ++k) {
			std::pair<int, int> cols = sun::get_column_k_k_plus_one(Y[i], k);
			EXPECT_EQ(cols.first, expected[i][k]);
			EXPECT_EQ(cols.second, expected[i][k+1]);
		}
	}
}


TEST(SYTUsage, GetAxialDistance)
{
	sun::Irrep alpha({3, 2});
	std::vector<sun::tbSYT> Y = sun::get_SYT(alpha);
	std::vector<int> cy = sun::get_column(Y[0], alpha);
	
	EXPECT_EQ(sun::get_axial_distance(Y[0], cy, 0, 1), 1);
	EXPECT_EQ(sun::get_axial_distance(Y[0], cy, 0, 2), -1);
	EXPECT_EQ(sun::get_axial_distance(Y[0], cy, 0, 3), 0);
	EXPECT_EQ(sun::get_axial_distance(Y[0], cy, 0, 4), -2);
	EXPECT_EQ(sun::get_axial_distance(Y[0], cy, 1, 2), -2);
	EXPECT_EQ(sun::get_axial_distance(Y[0], cy, 1, 3), -1);
	EXPECT_EQ(sun::get_axial_distance(Y[0], cy, 1, 4), -3);
	EXPECT_EQ(sun::get_axial_distance(Y[0], cy, 2, 3), 1);
	EXPECT_EQ(sun::get_axial_distance(Y[0], cy, 2, 4), -1);
	EXPECT_EQ(sun::get_axial_distance(Y[0], cy, 3, 4), -2);
}


TEST(SYTUsage, FullSimplify)
{
	sun::Irrep alpha({3, 2});
	std::vector<sun::tbSYT> Y = sun::get_SYT(alpha);
	std::vector<sun::tbSYT> ydev = Y;
	ydev.insert(ydev.end(), Y.begin(), Y.end());
	ydev.pop_back();
	
	std::vector<double> coeffs({1.0, 0.0, -1.0, 5.0, 3.14,
								2.0, 1.0,  1.0, 1.0});
	
	sun::fullsimplify_development(ydev, coeffs);
	
	EXPECT_EQ(ydev.size(), 4);
	EXPECT_EQ(ydev[0].value(), 68ULL);
	EXPECT_EQ(ydev[1].value(), 80ULL);
	EXPECT_EQ(ydev[2].value(), 272ULL);
	EXPECT_EQ(ydev[3].value(), 320ULL);
	EXPECT_NEAR(coeffs[0], 3.0, 1.0e-14);
	EXPECT_NEAR(coeffs[1], 1.0, 1.0e-14);
	EXPECT_NEAR(coeffs[2], 6.0, 1.0e-14);
	EXPECT_NEAR(coeffs[3], 3.14, 1.0e-14);
}


TEST(SYTUsage, DevelopConsecutiveNumber)
{
	sun::Irrep alpha({3, 2});
	std::vector<sun::tbSYT> Y = sun::get_SYT(alpha);
	std::vector<sun::tbSYT> ydev_in(1);
	ydev_in[0] = Y[0];
	std::vector<double> coeffs_in({1.0});
	std::vector<sun::tbSYT> ydev_out;
	std::vector<double> coeffs_out;
	
	const unsigned int k = 1;
	sun::develop_consecutive_number(alpha, ydev_in, coeffs_in, ydev_out, coeffs_out, k);
	
	EXPECT_EQ(ydev_out.size(), 2);
	EXPECT_EQ(ydev_out[0].value(), 68ULL);
	EXPECT_EQ(ydev_out[1].value(), 80ULL);
	EXPECT_NEAR(coeffs_out[0], 0.5, 1.0e-14);
	EXPECT_NEAR(coeffs_out[1], 0.8660254037844386, 1.0e-14);
}


TEST(SYTUsage, DevelopConsecutiveNumberInPlace)
{
	sun::Irrep alpha({3, 2});
	std::vector<sun::tbSYT> Y = sun::get_SYT(alpha);
	std::vector<sun::tbSYT> ydev(1);
	ydev[0] = Y[0];
	
	std::vector<double> coeffs({1.0});
	
	const unsigned int k = 1;
	sun::develop_consecutive_number_inplace(alpha, ydev, coeffs, k);
	
	EXPECT_EQ(ydev.size(), 2);
	EXPECT_EQ(ydev[0].value(), 68ULL);
	EXPECT_EQ(ydev[1].value(), 80ULL);
	EXPECT_NEAR(coeffs[0], 0.5, 1.0e-14);
	EXPECT_NEAR(coeffs[1], 0.8660254037844386, 1.0e-14);
}
