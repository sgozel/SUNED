// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <gtest/gtest.h>
#include <vector>
#include <utility>
#include <cstdint>

#include "testutils.h"
#include "sun/irrep/irrep.h"
#include "sun/syt/vsyt.h"
#include "sun/syt_usage/vsyt_usage.h"


typedef sun::vSYT<int8_t> SYT;


TEST(vSYTUsage, GetSYT)
{
	sun::Irrep alpha({3, 2});
	std::vector<SYT> Y = sun::get_SYT<int8_t>(alpha);
	
    expect_equal_vec(Y[0].container(), std::vector<int8_t>({0, 1, 0, 1, 0}));
    expect_equal_vec(Y[1].container(), std::vector<int8_t>({0, 0, 1, 1, 0}));
    expect_equal_vec(Y[2].container(), std::vector<int8_t>({0, 1, 0, 0, 1}));
    expect_equal_vec(Y[3].container(), std::vector<int8_t>({0, 0, 1, 0, 1}));
    expect_equal_vec(Y[4].container(), std::vector<int8_t>({0, 0, 0, 1, 1}));
}


TEST(vSYTUsage, OperatorLessThan)
{
    sun::Irrep alpha({3, 2});
	std::vector<SYT> Y = sun::get_SYT<int8_t>(alpha);
	
    EXPECT_EQ(Y[0]<Y[1], true);
    EXPECT_EQ(Y[1]<Y[2], true);
    EXPECT_EQ(Y[2]<Y[3], true);
    EXPECT_EQ(Y[3]<Y[4], true);
    
    EXPECT_EQ(Y[3]<Y[3], false);
    EXPECT_EQ(Y[4]<Y[3], false);
}


TEST(vSYTUsage, OperatorGreaterThan)
{
    sun::Irrep alpha({3, 2});
	std::vector<SYT> Y = sun::get_SYT<int8_t>(alpha);
	
    EXPECT_EQ(Y[1]>Y[0], true);
    EXPECT_EQ(Y[2]>Y[1], true);
    EXPECT_EQ(Y[3]>Y[2], true);
    EXPECT_EQ(Y[4]>Y[3], true);
    
    EXPECT_EQ(Y[3]>Y[0], true);
    EXPECT_EQ(Y[2]>Y[4], false);
    EXPECT_EQ(Y[1]>Y[1], false);
}


TEST(vSYTUsage, OperatorLessOrEqualThan)
{
    sun::Irrep alpha({3, 2});
	std::vector<SYT> Y = sun::get_SYT<int8_t>(alpha);
	
    EXPECT_EQ(Y[0]<=Y[1], true);
    EXPECT_EQ(Y[1]<=Y[2], true);
    EXPECT_EQ(Y[2]<=Y[3], true);
    EXPECT_EQ(Y[3]<=Y[4], true);
    
    EXPECT_EQ(Y[3]<=Y[0], false);
    EXPECT_EQ(Y[2]<=Y[4], true);
    EXPECT_EQ(Y[1]<=Y[1], true);
}


TEST(vSYTUsage, OperatorGreaterOrEqualThan)
{
    sun::Irrep alpha({3, 2});
	std::vector<SYT> Y = sun::get_SYT<int8_t>(alpha);
	
    EXPECT_EQ(Y[1]>=Y[0], true);
    EXPECT_EQ(Y[2]>=Y[1], true);
    EXPECT_EQ(Y[3]>=Y[2], true);
    EXPECT_EQ(Y[4]>=Y[3], true);
    
    EXPECT_EQ(Y[3]>=Y[0], true);
    EXPECT_EQ(Y[2]>=Y[4], false);
    EXPECT_EQ(Y[1]>=Y[1], true);
}


TEST(vSYTUsage, OperatorEqual)
{
    sun::Irrep alpha({3, 2});
	std::vector<SYT> Y = sun::get_SYT<int8_t>(alpha);
	
    EXPECT_EQ(Y[0]==Y[0], true);
    EXPECT_EQ(Y[0]==Y[1], false);
    EXPECT_EQ(Y[1]==Y[1], true);
    EXPECT_EQ(Y[3]==Y[4], false);
}


TEST(vSYTUsage, GetSYTParts)
{
	sun::Irrep alpha({3, 2});
	const UINT64 falpha = 5;
	std::vector<SYT> expected(falpha);
	expected[0] = SYT({0, 1, 0, 1, 0});
    expected[1] = SYT({0, 0, 1, 1, 0});
    expected[2] = SYT({0, 1, 0, 0, 1});
    expected[3] = SYT({0, 0, 1, 0, 1});
    expected[4] = SYT({0, 0, 0, 1, 1});
	
	std::vector<SYT> Y;
	
	UINT64 dim = 1;
	for (UINT64 from = 0; from < falpha; ++from) {
		Y = sun::get_SYT<int8_t>(alpha, from, dim);
		EXPECT_EQ(Y.size(), static_cast<size_t>(dim));
		EXPECT_EQ(Y[0], expected[from]);
	}
	
	dim = 2;
	for (UINT64 from = 0; from < falpha-1; ++from) {
		
		UINT64 dimeff = (from + dim > falpha ? dim-1 : dim);
		
		Y = sun::get_SYT<int8_t>(alpha, from, dimeff);
		
		EXPECT_EQ(Y.size(), static_cast<size_t>(dimeff));
		EXPECT_EQ(Y[0], expected[from]);
		if (dimeff == 2) {
			EXPECT_EQ(Y[1], expected[from+1]);
		}
	}
	
	Y = sun::get_SYT<int8_t>(alpha, 1, 3);
	EXPECT_EQ(Y.size(), 3);
	EXPECT_EQ(Y[0], expected[1]);
	EXPECT_EQ(Y[1], expected[2]);
	EXPECT_EQ(Y[2], expected[3]);
	
	Y = sun::get_SYT<int8_t>(alpha, 2, 3);
	EXPECT_EQ(Y.size(), 3);
	EXPECT_EQ(Y[0], expected[2]);
	EXPECT_EQ(Y[1], expected[3]);
	EXPECT_EQ(Y[2], expected[4]);
}


TEST(vSYTUsage, GetSYTPartsDimThrows)
{
	sun::Irrep alpha({3, 2});
	
	EXPECT_THROW(sun::get_SYT<int8_t>(alpha, 0, 6), std::runtime_error);
	EXPECT_THROW(sun::get_SYT<int8_t>(alpha, 1, 5), std::runtime_error);
	EXPECT_THROW(sun::get_SYT<int8_t>(alpha, 2, 4), std::runtime_error);
	EXPECT_THROW(sun::get_SYT<int8_t>(alpha, 3, 3), std::runtime_error);
	EXPECT_THROW(sun::get_SYT<int8_t>(alpha, 4, 2), std::runtime_error);
}


TEST(vSYTUsage, GetSYTPartsFromThrows)
{
	sun::Irrep alpha({3, 2});
	
	EXPECT_THROW(sun::get_SYT<int8_t>(alpha, 5, 1), std::runtime_error);
}


TEST(vSYTUsage, GetColumn)
{
	sun::Irrep alpha({3, 2});
	unsigned int Ns = alpha.n();
	
	std::vector<std::vector<int8_t>> expected(5, std::vector<int8_t>(Ns));
	expected[0] = {0, 0, 1, 1, 2};
	expected[1] = {0, 1, 0, 1, 2};
	expected[2] = {0, 0, 1, 2, 1};
	expected[3] = {0, 1, 0, 2, 1};
	expected[4] = {0, 1, 2, 0, 1};
	
	std::vector<SYT> Y = sun::get_SYT<int8_t>(alpha);
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


TEST(vSYTUsage, GetColumnkkplusone)
{
	sun::Irrep alpha({3, 2});
	unsigned int Ns = alpha.n();
	
	std::vector<std::vector<int>> expected(5, std::vector<int>(Ns));
	expected[0] = {0, 0, 1, 1, 2};
	expected[1] = {0, 1, 0, 1, 2};
	expected[2] = {0, 0, 1, 2, 1};
	expected[3] = {0, 1, 0, 2, 1};
	expected[4] = {0, 1, 2, 0, 1};
	
	std::vector<SYT> Y = sun::get_SYT<int8_t>(alpha);
	
	for (unsigned int i = 0; i < Y.size(); ++i) {
		for (unsigned int k = 0; k < Ns-1; ++k) {
			std::pair<int, int> cols = sun::get_column_k_k_plus_one(Y[i], k);
			EXPECT_EQ(cols.first, expected[i][k]);
			EXPECT_EQ(cols.second, expected[i][k+1]);
		}
	}
}


TEST(vSYTUsage, GetAxialDistance)
{
	sun::Irrep alpha({3, 2});
	std::vector<SYT> Y = sun::get_SYT<int8_t>(alpha);
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


TEST(vSYTUsage, FullSimplify)
{
	sun::Irrep alpha({3, 2});
	std::vector<SYT> Y = sun::get_SYT<int8_t>(alpha);
	std::vector<SYT> ydev = Y;
	ydev.insert(ydev.end(), Y.begin(), Y.end());
	ydev.pop_back();
	
	std::vector<double> coeffs({1.0, 0.0, -1.0, 5.0, 3.14,
								2.0, 1.0,  1.0, 1.0});
	
	sun::fullsimplify_development(ydev, coeffs);
	
	EXPECT_EQ(ydev.size(), 4);
	EXPECT_EQ(ydev[0], Y[0]);
	EXPECT_EQ(ydev[1], Y[1]);
	EXPECT_EQ(ydev[2], Y[3]);
	EXPECT_EQ(ydev[3], Y[4]);
	EXPECT_NEAR(coeffs[0], 3.0, 1.0e-14);
	EXPECT_NEAR(coeffs[1], 1.0, 1.0e-14);
	EXPECT_NEAR(coeffs[2], 6.0, 1.0e-14);
	EXPECT_NEAR(coeffs[3], 3.14, 1.0e-14);
}


TEST(vSYTUsage, DevelopConsecutiveNumber)
{
	sun::Irrep alpha({3, 2});
	std::vector<SYT> Y = sun::get_SYT<int8_t>(alpha);
	std::vector<SYT> ydev_in(1);
	ydev_in[0] = Y[0];
	std::vector<double> coeffs_in({1.0});
	std::vector<SYT> ydev_out;
	std::vector<double> coeffs_out;
	
	const unsigned int k = 1;
	sun::develop_consecutive_number(alpha, ydev_in, coeffs_in, ydev_out, coeffs_out, k);
	
	EXPECT_EQ(ydev_out.size(), 2);
	EXPECT_EQ(ydev_out[0], Y[0]);
	EXPECT_EQ(ydev_out[1], Y[1]);
	EXPECT_NEAR(coeffs_out[0], 0.5, 1.0e-14);
	EXPECT_NEAR(coeffs_out[1], 0.8660254037844386, 1.0e-14);
}


TEST(vSYTUsage, DevelopConsecutiveNumberInPlace)
{
	sun::Irrep alpha({3, 2});
	std::vector<SYT> Y = sun::get_SYT<int8_t>(alpha);	
	std::vector<SYT> ydev(1);
	ydev[0] = Y[0];
	
	std::vector<double> coeffs({1.0});
	
	const unsigned int k = 1;
	sun::develop_consecutive_number_inplace(alpha, ydev, coeffs, k);
	
	EXPECT_EQ(ydev.size(), 2);
	EXPECT_EQ(ydev[0], Y[0]);
	EXPECT_EQ(ydev[1], Y[1]);
	EXPECT_NEAR(coeffs[0], 0.5, 1.0e-14);
	EXPECT_NEAR(coeffs[1], 0.8660254037844386, 1.0e-14);
}
