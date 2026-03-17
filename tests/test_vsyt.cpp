// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <gtest/gtest.h>
#include <vector>
#include <cstdint>

#include "sun/syt/vsyt.h"
#include "testutils.h"

typedef sun::vSYT<int8_t> SYT;


TEST(vSYTConstruction, DefaultConstructor)
{
    SYT y;
    EXPECT_EQ(y.container().size(), 0);
}


TEST(vSYTConstruction, SizeConstructor)
{
    SYT y(12);
    EXPECT_EQ(y.container().size(), 12);
    EXPECT_EQ(y.n(), 12);
    for (unsigned int i = 0; i < y.n(); ++i) {
		EXPECT_EQ(y[i], 0);
	}
}


TEST(vSYTConstruction, VectorConstructor)
{
	std::vector<int> yv({1, 0, 0});
    SYT y(yv);
    EXPECT_EQ(y.n(), 3);
    for (unsigned int i = 0; i < y.n(); ++i) {
		EXPECT_EQ(y[i], yv[i]);
	}
}


TEST(vSYTConstruction, InititializerListConstructor)
{
    SYT y({0, 1, 0});
    EXPECT_EQ(y.n(), 3);
    EXPECT_EQ(y[0], 0);
    EXPECT_EQ(y[1], 1);
    EXPECT_EQ(y[2], 0);
}


TEST(vSYTSet, SetValue)
{
    SYT y1(3);
    y1.set(0, 0);
    y1.set(1, 0);
    y1.set(2, 1);
    expect_equal_vec(y1.container(), std::vector<int8_t>({0, 0, 1}));
    
    SYT y2(3);
    y2.set(0, 0);
    y2.set(1, 1);
    y2.set(2, 0);
    expect_equal_vec(y2.container(), std::vector<int8_t>({0, 1, 0}));
}


TEST(vSYTGet, GetValue)
{
    SYT y1(4);
    y1.set(0, 0);
    y1.set(1, 0);
    y1.set(2, 1);
    
    EXPECT_EQ(y1.get(0), 0);
    EXPECT_EQ(y1.get(1), 0);
    EXPECT_EQ(y1.get(2), 1);
    EXPECT_EQ(y1.get(3), 0);
    
    SYT y2(4);
    y2.set(0, 0);
    y2.set(1, 1);
    y2.set(2, 0);
    
    EXPECT_EQ(y2.get(0), 0);
    EXPECT_EQ(y2.get(1), 1);
    EXPECT_EQ(y2.get(2), 0);
    EXPECT_EQ(y2.get(3), 0);
}


TEST(vSYTExchange, TwoWayExchangeValue)
{
    SYT y1(3);
    y1.set(0, 0);
    y1.set(1, 0);
    y1.set(2, 1);
    expect_equal_vec(y1.container(), std::vector<int8_t>({0, 0, 1}));
    
    SYT y2(y1);
    y2.exchange(1, 2);
    expect_equal_vec(y2.container(), std::vector<int8_t>({0, 1, 0}));
}


TEST(vSYTExchange, DoubleTwoWayExchange)
{
    SYT y1(3);
    y1.set(0, 0);
    y1.set(1, 0);
    y1.set(2, 1);
    expect_equal_vec(y1.container(), std::vector<int8_t>({0, 0, 1}));
    
    SYT y2(y1);
    y2.exchange(1, 2);
    y2.exchange(1, 2);
    
    expect_equal_vec(y2.container(), std::vector<int8_t>({0, 0, 1}));
}


TEST(vSYTExchange, ThreeWayExchangeValue)
{
	// 0 1 3
	// 2 4
	// 5
    SYT y1(6);
    y1.set(0, 0);
    y1.set(1, 0);
    y1.set(2, 1);
    y1.set(3, 0);
    y1.set(4, 1);
    y1.set(5, 2);
    expect_equal_vec(y1.container(), std::vector<int8_t>({0, 0, 1, 0, 1, 2}));
    
    SYT y2(y1);
    y2.exchange(3, 4, 5);
    
    EXPECT_EQ(y2.get(3), 2);
    EXPECT_EQ(y2.get(4), 0);
    EXPECT_EQ(y2.get(5), 1);
    
    SYT y3(y2);
    y3.exchange(3, 4, 5);
    
    EXPECT_EQ(y3.get(3), 1);
    EXPECT_EQ(y3.get(4), 2);
    EXPECT_EQ(y3.get(5), 0);
    
    SYT y4(y3);
    y4.exchange(3, 4, 5);
    EXPECT_EQ(y4.container(), std::vector<int8_t>({0, 0, 1, 0, 1, 2})); // back to original SYT
}


TEST(vSYTOperatorComparison, OperatorEqual)
{
	// 0 1 3
	// 2 4
	// 5
	SYT y1(6);
    y1.set(0, 0);
    y1.set(1, 0);
    y1.set(2, 1);
    y1.set(3, 0);
    y1.set(4, 1);
    y1.set(5, 2);
    
    SYT y2(y1);
    EXPECT_EQ(y1==y2, true);
    y2.exchange(1, 2);
    EXPECT_EQ(y1==y2, false);
    y2.exchange(3, 4, 5);
    EXPECT_EQ(y1==y2, false);
    y2.exchange(3, 4, 5);
    EXPECT_EQ(y1==y2, false);
    y2.exchange(1, 2);
    EXPECT_EQ(y1==y2, false);
    y2.exchange(3, 4, 5);
    EXPECT_EQ(y1==y2, true);
}


TEST(vSYTOperatorComparison, OperatorLessThan)
{
	SYT y1(3);
    y1.set(0, 0);
    y1.set(1, 0);
    y1.set(2, 1);
    expect_equal_vec(y1.container(), std::vector<int8_t>({0, 0, 1}));
    
    SYT y2(3);
    y2.set(0, 0);
    y2.set(1, 1);
    y2.set(2, 0);
    expect_equal_vec(y2.container(), std::vector<int8_t>({0, 1, 0}));
    
    EXPECT_EQ(y1<y2, false);
    EXPECT_EQ(y2<y1, true);
    EXPECT_EQ(y1<y1, false);
    EXPECT_EQ(y2<y2, false);
}


TEST(vSYTOperatorComparison, OperatorLessOrEqualThan)
{
	SYT y1(3);
    y1.set(0, 0);
    y1.set(1, 0);
    y1.set(2, 1);
    expect_equal_vec(y1.container(), std::vector<int8_t>({0, 0, 1}));
    
    SYT y2(3);
    y2.set(0, 0);
    y2.set(1, 1);
    y2.set(2, 0);
    expect_equal_vec(y2.container(), std::vector<int8_t>({0, 1, 0}));
    
    EXPECT_EQ(y1<=y2, false);
    EXPECT_EQ(y2<=y1, true);
    EXPECT_EQ(y1<=y1, true);
    EXPECT_EQ(y2<=y2, true);
}


TEST(vSYTOperatorComparison, OperatorGreaterThan)
{
	SYT y1(3);
    y1.set(0, 0);
    y1.set(1, 0);
    y1.set(2, 1);
    expect_equal_vec(y1.container(), std::vector<int8_t>({0, 0, 1}));
    
    SYT y2(3);
    y2.set(0, 0);
    y2.set(1, 1);
    y2.set(2, 0);
    expect_equal_vec(y2.container(), std::vector<int8_t>({0, 1, 0}));
    
    EXPECT_EQ(y1>y2, true);
    EXPECT_EQ(y2>y1, false);
    EXPECT_EQ(y1>y1, false);
    EXPECT_EQ(y2>y2, false);
}


TEST(vSYTOperatorComparison, OperatorGreaterOrEqualThan)
{
	SYT y1(3);
    y1.set(0, 0);
    y1.set(1, 0);
    y1.set(2, 1);
    expect_equal_vec(y1.container(), std::vector<int8_t>({0, 0, 1}));
    
    SYT y2(3);
    y2.set(0, 0);
    y2.set(1, 1);
    y2.set(2, 0);
    expect_equal_vec(y2.container(), std::vector<int8_t>({0, 1, 0}));
    
    EXPECT_EQ(y1>=y2, true);
    EXPECT_EQ(y2>=y1, false);
    EXPECT_EQ(y1>=y1, true);
    EXPECT_EQ(y2>=y2, true);
}
