// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <gtest/gtest.h>

#include "sun/syt/bsyt.h"


TEST(bSYTConstruction, DefaultConstructor)
{
    sun::tbSYT y;
    EXPECT_EQ(y.value(), 0ULL);
}


TEST(bSYTSet, SetValue)
{
    sun::bSYT<uint64_t, 2> y1;
    y1.set(0, 0);
    y1.set(1, 0);
    y1.set(2, 1);
    EXPECT_EQ(y1.value(), 16ULL);
    
    sun::bSYT<uint64_t, 2> y2;
    y2.set(0, 0);
    y2.set(1, 1);
    y2.set(2, 0);
    EXPECT_EQ(y2.value(), 4ULL);
}


TEST(bSYTGet, GetValue)
{
    sun::bSYT<uint64_t, 2> y1;
    y1.set(0, 0);
    y1.set(1, 0);
    y1.set(2, 1);
    
    EXPECT_EQ(y1.get(0), 0);
    EXPECT_EQ(y1.get(1), 0);
    EXPECT_EQ(y1.get(2), 1);
    EXPECT_EQ(y1.get(3), 0);
    
    sun::bSYT<uint64_t, 2> y2;
    y2.set(0, 0);
    y2.set(1, 1);
    y2.set(2, 0);
    
    EXPECT_EQ(y2.get(0), 0);
    EXPECT_EQ(y2.get(1), 1);
    EXPECT_EQ(y2.get(2), 0);
    EXPECT_EQ(y2.get(3), 0);
}


TEST(bSYTExchange, TwoWayExchangeValue)
{
    sun::bSYT<uint64_t, 2> y1;
    y1.set(0, 0);
    y1.set(1, 0);
    y1.set(2, 1);
    
    EXPECT_EQ(y1.value(), 16ULL);
    
    sun::bSYT<uint64_t, 2> y2(y1);
    y2.exchange(1, 2);
    
    EXPECT_EQ(y2.value(), 4ULL);
}


TEST(bSYTExchange, DoubleTwoWayExchange)
{
    sun::bSYT<uint64_t, 2> y1;
    y1.set(0, 0);
    y1.set(1, 0);
    y1.set(2, 1);
    
    EXPECT_EQ(y1.value(), 16ULL);
    
    sun::bSYT<uint64_t, 2> y2(y1);
    y2.exchange(1, 2);
    y2.exchange(1, 2);
    
    EXPECT_EQ(y2.value(), 16ULL);
}


TEST(bSYTExchange, ThreeWayExchangeValue)
{
	// 0 1 3
	// 2 4
	// 5
    sun::bSYT<uint64_t, 2> y1;
    y1.set(0, 0);
    y1.set(1, 0);
    y1.set(2, 1);
    y1.set(3, 0);
    y1.set(4, 1);
    y1.set(5, 2);
    
    EXPECT_EQ(y1.value(), 2320ULL);
    
    sun::bSYT<uint64_t, 2> y2(y1);
    y2.exchange(3, 4, 5);
    
    EXPECT_EQ(y2.value(), 1168ULL);
    EXPECT_EQ(y2.get(3), 2);
    EXPECT_EQ(y2.get(4), 0);
    EXPECT_EQ(y2.get(5), 1);
    
    sun::bSYT<uint64_t, 2> y3(y2);
    y3.exchange(3, 4, 5);
    
    EXPECT_EQ(y3.value(), 592ULL);
    EXPECT_EQ(y3.get(3), 1);
    EXPECT_EQ(y3.get(4), 2);
    EXPECT_EQ(y3.get(5), 0);
    
    sun::bSYT<uint64_t, 2> y4(y3);
    y4.exchange(3, 4, 5);
    EXPECT_EQ(y4.value(), 2320ULL); // back to original SYT
}


TEST(bSYTOperatorComparison, OperatorEqual)
{
	// 0 1 3
	// 2 4
	// 5
	sun::bSYT<uint64_t, 2> y1;
    y1.set(0, 0);
    y1.set(1, 0);
    y1.set(2, 1);
    y1.set(3, 0);
    y1.set(4, 1);
    y1.set(5, 2);
    
    sun::bSYT<uint64_t, 2> y2(y1);
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


TEST(bSYTOperatorComparison, OperatorLessThan)
{
	sun::bSYT<uint64_t, 2> y1;
    y1.set(0, 0);
    y1.set(1, 0);
    y1.set(2, 1);
    EXPECT_EQ(y1.value(), 16ULL);
    
    sun::bSYT<uint64_t, 2> y2;
    y2.set(0, 0);
    y2.set(1, 1);
    y2.set(2, 0);
    EXPECT_EQ(y2.value(), 4ULL);
    
    EXPECT_EQ(y1<y2, false);
    EXPECT_EQ(y2<y1, true);
    EXPECT_EQ(y1<y1, false);
    EXPECT_EQ(y2<y2, false);
}


TEST(bSYTOperatorComparison, OperatorLessOrEqualThan)
{
	sun::bSYT<uint64_t, 2> y1;
    y1.set(0, 0);
    y1.set(1, 0);
    y1.set(2, 1);
    EXPECT_EQ(y1.value(), 16ULL);
    
    sun::bSYT<uint64_t, 2> y2;
    y2.set(0, 0);
    y2.set(1, 1);
    y2.set(2, 0);
    EXPECT_EQ(y2.value(), 4ULL);
    
    EXPECT_EQ(y1<=y2, false);
    EXPECT_EQ(y2<=y1, true);
    EXPECT_EQ(y1<=y1, true);
    EXPECT_EQ(y2<=y2, true);
}


TEST(bSYTOperatorComparison, OperatorGreaterThan)
{
	sun::bSYT<uint64_t, 2> y1;
    y1.set(0, 0);
    y1.set(1, 0);
    y1.set(2, 1);
    EXPECT_EQ(y1.value(), 16ULL);
    
    sun::bSYT<uint64_t, 2> y2;
    y2.set(0, 0);
    y2.set(1, 1);
    y2.set(2, 0);
    EXPECT_EQ(y2.value(), 4ULL);
    
    EXPECT_EQ(y1>y2, true);
    EXPECT_EQ(y2>y1, false);
    EXPECT_EQ(y1>y1, false);
    EXPECT_EQ(y2>y2, false);
}


TEST(bSYTOperatorComparison, OperatorGreaterOrEqualThan)
{
	sun::bSYT<uint64_t, 2> y1;
    y1.set(0, 0);
    y1.set(1, 0);
    y1.set(2, 1);
    EXPECT_EQ(y1.value(), 16ULL);
    
    sun::bSYT<uint64_t, 2> y2;
    y2.set(0, 0);
    y2.set(1, 1);
    y2.set(2, 0);
    EXPECT_EQ(y2.value(), 4ULL);
    
    EXPECT_EQ(y1>=y2, true);
    EXPECT_EQ(y2>=y1, false);
    EXPECT_EQ(y1>=y1, true);
    EXPECT_EQ(y2>=y2, true);
}


TEST(bSYTCheck, ContainerBitExtent)
{
	sun::bSYT<uint64_t, 2> y;
	constexpr unsigned int extent = y.bitextent();
	EXPECT_EQ(extent, 2);
	
	constexpr unsigned int maxN = 1<<extent;
	unsigned int n = 4;
	
	EXPECT_NO_THROW(sun::tbSYT::check(maxN-1, n));
	EXPECT_NO_THROW(sun::tbSYT::check(maxN, n));
	EXPECT_THROW(sun::tbSYT::check(maxN+1, n), std::runtime_error);
	EXPECT_THROW(sun::tbSYT::check(maxN+2, n), std::runtime_error);
}
