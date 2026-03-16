// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include <gtest/gtest.h>


template<typename V>
static void expect_equal_vec(const V& a,
							 const V& b)
{
    ASSERT_EQ(a.size(), b.size());
    for (size_t i = 0; i < a.size(); ++i) {
        EXPECT_EQ(a[i], b[i]) << "  at index " << i;
	}
}


template<typename V>
static void expect_near_vec(const V& a,
                            const V& b,
                            double tol = 1e-12)
{
    ASSERT_EQ(a.size(), b.size());
    for (size_t i = 0; i < a.size(); ++i) {
        EXPECT_NEAR(a[i], b[i], tol) << "  at index " << i;
	}
}

#endif

