// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <gtest/gtest.h>
#include <climits>
#include <stdexcept>

#include "common/datatypes.h"

TEST(Datatypes, size_bytes)
{
	EXPECT_EQ(sizeof(UCHAR), 1);
	EXPECT_EQ(sizeof(INT8), 1);
    EXPECT_EQ(sizeof(UINT8), 1);
	EXPECT_EQ(sizeof(INT16), 2);
    EXPECT_EQ(sizeof(UINT16), 2);
    EXPECT_EQ(sizeof(INT32), 4);
    EXPECT_EQ(sizeof(UINT32), 4);
    EXPECT_EQ(sizeof(INT64), 8);
    EXPECT_EQ(sizeof(UINT64), 8);
}

TEST(Datatypes, size_bits)
{
	static_assert(sizeof(INT8) * CHAR_BIT == 8);
    EXPECT_EQ(sizeof(INT8) * CHAR_BIT, 8u);
	static_assert(sizeof(INT16) * CHAR_BIT == 16);
    EXPECT_EQ(sizeof(INT16) * CHAR_BIT, 16u);
	static_assert(sizeof(INT32) * CHAR_BIT == 32);
    EXPECT_EQ(sizeof(INT32) * CHAR_BIT, 32u);
	static_assert(sizeof(UINT32) * CHAR_BIT == 32);
    EXPECT_EQ(sizeof(UINT32) * CHAR_BIT, 32u);
    static_assert(sizeof(INT64) * CHAR_BIT == 64);
    EXPECT_EQ(sizeof(INT64) * CHAR_BIT, 64u);
    static_assert(sizeof(UINT64) * CHAR_BIT == 64);
    EXPECT_EQ(sizeof(UINT64) * CHAR_BIT, 64u);
}
