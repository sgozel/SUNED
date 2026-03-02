// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <stdexcept>

#include "tmatrix/tmatrix.h"
#include "testutils.h"


using namespace lanczos;

TEST(TmatrixConstruction, DefaultConstructor)
{
    Tmatrix t;
    EXPECT_EQ(t.size(), 0u);
}


TEST(TmatrixConstruction, VectorConstructor)
{
    std::vector<double> alpha = {1.0, 2.0, 3.0};
    std::vector<double> beta  = {0.5, 0.5, 0.0};
    Tmatrix t(alpha, beta);
    EXPECT_EQ(t.size(), 3u);
    expect_near_vec(t.get_alpha(), alpha);
    expect_near_vec(t.get_beta(),  beta);
}


TEST(TmatrixConstruction, MismatchedSizesAborts)
{
    std::vector<double> alpha = {1.0, 2.0};
    std::vector<double> beta  = {0.5};
    EXPECT_THROW(Tmatrix(alpha, beta), std::runtime_error);
}


TEST(TmatrixMutation, PushBack)
{
    Tmatrix t;
    t.push_back(1.0, 0.5);
    EXPECT_EQ(t.size(), 1u);
    t.push_back(2.0, 0.3);
    EXPECT_EQ(t.size(), 2u);
    EXPECT_DOUBLE_EQ(t.get_alpha()[0], 1.0);
    EXPECT_DOUBLE_EQ(t.get_alpha()[1], 2.0);
    EXPECT_DOUBLE_EQ(t.get_beta()[0],  0.5);
    EXPECT_DOUBLE_EQ(t.get_beta()[1],  0.3);
}


TEST(TmatrixMutation, PopBack)
{
    Tmatrix t;
    t.push_back(1.0, 0.5);
    t.push_back(2.0, 0.3);
    t.pop_back();
    EXPECT_EQ(t.size(), 1u);
    EXPECT_DOUBLE_EQ(t.get_alpha()[0], 1.0);
    EXPECT_DOUBLE_EQ(t.get_beta()[0], 0.5);
}


TEST(TmatrixMutation, GetLastBeta)
{
    Tmatrix t;
    t.push_back(1.0, 0.7);
    t.push_back(2.0, 0.4);
    EXPECT_DOUBLE_EQ(t.get_last_beta(), 0.4);
}


TEST(TmatrixEigenvalues, OneDimensional)
{
    Tmatrix t;
    t.push_back(3.14, 0.0);
    auto ev = t.eigenvalues();
    ASSERT_EQ(ev.size(), 1u);
    EXPECT_NEAR(ev[0], 3.14, 1e-12);
}


TEST(TmatrixEigenvalues, TwoDimensional)
{
    const double a0 = 2.0, a1 = 4.0, b = 1.0;
    const double mid   = (a0 + a1) / 2.0;
    const double half  = (a0 - a1) / 2.0;
    const double lam1  = mid - std::sqrt(half*half + b*b);
    const double lam2  = mid + std::sqrt(half*half + b*b);

    Tmatrix t;
    t.push_back(a0, b);
    t.push_back(a1, 0.0);

    auto ev = t.eigenvalues();
    ASSERT_EQ(ev.size(), 2u);
    EXPECT_NEAR(ev[0], lam1, 1e-12);
    EXPECT_NEAR(ev[1], lam2, 1e-12);
}


TEST(TmatrixEigenvalues, DiagonalIdentity)
{
    const unsigned int N  = 5;
    const double d = 7.0;
    Tmatrix t;
    for (unsigned int i = 0; i < N; ++i) {
        t.push_back(d, 0.0);
	}
    auto ev = t.eigenvalues();
    ASSERT_EQ(ev.size(), static_cast<size_t>(N));
    for (auto e : ev) {
        EXPECT_NEAR(e, d, 1e-12);
	}
}


TEST(TmatrixEigenvalues, Laplacian1D)
{
    const unsigned int N = 6;
    Tmatrix t;
    for (unsigned int i = 0; i < N; ++i) {
        t.push_back(2.0, -1.0);
	}

    auto ev = t.eigenvalues();
    ASSERT_EQ(ev.size(), static_cast<size_t>(N));

    std::vector<double> expected(N);
    for (unsigned int k = 1; k <= N; ++k) {
        expected[k-1] = 2.0 - 2.0 * std::cos(k * M_PI / (N + 1));
	}
    std::sort(expected.begin(), expected.end());

    for (unsigned int i = 0; i < N; ++i) {
        EXPECT_NEAR(ev[i], expected[i], 1e-12) << "  at index " << i;
	}
}


TEST(TmatrixEigenvalues, CachingConsistency)
{
    Tmatrix t;
    t.push_back(1.0, 0.5);
    t.push_back(3.0, 0.0);

    auto ev1 = t.eigenvalues();
    auto ev2 = t.eigenvalues();
    ASSERT_EQ(ev1.size(), ev2.size());
    for (size_t i = 0; i < ev1.size(); ++i) {
        EXPECT_DOUBLE_EQ(ev1[i], ev2[i]);
	}
}


TEST(TmatrixEigenvalues, ZeroDimensionalThrows)
{
    Tmatrix t;
    EXPECT_THROW(t.eigenvalues(), std::runtime_error);
}


TEST(TmatrixEigenvectors, OneDimensional)
{
    Tmatrix t;
    t.push_back(5.0, 0.0);
    auto vecs = t.eigenvectors();
    ASSERT_EQ(vecs.size(), 1u);
    ASSERT_EQ(vecs[0].size(), 1u);
    EXPECT_NEAR(std::abs(vecs[0][0]), 1.0, 1e-12);
}


TEST(TmatrixEigenvectors, Orthonormality)
{
    const unsigned int N = 5;
    Tmatrix t;
    for (unsigned int i = 0; i < N; ++i) {
        t.push_back(2.0, -1.0);
	}

    auto vecs = t.eigenvectors();
    ASSERT_EQ(vecs.size(), static_cast<size_t>(N));

    for (unsigned int i = 0; i < N; ++i) {
        for (unsigned int j = 0; j < N; ++j) {
            double dot = 0.0;
            for (unsigned int k = 0; k < N; ++k) {
                dot += vecs[i][k] * vecs[j][k];
			}
            EXPECT_NEAR(dot, (i == j) ? 1.0 : 0.0, 1e-10)
                << "  inner product of vec " << i << " and " << j;
        }
    }
}


TEST(TmatrixEigenvectors, ResidualCheck)
{
    const unsigned int N = 4;
    std::vector<double> alpha = {2.0, 3.0, 1.0, 4.0};
    std::vector<double> beta  = {0.5, 0.7, 0.3, 0.0};

    Tmatrix t(alpha, beta);
    auto [eigvals, eigvecs] = t.eig();

    ASSERT_EQ(eigvals.size(), static_cast<size_t>(N));
    ASSERT_EQ(eigvecs.size(), static_cast<size_t>(N));

    for (unsigned int i = 0; i < N; ++i) {
        const auto& v = eigvecs[i];
        double lam = eigvals[i];

        for (unsigned int k = 0; k < N; ++k) {
            double tv = alpha[k] * v[k];
            if (k > 0) {
				tv += beta[k-1] * v[k-1];
			}
            if (k < N-1) {
				tv += beta[k] * v[k+1];
			}
            EXPECT_NEAR(tv, lam * v[k], 1e-9)
                << "  eigenvector " << i << ", component " << k;
        }
    }
}


TEST(TmatrixEig, ConsistentWithEigenvalues)
{
    Tmatrix t;
    t.push_back(2.0, -1.0);
    t.push_back(3.0, -1.0);
    t.push_back(1.0,  0.0);

    auto ev_only  = t.eigenvalues();
    auto [ev, vv] = t.eig();

    expect_near_vec(ev_only, ev);
}


TEST(TmatrixEig, ZeroDimensionalThrows)
{
    Tmatrix t;
    EXPECT_THROW(t.eig(), std::runtime_error);
}
