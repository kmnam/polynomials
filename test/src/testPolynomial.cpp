#define BOOST_TEST_MODULE testPolynomial
#define BOOST_TEST_DYN_LINK
#include <iostream>
#include <Eigen/Dense>
#include <boost/test/included/unit_test.hpp>
#include "../../include/polynomial/polynomial.hpp"

/*
 * Test module for the Polynomial class.
 *
 * Authors:
 *     Kee-Myoung Nam, Department of Systems Biology, Harvard Medical School
 * Last updated:
 *     11/21/2019
 */
typedef Polynomial<double> PolyDouble;

BOOST_AUTO_TEST_CASE(testInitialize)
{
    /*
     * Check that both constructors correctly initialize. 
     */
    PolyDouble f;
    BOOST_TEST(f.degree() == 0);
    BOOST_TEST(f.coefficients().size() == 1);
    BOOST_TEST(f.coefficients()(0) == 0.0);

    VectorXd coefs(3);
    coefs << 1.0, 2.0, 3.0;
    PolyDouble g(coefs);
    BOOST_TEST(g.degree() == 2);
    BOOST_TEST(g.coefficients().size() == 3);
    BOOST_TEST(g.coefficients()(0) == 1.0);
    BOOST_TEST(g.coefficients()(1) == 2.0);
    BOOST_TEST(g.coefficients()(2) == 3.0);
}

BOOST_AUTO_TEST_CASE(testInitializeFromDouble)
{
    /*
     * Check that the constructor with double argument initializes correctly.
     */
    PolyDouble f(1.0);
    BOOST_TEST(f.degree() == 0);
    BOOST_TEST(f.coefficients().size() == 1);
    BOOST_TEST(f.coefficients()(0) == 1.0);
}

BOOST_AUTO_TEST_CASE(testEval)
{
    /*
     * Check that all three eval() methods are correct.
     */
    // Instantiate a quintic polynomial
    VectorXd coefs(6);
    coefs << 5.0, 6.0, 7.0, 8.0, 9.0, 10.0;
    PolyDouble p(coefs);

    // Evaluate the polynomial at x = 1
    double px = p.eval(1.0);
    BOOST_TEST(px == 45.0);

    // Evaluate the polynomial at z = 2 + 3i
    std::complex<double> z(2.0, 3.0);
    std::complex<double> pz = p.eval(z);
    BOOST_TEST((pz.real() == -237.0 && pz.imag() == -6876.0));

    // Evaluate the polynomial at x = 1, 2, 4
    VectorXd v(3);
    v << -1, 2, 4;
    VectorXd pv = p.eval(v);
    BOOST_TEST((pv(0) == -3.0 && pv(1) == 573.0 && pv(2) == 13197.0));
}

BOOST_AUTO_TEST_CASE(testAdd)
{
    /*
     * Check that addition is correct. 
     */
    // Case 1: Add a cubic and a quintic (no canceling)
    VectorXd coefs_p(4);
    VectorXd coefs_q(6);
    coefs_p << 1.0, 2.0, 3.0, 4.0;
    coefs_q << 5.0, 6.0, 7.0, 8.0, 9.0, 10.0;
    PolyDouble p(coefs_p), q(coefs_q);
    PolyDouble r = p + q;
    BOOST_TEST(r.degree() == 5);
    BOOST_TEST(r.coefficients().size() == 6);
    BOOST_TEST(r.coefficients()(0) == 6.0);
    BOOST_TEST(r.coefficients()(1) == 8.0);
    BOOST_TEST(r.coefficients()(2) == 10.0);
    BOOST_TEST(r.coefficients()(3) == 12.0);
    BOOST_TEST(r.coefficients()(4) == 9.0);
    BOOST_TEST(r.coefficients()(5) == 10.0);
}

BOOST_AUTO_TEST_CASE(testSub)
{
    /*
     * Check that subtraction is correct.
     */
    // Case 1: Subtract a quintic from a cubic (no canceling)
    VectorXd coefs_p(4);
    VectorXd coefs_q(6);
    coefs_p << 1.0, 2.0, 3.0, 4.0;
    coefs_q << 5.0, 6.0, 7.0, 8.0, 9.0, 10.0;
    PolyDouble p(coefs_p), q(coefs_q);
    PolyDouble r = p - q;
    BOOST_TEST(r.degree() == 5);
    BOOST_TEST(r.coefficients().size() == 6);
    BOOST_TEST(r.coefficients()(0) == -4.0);
    BOOST_TEST(r.coefficients()(1) == -4.0);
    BOOST_TEST(r.coefficients()(2) == -4.0);
    BOOST_TEST(r.coefficients()(3) == -4.0);
    BOOST_TEST(r.coefficients()(4) == -9.0);
    BOOST_TEST(r.coefficients()(5) == -10.0);
}

BOOST_AUTO_TEST_CASE(testMul)
{
    /*
     * Check that multiplication is correct.
     */
    // Case 1: Multiply a cubic by a quintic
    VectorXd coefs_p(4);
    VectorXd coefs_q(6);
    coefs_p << 1.0, 2.0, 3.0, 4.0;
    coefs_q << 5.0, 6.0, 7.0, 8.0, 9.0, 10.0;
    PolyDouble p(coefs_p), q(coefs_q);
    PolyDouble r = p * q;
    BOOST_TEST(r.degree() == 8);
    BOOST_TEST(r.coefficients().size() == 9);
    BOOST_TEST(r.coefficients()(0) == 5.0);
    BOOST_TEST(r.coefficients()(1) == 16.0);
    BOOST_TEST(r.coefficients()(2) == 34.0);
    BOOST_TEST(r.coefficients()(3) == 60.0);
    BOOST_TEST(r.coefficients()(4) == 70.0);
    BOOST_TEST(r.coefficients()(5) == 80.0);
    BOOST_TEST(r.coefficients()(6) == 79.0);
    BOOST_TEST(r.coefficients()(7) == 66.0);
    BOOST_TEST(r.coefficients()(8) == 40.0);
}

BOOST_AUTO_TEST_CASE(testQuadraticRoots)
{
    /*
     * Obtain the roots of a quadratic.
     */
    using std::abs;
    using std::sqrt;

    // Case 1: No repeated roots, both real
    VectorXd coefs_p(3);
    coefs_p << 2.0, 3.0, 1.0;    // x^2 + 3x + 2
    PolyDouble p(coefs_p);
    Matrix<std::complex<double>, Dynamic, 1> roots_p = p.roots();
    BOOST_TEST(roots_p.size() == 2);
    BOOST_TEST((abs(roots_p(0).real() + 2) < 1e-10 || abs(roots_p(0).real() + 1) < 1e-10));
    BOOST_TEST((abs(roots_p(1).real() + 2) < 1e-10 || abs(roots_p(1).real() + 1) < 1e-10));
    BOOST_TEST((abs(roots_p(0).imag()) < 1e-10 && abs(roots_p(1).imag()) < 1e-10));

    // Case 2: No repeated roots, both complex
    VectorXd coefs_q(3);
    coefs_q << 5.0, 0.0, 1.0;    // x^2 + 5
    PolyDouble q(coefs_q);
    Matrix<std::complex<double>, Dynamic, 1> roots_q = q.roots();
    BOOST_TEST(roots_q.size() == 2);
    BOOST_TEST((abs(roots_q(0).real()) < 1e-10 && abs(roots_q(1).real()) < 1e-10));
    BOOST_TEST((abs(roots_q(0).imag() + sqrt(5)) < 1e-10 || abs(roots_q(0).imag() - sqrt(5)) < 1e-10));
    BOOST_TEST((abs(roots_q(1).imag() + sqrt(5)) < 1e-10 || abs(roots_q(1).imag() - sqrt(5)) < 1e-10));
}
