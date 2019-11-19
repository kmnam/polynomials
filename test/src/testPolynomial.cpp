#define BOOST_TEST_MODULE testPolynomial
#define BOOST_TEST_DYN_LINK
#include <iostream>
#include <Eigen/Dense>
#include <boost/test/unit_test.hpp>
#include "../../include/polynomial.hpp"

/*
 * Test module for the Polynomial class.
 *
 * Authors:
 *     Kee-Myoung Nam, Department of Systems Biology, Harvard Medical School
 * Last updated:
 *     11/19/2019
 */
typedef Polynomial<double> PolyD;

BOOST_AUTO_TEST_CASE(testInitialize)
{
    PolyD f;
    BOOST_TEST(f.degree() == 0);
    BOOST_TEST(f.coefficients().size() == 1);
    BOOST_TEST(f.coefficients()(0) == 0.0);

    VectorXd coefs(3);
    coefs << 1.0, 2.0, 3.0;
    PolyD g(coefs);
    BOOST_TEST(g.degree() == 2);
    BOOST_TEST(g.coefficients().size() == 3);
    BOOST_TEST(g.coefficients()(0) == 1.0);
    BOOST_TEST(g.coefficients()(1) == 2.0);
    BOOST_TEST(g.coefficients()(2) == 3.0);
}

BOOST_AUTO_TEST_CASE(testInitializeFromDouble)
{
    PolyD f(1.0);
    BOOST_TEST(f.degree() == 0);
    BOOST_TEST(f.coefficients().size() == 1);
    BOOST_TEST(f.coefficients()(0) == 1.0);
}

BOOST_AUTO_TEST_CASE(testSum)
{
    // Add a cubic and a quintic
    VectorXd coefs_p(4);
    VectorXd coefs_q(6);
    coefs_p << 1.0, 2.0, 3.0, 4.0;
    coefs_q << 5.0, 6.0, 7.0, 8.0, 9.0, 10.0;
    PolyD p(coefs_p), q(coefs_q);
    PolyD sum = p + q;
    BOOST_TEST(sum.degree() == 5);
    BOOST_TEST(sum.coefficients().size() == 6);
    BOOST_TEST(sum.coefficients()(0) == 6.0);
    BOOST_TEST(sum.coefficients()(1) == 8.0);
    BOOST_TEST(sum.coefficients()(2) == 10.0);
    BOOST_TEST(sum.coefficients()(3) == 12.0);
    BOOST_TEST(sum.coefficients()(4) == 9.0);
    BOOST_TEST(sum.coefficients()(5) == 10.0);
}
