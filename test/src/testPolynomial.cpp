#define BOOST_TEST_MODULE testPolynomial
#define BOOST_TEST_DYN_LINK
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "../../include/polynomial/polynomial.hpp"

/*
 * Test module for the Polynomial class.
 *
 * Authors:
 *     Kee-Myoung Nam, Department of Systems Biology, Harvard Medical School
 * Last updated:
 *     12/16/2020
 */
BOOST_AUTO_TEST_CASE(testInit)
{
    /*
     * Check that all constructors correctly initialize. 
     */
    typedef number<mpfr_float_backend<20> > mpfr_20;

    // Testing default (empty) constructor
    Polynomial<20> f;
    BOOST_TEST(f.degree() == 0);
    BOOST_TEST(f.coefficients().size() == 1);
    BOOST_TEST(f.coefficients()[0] == 0.0);

    // Testing constructor with scalar argument
    Polynomial<20> g(5);
    BOOST_TEST(g.degree() == 0);
    BOOST_TEST(g.coefficients().size() == 1);
    BOOST_TEST(g.coefficients()[0] == 5.0);

    // Testing constructor with std::vector argument
    std::vector<mpfr_20> coefs = {1.0, 2.0, 3.0};
    Polynomial<20> h(coefs);
    BOOST_TEST(h.degree() == 2);
    std::vector<mpfr_20> h_coefs = h.coefficients();
    BOOST_TEST(h_coefs.size() == 3);
    BOOST_TEST(h_coefs[0] == 1.0);
    BOOST_TEST(h_coefs[1] == 2.0);
    BOOST_TEST(h_coefs[2] == 3.0);

    // Testing constructor with double scalar argument
    double value = 17;
    Polynomial<20> k(value);
    BOOST_TEST(k.degree() == 0);
    BOOST_TEST(k.coefficients().size() == 1);
    BOOST_TEST(k.coefficients()[0] == 17.0);

    // Testing constructor with std::vector<double> argument
    std::vector<double> coefs_dbl = {71.0, 42.0, 99.0};
    Polynomial<20> m(coefs_dbl);
    BOOST_TEST(m.degree() == 2);
    std::vector<mpfr_20> m_coefs = m.coefficients();
    BOOST_TEST(m_coefs.size() == 3);
    BOOST_TEST(m_coefs[0] == 71.0);
    BOOST_TEST(m_coefs[1] == 42.0);
    BOOST_TEST(m_coefs[2] == 99.0);
}

BOOST_AUTO_TEST_CASE(testEval)
{
    /*
     * Test all four eval() methods.
     */
    typedef number<mpfr_float_backend<20> > mpfr_20;
    typedef number<mpc_complex_backend<20> > mpc_20;
    typedef number<mpfr_float_backend<30> > mpfr_30;
    typedef number<mpc_complex_backend<30> > mpc_30;
    typedef number<mpfr_float_backend<40> > mpfr_40;
    typedef number<mpc_complex_backend<40> > mpc_40;

    // Evaluate with constant precision (N = M = P = 20)
    std::vector<mpfr_20> coefs = {1.0, 2.0, 3.0};
    Polynomial<20> f(coefs);
    mpc_20 z(5.0, 0.0);
    BOOST_TEST(f.eval(z) == 86.0);
    BOOST_TEST(f.eval<20>(z) == 86.0);
    BOOST_TEST((f.eval<20, 20>(z)) == 86.0);
    mpfr_20 x = 5.0;
    BOOST_TEST(f.eval(x) == 86.0);
    BOOST_TEST(f.eval<20>(x) == 86.0);
    BOOST_TEST((f.eval<20, 20>(x)) == 86.0);

    // Evaluate with input value of precision = 30
    mpc_30 a(5.0, 0.0);
    BOOST_TEST(f.eval<30>(a) == 86.0);
    BOOST_TEST((f.eval<30, 20>(a)) == 86.0);
    mpfr_30 y = 5.0;
    BOOST_TEST(f.eval<30>(y) == 86.0);
    BOOST_TEST((f.eval<30, 20>(y)) == 86.0);

    // Evaluate with input value of precision = 30 and output value of precision = 40
    BOOST_TEST((f.eval<30, 40>(a)) == 86.0);
    BOOST_TEST((f.eval<30, 40>(y)) == 86.0);
}
