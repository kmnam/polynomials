#define BOOST_TEST_MODULE testPolynomial
#define BOOST_TEST_DYN_LINK
#include <iostream>
#include <string>
#include <sstream>
#include <limits>
#include <utility>
#include <tuple>
#include <boost/random.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpc.hpp>
#include <boost/test/included/unit_test.hpp>
#include "../../include/polynomial/polynomial.hpp"

/*
 * Test module for the Polynomial class.
 *
 * Authors:
 *     Kee-Myoung Nam, Department of Systems Biology, Harvard Medical School
 * Last updated:
 *     5/8/2021
 */
using boost::multiprecision::number; 
using boost::multiprecision::mpfr_float_backend; 
using boost::multiprecision::mpc_complex_backend; 

template <unsigned N>
std::tuple<bool, std::string, int> getNumberAsString(number<mpfr_float_backend<N> > num)
{
    /*
     * Return a maximally precise string representation of the given number. 
     */
    std::stringstream ss;

    // Feed number into stringstream with maximum precision
    ss << std::setprecision(std::numeric_limits<number<mpfr_float_backend<N> > >::max_digits10)
       << std::scientific;
    ss << num;

    // Get output string and parse sign, mantissa, and exponent 
    std::string s = ss.str();
    bool sign;       // Returns true for zero
    int exponent;
    std::string mantissa; 
    if (s.rfind("-", 0) == 0)    // Negative number 
    {
        sign = false;  
        mantissa = s.substr(1, s.find("e"));
        s.erase(0, s.find("e") + 1); 
        exponent = std::stoi(s); 
    }
    else                         // Positive number or zero 
    {
        sign = true;
        mantissa = s.substr(0, s.find("e")); 
        s.erase(0, s.find("e") + 1); 
        exponent = std::stoi(s); 
    }

    return std::make_tuple(sign, mantissa, exponent);
}

template <unsigned N>
std::vector<number<mpfr_float_backend<N> > > logRandomCoefs(const unsigned deg,
                                                            boost::random::mt19937& rng,
                                                            boost::random::uniform_real_distribution<>& dist)
{
    /*
     * Return a vector of randomly sampled coefficients for a polynomial 
     * of given degree, from a log-uniform distribution with the given
     * min/max exponents.
     *
     * The numbers are sampled first as doubles, then cast to the given 
     * precision (N).  
     */
    std::vector<number<mpfr_float_backend<N> > > coefs;
    for (int i = 0; i < deg + 1; ++i)
    {
        double c = dist(rng);
        number<mpfr_float_backend<N> > log_coef(c);
        coefs.push_back(boost::multiprecision::pow(10.0, log_coef));
    }

    return coefs; 
}

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

    // Testing precision of stored coefficient 
    std::tuple<bool, std::string, int> coef = getNumberAsString<20>(g.coefficients()[0]);
    BOOST_TEST(std::get<0>(coef)); 
    BOOST_TEST(std::get<1>(coef).size() > 20);
    BOOST_TEST(std::get<1>(coef).rfind("5.0000000000000000000", 0) == 0);
    BOOST_TEST(std::get<2>(coef) == 0);

    // Testing constructor with std::vector argument
    std::vector<mpfr_20> coefs = {1.0, 2.0, 3.0};
    Polynomial<20> h(coefs);
    BOOST_TEST(h.degree() == 2);
    std::vector<mpfr_20> h_coefs = h.coefficients();
    BOOST_TEST(h_coefs.size() == 3);
    BOOST_TEST(h_coefs[0] == 1.0);
    BOOST_TEST(h_coefs[1] == 2.0);
    BOOST_TEST(h_coefs[2] == 3.0);

    // Testing precision of stored coefficients
    for (int i = 0; i < 3; ++i)
    {
        coef = getNumberAsString<20>(h_coefs[i]);
        BOOST_TEST(std::get<0>(coef));
        BOOST_TEST(std::get<1>(coef).size() > 20);
        if (i == 0)      BOOST_TEST(std::get<1>(coef).rfind("1.0000000000000000000", 0) == 0);
        else if (i == 1) BOOST_TEST(std::get<1>(coef).rfind("2.0000000000000000000", 0) == 0);
        else             BOOST_TEST(std::get<1>(coef).rfind("3.0000000000000000000", 0) == 0);
        BOOST_TEST(std::get<2>(coef) == 0);
    }

    // Testing constructor with double scalar argument
    double value = 17;
    Polynomial<20> k(value);
    BOOST_TEST(k.degree() == 0);
    BOOST_TEST(k.coefficients().size() == 1);
    BOOST_TEST(k.coefficients()[0] == 17.0);

    // Testing precision of stored coefficient 
    coef = getNumberAsString<20>(k.coefficients()[0]);
    BOOST_TEST(std::get<0>(coef));
    BOOST_TEST(std::get<1>(coef).size() > 20);
    BOOST_TEST(std::get<1>(coef).rfind("1.7000000000000000000", 0) == 0);
    BOOST_TEST(std::get<2>(coef) == 1);

    // Testing constructor with std::vector<double> argument
    std::vector<double> coefs_dbl = {71.0, 42.0, 99.0};
    Polynomial<20> m(coefs_dbl);
    BOOST_TEST(m.degree() == 2);
    std::vector<mpfr_20> m_coefs = m.coefficients();
    BOOST_TEST(m_coefs.size() == 3);
    BOOST_TEST(m_coefs[0] == 71.0);
    BOOST_TEST(m_coefs[1] == 42.0);
    BOOST_TEST(m_coefs[2] == 99.0);

    // Testing precision of stored coefficients
    for (int i = 0; i < 3; ++i)
    {
        coef = getNumberAsString<20>(m_coefs[i]);
        BOOST_TEST(std::get<0>(coef));
        BOOST_TEST(std::get<1>(coef).size() > 20);
        if      (i == 0) BOOST_TEST(std::get<1>(coef).rfind("7.1000000000000000000", 0) == 0);
        else if (i == 1) BOOST_TEST(std::get<1>(coef).rfind("4.2000000000000000000", 0) == 0);
        else             BOOST_TEST(std::get<1>(coef).rfind("9.9000000000000000000", 0) == 0);
        BOOST_TEST(std::get<2>(coef) == 1);
    }
}

BOOST_AUTO_TEST_CASE(testEval)
{
    /*
     * Test all scalar eval() methods.
     */
    typedef number<mpfr_float_backend<20> > mpfr_20;
    typedef number<mpc_complex_backend<20> > mpc_20;
    typedef number<mpfr_float_backend<30> > mpfr_30;
    typedef number<mpc_complex_backend<30> > mpc_30;
    typedef number<mpfr_float_backend<40> > mpfr_40;
    typedef number<mpc_complex_backend<40> > mpc_40;

    // Evaluate with constant precision (N = M = P = 20) and complex input
    // and complex output 
    std::vector<mpfr_20> coefs = {1.0, 2.0, 3.0};
    Polynomial<20> f(coefs);
    mpc_20 z(5.0, 0.0);
    mpc_20 z0 = f.eval(z);
    std::tuple<bool, std::string, int> zc0 = getNumberAsString<20>(z0.real());
    mpc_20 z1 = f.eval<20, 20>(z);
    std::tuple<bool, std::string, int> zc1 = getNumberAsString<20>(z1.real());
    BOOST_TEST(z0 == 86.0);
    BOOST_TEST(std::get<0>(zc0));
    BOOST_TEST(std::get<1>(zc0).size() > 20);
    BOOST_TEST(std::get<1>(zc0).rfind("8.6000000000000000000", 0) == 0);
    BOOST_TEST(std::get<2>(zc0) == 1);
    BOOST_TEST(z1 == 86.0);
    BOOST_TEST(std::get<0>(zc1));
    BOOST_TEST(std::get<1>(zc1).size() > 20);
    BOOST_TEST(std::get<1>(zc1).rfind("8.6000000000000000000", 0) == 0);
    BOOST_TEST(std::get<2>(zc1) == 1);

    // Evaluate with constant precision (20) and std::complex<double> input
    std::complex<double> zd(5.0, 0.0);
    mpc_20 zd0 = f.eval(zd);
    std::tuple<bool, std::string, int> zdc0 = getNumberAsString<20>(zd0.real());
    mpc_20 zd1 = f.eval<20>(zd);
    std::tuple<bool, std::string, int> zdc1 = getNumberAsString<20>(zd1.real());
    BOOST_TEST(zd0 == 86.0);
    BOOST_TEST(std::get<0>(zdc0));
    BOOST_TEST(std::get<1>(zdc0).size() > 20);
    BOOST_TEST(std::get<1>(zdc0).rfind("8.6000000000000000000", 0) == 0);
    BOOST_TEST(std::get<2>(zdc0) == 1);
    BOOST_TEST(zd1 == 86.0);
    BOOST_TEST(std::get<0>(zdc1)); 
    BOOST_TEST(std::get<1>(zdc1).size() > 20);
    BOOST_TEST(std::get<1>(zdc1).rfind("8.6000000000000000000", 0) == 0);
    BOOST_TEST(std::get<2>(zdc1) == 1);

    // Evaluate with constant precision (20) and real input and real output 
    mpfr_20 x = 5.0;
    mpfr_20 x0 = f.eval(x);
    std::tuple<bool, std::string, int> xc0 = getNumberAsString<20>(x0);
    mpfr_20 x1 = f.eval<20, 20>(x);
    std::tuple<bool, std::string, int> xc1 = getNumberAsString<20>(x1);
    BOOST_TEST(x0 == 86.0);
    BOOST_TEST(std::get<0>(xc0)); 
    BOOST_TEST(std::get<1>(xc0).size() > 20);
    BOOST_TEST(std::get<1>(xc0).rfind("8.6000000000000000000", 0) == 0);
    BOOST_TEST(std::get<2>(xc0) == 1);
    BOOST_TEST(x1 == 86.0);
    BOOST_TEST(std::get<0>(xc1)); 
    BOOST_TEST(std::get<1>(xc1).size() > 20);
    BOOST_TEST(std::get<1>(xc1).rfind("8.6000000000000000000", 0) == 0);
    BOOST_TEST(std::get<2>(xc1) == 1);

    // Evaluate with constant precision (20) and double input
    double xd = 5.0;
    mpfr_20 xd0 = f.eval(xd);
    std::tuple<bool, std::string, int> xdc0 = getNumberAsString<20>(xd0.real());
    mpfr_20 xd1 = f.eval<20>(xd);
    std::tuple<bool, std::string, int> xdc1 = getNumberAsString<20>(xd1.real());
    BOOST_TEST(xd0 == 86.0);
    BOOST_TEST(std::get<0>(xdc0)); 
    BOOST_TEST(std::get<1>(xdc0).size() > 20);
    BOOST_TEST(std::get<1>(xdc0).rfind("8.6000000000000000000", 0) == 0);
    BOOST_TEST(std::get<2>(xdc0) == 1);
    BOOST_TEST(xd1 == 86.0);
    BOOST_TEST(std::get<0>(xdc1)); 
    BOOST_TEST(std::get<1>(xdc1).size() > 20);
    BOOST_TEST(std::get<1>(xdc1).rfind("8.6000000000000000000", 0) == 0);
    BOOST_TEST(std::get<2>(xdc1) == 1);

    // Evaluate with input precision = 30 and output precision = 20
    mpc_30 a(5.0, 0.0);
    mpc_20 a0 = f.eval<30, 20>(a);
    std::tuple<bool, std::string, int> ac0 = getNumberAsString<20>(a0.real());
    mpfr_30 y = 5.0;
    mpfr_20 y0 = f.eval<30, 20>(y);
    std::tuple<bool, std::string, int> yc0 = getNumberAsString<20>(y0);
    BOOST_TEST(a0 == 86.0);
    BOOST_TEST(std::get<0>(ac0));
    BOOST_TEST(std::get<1>(ac0).size() > 20);
    BOOST_TEST(std::get<1>(ac0).rfind("8.6000000000000000000", 0) == 0);
    BOOST_TEST(std::get<2>(ac0) == 1);
    BOOST_TEST(y0 == 86.0);
    BOOST_TEST(std::get<0>(yc0));
    BOOST_TEST(std::get<1>(yc0).size() > 20);
    BOOST_TEST(std::get<1>(yc0).rfind("8.6000000000000000000", 0) == 0);
    BOOST_TEST(std::get<2>(yc0) == 1);

    // Evaluate with input precision = 30 and output precision = 40
    mpc_40 a2 = f.eval<30, 40>(a);
    std::tuple<bool, std::string, int> ac2 = getNumberAsString<40>(a2.real());
    mpfr_40 y2 = f.eval<30, 40>(y);
    std::tuple<bool, std::string, int> yc2 = getNumberAsString<40>(y2);
    BOOST_TEST(a2 == 86.0);
    BOOST_TEST(std::get<0>(ac2)); 
    BOOST_TEST(std::get<1>(ac2).size() > 40);
    BOOST_TEST(std::get<1>(ac2).rfind("8.600000000000000000000000000000000000000", 0) == 0);
    BOOST_TEST(std::get<2>(ac2) == 1);
    BOOST_TEST(y2 == 86.0);
    BOOST_TEST(std::get<0>(yc2)); 
    BOOST_TEST(std::get<1>(yc2).size() > 40);
    BOOST_TEST(std::get<1>(yc2).rfind("8.600000000000000000000000000000000000000", 0) == 0);
    BOOST_TEST(std::get<2>(yc2) == 1);
}

BOOST_AUTO_TEST_CASE(testOperators)
{
    /*
     * Test the addition, subtraction, and multiplication operators.  
     */
    typedef number<mpfr_float_backend<20> > mpfr_20;
    typedef number<mpfr_float_backend<30> > mpfr_30;
    typedef number<mpfr_float_backend<40> > mpfr_40;

    // Define four polynomials of varying precisions and degrees 
    std::vector<mpfr_20> p_coefs;           // 41.5 - 74.6*x + 19.6*x^2
    p_coefs.push_back(mpfr_20("41.5"));
    p_coefs.push_back(mpfr_20("-74.6"));
    p_coefs.push_back(mpfr_20("19.6"));
    Polynomial<20> p(p_coefs);
    for (int i = 0; i < 3; ++i)
    {
        std::tuple<bool, std::string, int> coef = getNumberAsString<20>(p.coefficients()[i]);
        if (i == 0)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("4.1500000000000000000", 0) == 0 || std::get<1>(coef).rfind("4.1499999999999999999", 0) == 0));
        }
        else if (i == 1)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("7.4600000000000000000", 0) == 0 || std::get<1>(coef).rfind("7.4599999999999999999", 0) == 0));
        }
        else if (i == 2)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.9600000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.9599999999999999999", 0) == 0));
        }
        BOOST_TEST(std::get<2>(coef) == 1); 
    }
    std::vector<mpfr_30> q_coefs;           // 55.9 + 58.5*x - 96.5*x^2 + 83.4*x^3
    q_coefs.push_back(mpfr_30("55.9"));
    q_coefs.push_back(mpfr_30("58.5"));
    q_coefs.push_back(mpfr_30("-96.5"));
    q_coefs.push_back(mpfr_30("83.4"));
    Polynomial<30> q(q_coefs);
    for (int i = 0; i < 4; ++i)
    {
        std::tuple<bool, std::string, int> coef = getNumberAsString<30>(q.coefficients()[i]);
        if (i == 0)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("5.59000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("5.589999999999999999999999999999", 0) == 0));
        }
        else if (i == 1)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("5.85000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("5.849999999999999999999999999999", 0) == 0));
        }
        else if (i == 2)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("9.65000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("9.649999999999999999999999999999", 0) == 0));
        }
        else if (i == 3)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("8.34000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("8.339999999999999999999999999999", 0) == 0));
        }
        BOOST_TEST(std::get<2>(coef) == 1);
    }
    std::vector<mpfr_40> r_coefs;           // 12.9 + 12.1*x + 0.169*x^2 - 5.72*x^3 + 30.0*x^4 + 56.4*x^5
    r_coefs.push_back(mpfr_40("12.9"));
    r_coefs.push_back(mpfr_40("12.1"));
    r_coefs.push_back(mpfr_40("0.169"));
    r_coefs.push_back(mpfr_40("-5.72"));
    r_coefs.push_back(mpfr_40("30.0"));
    r_coefs.push_back(mpfr_40("56.4")); 
    Polynomial<40> r(r_coefs);
    for (int i = 0; i < 6; ++i)
    {
        std::tuple<bool, std::string, int> coef = getNumberAsString<40>(r.coefficients()[i]);
        if (i == 0)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.290000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.289999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1); 
        }
        else if (i == 1)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.210000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.209999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 2)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.690000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.689999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == -1); 
        }
        else if (i == 3)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("5.720000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("5.719999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 0);
        }
        else if (i == 4)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("3.000000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("2.999999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1); 
        }
        else if (i == 5)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("5.640000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("5.639999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1); 
        }
    }
    std::vector<mpfr_40> s_coefs;           // 5.25 + 5.88*x - 8.66*x^2 + 74.7*x^3 - 22.5*x^4
    s_coefs.push_back(mpfr_40("5.25"));
    s_coefs.push_back(mpfr_40("5.88"));
    s_coefs.push_back(mpfr_40("-8.66"));
    s_coefs.push_back(mpfr_40("74.7"));
    s_coefs.push_back(mpfr_40("-22.5"));
    Polynomial<40> s(s_coefs);
    for (int i = 0; i < 5; ++i)
    {
        std::tuple<bool, std::string, int> coef = getNumberAsString<40>(s.coefficients()[i]);
        if (i == 0)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("5.250000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("5.249999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 0); 
        }
        else if (i == 1)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("5.880000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("5.879999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 0); 
        }
        else if (i == 2)
        {
            BOOST_TEST(!std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("8.660000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("8.659999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 0); 
        }
        else if (i == 3)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("7.470000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("7.469999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 4)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("2.250000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("2.249999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
    }

    // Add p and q to get a result of precision 40
    Polynomial<40> p_plus_q = p.operator+<30, 40>(q);
    std::vector<mpfr_40> sum_coefs = p_plus_q.coefficients();
    BOOST_TEST(p_plus_q.degree() == 3); 
    BOOST_TEST(sum_coefs.size() == 4);
    for (int i = 0; i < 4; ++i)
    {
        std::tuple<bool, std::string, int> coef = getNumberAsString<40>(sum_coefs[i]);
        if (i == 0)
        {
            BOOST_TEST(std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("9.740000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("9.739999999999999999999999999999999999999", 0) == 0));
        }
        else if (i == 1)
        {
            BOOST_TEST(!std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("1.610000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.609999999999999999999999999999999999999", 0) == 0));
        }
        else if (i == 2)
        {
            BOOST_TEST(!std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("7.690000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("7.689999999999999999999999999999999999999", 0) == 0));
        }
        else if (i == 3)
        {
            BOOST_TEST(std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("8.3400000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("8.339999999999999999999999999999999999999", 0) == 0));
        }
        BOOST_TEST(std::get<2>(coef) == 1);
    }

    // Add r and s (both precision 40, output precision 40)
    Polynomial<40> r_plus_s = r + s;
    sum_coefs = r_plus_s.coefficients();
    BOOST_TEST(r_plus_s.degree() == 5);
    BOOST_TEST(sum_coefs.size() == 6);
    for (int i = 0; i < 6; ++i)
    {
        std::tuple<bool, std::string, int> coef = getNumberAsString<40>(sum_coefs[i]);
        if (i == 0)
        {     
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.815000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.814999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 1)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.798000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.797999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 2)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("8.491000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("8.490999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 0);
        }
        else if (i == 3)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("6.898000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("6.897999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 4)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("7.500000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("7.499999999999999999999999999999999999999", 0) == 0)); 
            BOOST_TEST(std::get<2>(coef) == 0);
        }
        else if (i == 5)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("5.640000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("5.639999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
    } 

    // Subtract q from p to get a result of precision 40
    Polynomial<40> p_minus_q = p.operator-<30, 40>(q);
    std::vector<mpfr_40> diff_coefs = p_minus_q.coefficients();
    BOOST_TEST(p_minus_q.degree() == 3); 
    BOOST_TEST(diff_coefs.size() == 4);
    for (int i = 0; i < 4; ++i)
    {
        std::tuple<bool, std::string, int> coef = getNumberAsString<40>(diff_coefs[i]);
        if (i == 0)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.440000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.439999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 1)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.331000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.330999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 2);
        }
        else if (i == 2)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.161000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.160999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 2);
        }
        else if (i == 3)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("8.340000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("8.339999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
    }

    // Subtract s from r (both precision 40, output precision 40)
    Polynomial<40> r_minus_s = r - s;
    diff_coefs = r_minus_s.coefficients();
    BOOST_TEST(r_minus_s.degree() == 5);
    BOOST_TEST(diff_coefs.size() == 6);
    for (int i = 0; i < 6; ++i)
    {
        std::tuple<bool, std::string, int> coef = getNumberAsString<40>(diff_coefs[i]);
        if (i == 0)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("7.650000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("7.649999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 0);
        }
        else if (i == 1)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("6.220000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("6.219999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 0);
        }
        else if (i == 2)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("8.829000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("8.828999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 0);
        }
        else if (i == 3)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("8.042000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("8.041999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 4)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("5.250000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("5.249999999999999999999999999999999999999", 0) == 0)); 
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 5)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("5.640000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("5.639999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
    } 

    // Multiply p and q to get a result of precision 40
    // Result: 2319.85 - 1742.39*x - 7273.21*x^2 + 11806.6*x^3 - 8113.04*x^4 + 1634.64*x^5
    Polynomial<40> p_times_q = p.operator*<30, 40>(q);
    std::vector<mpfr_40> prod_coefs = p_times_q.coefficients();
    BOOST_TEST(p_times_q.degree() == 5);
    BOOST_TEST(prod_coefs.size() == 6);
    for (int i = 0; i < 6; ++i)
    {
        std::tuple<bool, std::string, int> coef = getNumberAsString<40>(prod_coefs[i]);
        if (i == 0)
        {     
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("2.319850000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("2.319849999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 3);
        }
        else if (i == 1)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.742390000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.742389999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 3);
        }
        else if (i == 2)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("7.273210000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("7.273209999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 3);
        }
        else if (i == 3)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.180660000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.180659999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 4);
        }
        else if (i == 4)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("8.113040000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("8.113039999999999999999999999999999999999", 0) == 0)); 
            BOOST_TEST(std::get<2>(coef) == 3);
        }
        else if (i == 5)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.634640000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.634639999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 3);
        }
    }

    // Multiply r and s (both precision 40, output precision 40)
    // Result: 67.725 + 139.377*x - 39.67875*x^2 + 829.80772*x^3 + 736.02286*x^4
    //       + 262.4095*x^5 - 359.2545*x^6 + 1881.276*x^7 + 3538.08*x^8 - 1269*x^9
    Polynomial<40> r_times_s = r * s;
    prod_coefs = r_times_s.coefficients();
    BOOST_TEST(r_times_s.degree() == 9);
    BOOST_TEST(prod_coefs.size() == 10);
    for (int i = 0; i < 10; ++i)
    {
        std::tuple<bool, std::string, int> coef = getNumberAsString<40>(prod_coefs[i]);
        if (i == 0)
        {     
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("6.772500000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("6.772499999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 1)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.393770000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.393769999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 2);
        }
        else if (i == 2)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("3.967875000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("3.967874999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 3)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("8.298077200000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("8.298077199999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 2);
        }
        else if (i == 4)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("7.360228600000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("7.360228599999999999999999999999999999999", 0) == 0)); 
            BOOST_TEST(std::get<2>(coef) == 2);
        }
        else if (i == 5)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("2.624095000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("2.624094999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 2);
        }
        else if (i == 6)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("3.592545000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("3.592544999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 2);
        }
        else if (i == 7)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.881276000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.881275999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 3);
        }
        else if (i == 8)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("3.538080000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("3.538079999999999999999999999999999999999", 0) == 0)); 
            BOOST_TEST(std::get<2>(coef) == 3);
        }
        else if (i == 9)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.269000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.268999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 3);
        }
    }
}

BOOST_AUTO_TEST_CASE(testScalarOperators)
{
    /*
     * Test the scalar addition, subtraction, multiplication, and division
     * operators.  
     */
    typedef number<mpfr_float_backend<20> > mpfr_20;
    typedef number<mpfr_float_backend<30> > mpfr_30;
    typedef number<mpfr_float_backend<40> > mpfr_40;

    // Define p and r as before 
    std::vector<mpfr_20> p_coefs;           // 41.5 - 74.6*x + 19.6*x^2
    p_coefs.push_back(mpfr_20("41.5"));
    p_coefs.push_back(mpfr_20("-74.6"));
    p_coefs.push_back(mpfr_20("19.6"));
    Polynomial<20> p(p_coefs);
    std::vector<mpfr_40> r_coefs;           // 12.9 + 12.1*x + 0.169*x^2 - 5.72*x^3 + 30.0*x^4 + 56.4*x^5
    r_coefs.push_back(mpfr_40("12.9"));
    r_coefs.push_back(mpfr_40("12.1"));
    r_coefs.push_back(mpfr_40("0.169"));
    r_coefs.push_back(mpfr_40("-5.72"));
    r_coefs.push_back(mpfr_40("30.0"));
    r_coefs.push_back(mpfr_40("56.4")); 
    Polynomial<40> r(r_coefs);

    // Define two scalars, a and b
    mpfr_30 a("-72.6");
    mpfr_40 b("-3.73");

    // Add a to p to obtain result of precision 40
    Polynomial<40> p_plus_a = p.operator+<30, 40>(a);
    std::vector<mpfr_40> sum_coefs = p_plus_a.coefficients();
    BOOST_TEST(p_plus_a.degree() == 2);
    BOOST_TEST(sum_coefs.size() == 3);
    for (int i = 0; i < 3; ++i)
    {
        std::tuple<bool, std::string, int> coef = getNumberAsString<40>(sum_coefs[i]);
        if (i == 0)
        {
            BOOST_TEST(!std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("3.110000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("3.109999999999999999999999999999999999999", 0) == 0));
        }
        else if (i == 1)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("7.460000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("7.459999999999999999999999999999999999999", 0) == 0));
        }
        else if (i == 2)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.960000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.959999999999999999999999999999999999999", 0) == 0));
        }
        BOOST_TEST(std::get<2>(coef) == 1);
    }

    // Add b to r (both precision 40, output precision 40)
    Polynomial<40> r_plus_b = r + b;
    sum_coefs = r_plus_b.coefficients();
    BOOST_TEST(r_plus_b.degree() == 5);
    BOOST_TEST(sum_coefs.size() == 6);
    for (int i = 0; i < 6; ++i)
    {
        std::tuple<bool, std::string, int> coef = getNumberAsString<40>(sum_coefs[i]);
        if (i == 0)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("9.170000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("9.169999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 0);
        }
        else if (i == 1)
        {
            BOOST_TEST(std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("1.210000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.209999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 2)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.690000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.689999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == -1);
        }
        else if (i == 3)
        {
            BOOST_TEST(!std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("5.720000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("5.719999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 0);
        }
        else if (i == 4)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("3.000000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("2.999999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 5)
        {
            BOOST_TEST(std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("5.640000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("5.639999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
    }

    // Add r to b (both precision 40, output precision 40)
    Polynomial<40> b_plus_r = b + r;
    sum_coefs = b_plus_r.coefficients();
    BOOST_TEST(b_plus_r.degree() == 5);
    BOOST_TEST(sum_coefs.size() == 6);
    for (int i = 0; i < 6; ++i)
    {
        std::tuple<bool, std::string, int> coef = getNumberAsString<40>(sum_coefs[i]);
        if (i == 0)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("9.170000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("9.169999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 0);
        }
        else if (i == 1)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.210000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.209999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 2)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.690000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.689999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == -1);
        }
        else if (i == 3)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("5.720000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("-5.719999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 0);
        }
        else if (i == 4)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("3.000000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("2.999999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 5)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("5.640000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("5.639999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
    }

    // Subtract a from p to obtain result of precision 40
    Polynomial<40> p_minus_a = p.operator-<30, 40>(a);
    std::vector<mpfr_40> diff_coefs = p_minus_a.coefficients();
    BOOST_TEST(p_minus_a.degree() == 2);
    BOOST_TEST(diff_coefs.size() == 3);
    for (int i = 0; i < 3; ++i)
    {
        std::tuple<bool, std::string, int> coef = getNumberAsString<40>(diff_coefs[i]);
        if (i == 0)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.141000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.140999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 2);
        }
        else if (i == 1)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("7.460000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("7.459999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 2)
        {
            BOOST_TEST(std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("1.960000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.959999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
    }

    // Subtract b from r (both precision 40, output precision 40)
    Polynomial<40> r_minus_b = r - b;
    diff_coefs = r_minus_b.coefficients();
    BOOST_TEST(r_minus_b.degree() == 5);
    BOOST_TEST(diff_coefs.size() == 6);
    for (int i = 0; i < 6; ++i)
    {
        std::tuple<bool, std::string, int> coef = getNumberAsString<40>(diff_coefs[i]);
        if (i == 0)
        {
            BOOST_TEST(std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("1.663000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.662999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 1)
        {
            BOOST_TEST(std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("1.210000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.209999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 2)
        {
            BOOST_TEST(std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("1.690000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.689999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == -1);
        }
        else if (i == 3)
        {
            BOOST_TEST(!std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("5.720000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("5.719999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 0);
        }
        else if (i == 4)
        {
            BOOST_TEST(std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("3.000000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("2.999999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 5)
        {
            BOOST_TEST(std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("5.640000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("5.639999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
    }

    // Subtract r from b (both precision 40, output precision 40)
    // (Note that this also tests the negation operator)
    Polynomial<40> b_minus_r = b - r;
    diff_coefs = b_minus_r.coefficients();
    BOOST_TEST(b_minus_r.degree() == 5);
    BOOST_TEST(diff_coefs.size() == 6);
    for (int i = 0; i < 6; ++i)
    {
        std::tuple<bool, std::string, int> coef = getNumberAsString<40>(diff_coefs[i]);
        if (i == 0)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.663000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.662999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 1)
        {
            BOOST_TEST(!std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("1.210000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.209999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 2)
        {
            BOOST_TEST(!std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("1.690000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.689999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == -1);
        }
        else if (i == 3)
        {
            BOOST_TEST(std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("5.720000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("5.719999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 0);
        }
        else if (i == 4)
        {
            BOOST_TEST(!std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("3.000000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("2.999999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 5)
        {
            BOOST_TEST(!std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("5.640000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("5.639999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
    }

    // Multiply p by a to obtain result of precision 40
    // Result: -3012.9 + 5415.96*x - 1422.96*x^2
    Polynomial<40> p_times_a = p.operator*<30, 40>(a);
    std::vector<mpfr_40> prod_coefs = p_times_a.coefficients();
    BOOST_TEST(p_times_a.degree() == 2);
    BOOST_TEST(prod_coefs.size() == 3);
    for (int i = 0; i < 3; ++i)
    {
        std::tuple<bool, std::string, int> coef = getNumberAsString<40>(prod_coefs[i]);
        if (i == 0)
        {
            BOOST_TEST(!std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("3.012900000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("3.012899999999999999999999999999999999999", 0) == 0));
        }
        else if (i == 1)
        {
            BOOST_TEST(std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("5.415960000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("5.415959999999999999999999999999999999999", 0) == 0));
        }
        else if (i == 2)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.422960000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.422959999999999999999999999999999999999", 0) == 0));
        }
        BOOST_TEST(std::get<2>(coef) == 3);
    }

    // Multiply r by b (both precision 40, output precision 40)
    // Result: -48.117 - 45.133*x - 0.63037*x^2 + 21.3356*x^3 - 111.9*x^4 - 210.372*x^5
    Polynomial<40> r_times_b = r * b;
    prod_coefs = r_times_b.coefficients();
    BOOST_TEST(r_times_b.degree() == 5);
    BOOST_TEST(prod_coefs.size() == 6);
    for (int i = 0; i < 6; ++i)
    {
        std::tuple<bool, std::string, int> coef = getNumberAsString<40>(prod_coefs[i]);
        if (i == 0)
        {
            BOOST_TEST(!std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("4.811700000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("4.811699999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 1)
        {
            BOOST_TEST(!std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("4.513300000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("4.513299999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 2)
        {
            BOOST_TEST(!std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("6.303700000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("6.303699999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == -1);
        }
        else if (i == 3)
        {
            BOOST_TEST(std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("2.133560000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("2.133559999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 4)
        {
            BOOST_TEST(!std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("1.119000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.118999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 2);
        }
        else if (i == 5)
        {
            BOOST_TEST(!std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("2.103720000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("2.103719999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 2);
        }
    }

    // Multiply b by r (both precision 40, output precision 40)
    // Result: -48.117 - 45.133*x - 0.63037*x^2 + 21.3356*x^3 - 111.9*x^4 - 210.372*x^5
    Polynomial<40> b_times_r = b * r;
    prod_coefs = b_times_r.coefficients();
    BOOST_TEST(b_times_r.degree() == 5);
    BOOST_TEST(prod_coefs.size() == 6);
    for (int i = 0; i < 6; ++i)
    {
        std::tuple<bool, std::string, int> coef = getNumberAsString<40>(prod_coefs[i]);
        if (i == 0)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("4.811700000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("4.811699999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 1)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("4.513300000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("4.513299999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 2)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("6.303700000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("6.303699999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == -1);
        }
        else if (i == 3)
        {
            BOOST_TEST(std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("2.133560000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("2.133559999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
        else if (i == 4)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("1.119000000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("1.118999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 2);
        }
        else if (i == 5)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("2.103720000000000000000000000000000000000", 0) == 0 || std::get<1>(coef).rfind("2.103719999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 2);
        }
    }

    // Divide p by a to obtain result of precision 40
    // Result: -0.5716253443526170798898071625344352617080
    //         + 1.027548209366391184573002754820936639118*x
    //         - 0.2699724517906336088154269972451790633609*x^2
    Polynomial<40> p_div_a = p.operator/<30, 40>(a);
    std::vector<mpfr_40> quot_coefs = p_div_a.coefficients();
    BOOST_TEST(p_div_a.degree() == 2);
    BOOST_TEST(quot_coefs.size() == 3);
    for (int i = 0; i < 3; ++i)
    {
        std::tuple<bool, std::string, int> coef = getNumberAsString<40>(quot_coefs[i]);
        if (i == 0)
        {
            BOOST_TEST(!std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("5.716253443526170798898071625344352617080", 0) == 0 || std::get<1>(coef).rfind("5.716253443526170798898071625344352617079", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == -1);
        }
        else if (i == 1)
        {
            BOOST_TEST(std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("1.027548209366391184573002754820936639118", 0) == 0 || std::get<1>(coef).rfind("1.027548209366391184573002754820936639117", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 0);
        }
        else if (i == 2)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("2.699724517906336088154269972451790633609", 0) == 0 || std::get<1>(coef).rfind("2.699724517906336088154269972451790633608", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == -1);
        }
    }

    // Divide r by b (both precision 40, output precision 40)
    // Result: -3.458445040214477211796246648793565683646
    //         - 3.243967828418230563002680965147453083110*x
    //         - 0.04530831099195710455764075067024128686327*x^2
    //         + 1.533512064343163538873994638069705093834*x^3
    //         - 8.042895442359249329758713136729222520107*x^4
    //         - 15.12064343163538873994638069705093833780*x^5
    Polynomial<40> r_div_b = r / b;
    quot_coefs = r_div_b.coefficients();
    BOOST_TEST(r_div_b.degree() == 5);
    BOOST_TEST(quot_coefs.size() == 6);
    for (int i = 0; i < 6; ++i)
    {
        std::tuple<bool, std::string, int> coef = getNumberAsString<40>(quot_coefs[i]);
        if (i == 0)
        {
            BOOST_TEST(!std::get<0>(coef));
            BOOST_TEST((std::get<1>(coef).rfind("3.458445040214477211796246648793565683646", 0) == 0 || std::get<1>(coef).rfind("3.458445040214477211796246648793565683645", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 0);
        }
        else if (i == 1)
        {
            BOOST_TEST(!std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("3.243967828418230563002680965147453083110", 0) == 0 || std::get<1>(coef).rfind("3.243967828418230563002680965147453083109", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 0);
        }
        else if (i == 2)
        {
            BOOST_TEST(!std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("4.530831099195710455764075067024128686327", 0) == 0 || std::get<1>(coef).rfind("4.530831099195710455764075067024128686326", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == -2);
        }
        else if (i == 3)
        {
            BOOST_TEST(std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("1.533512064343163538873994638069705093834", 0) == 0 || std::get<1>(coef).rfind("1.533512064343163538873994638069705093833", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 0);
        }
        else if (i == 4)
        {
            BOOST_TEST(!std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("8.042895442359249329758713136729222520107", 0) == 0 || std::get<1>(coef).rfind("8.042895442359249329758713136729222520106", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 0);
        }
        else if (i == 5)
        {
            BOOST_TEST(!std::get<0>(coef)); 
            BOOST_TEST((std::get<1>(coef).rfind("1.512064343163538873994638069705093833780", 0) == 0 || std::get<1>(coef).rfind("1.512064343163538873994638069705093833779", 0) == 0));
            BOOST_TEST(std::get<2>(coef) == 1);
        }
    }
}

BOOST_AUTO_TEST_CASE(testSolveQuadratic)
{
    /*
     * Get the roots of a quadratic using the quadratic formula. 
     */
    typedef number<mpfr_float_backend<20> >  mpfr_20;
    typedef number<mpc_complex_backend<20> > mpc_20;
    typedef number<mpfr_float_backend<30> >  mpfr_30;
    typedef number<mpc_complex_backend<30> > mpc_30;

    // Define p as before 
    std::vector<mpfr_20> p_coefs;           // 41.5 - 74.6*x + 19.6*x^2
    p_coefs.push_back(mpfr_20("41.5"));
    p_coefs.push_back(mpfr_20("-74.6"));
    p_coefs.push_back(mpfr_20("19.6"));
    Polynomial<20> p(p_coefs);
    unsigned max_iter = 100; 
    double tol = std::numeric_limits<double>::epsilon(); 
    boost::random::mt19937 rng;
    boost::random::uniform_real_distribution<> dist(0, 1);

    // Get roots with same precision (20) as coefficients
    // Solutions: 0.67656414525108398788 and 3.1295583037285078489
    std::pair<std::vector<mpc_20>, bool> roots = p.roots(max_iter, tol, tol, rng, dist);
    std::sort(roots.first.begin(), roots.first.end(), [](mpc_20 a, mpc_20 b){ return (a.real() < b.real()); });
    for (int i = 0; i < 2; ++i)
    {
        BOOST_TEST(roots.first[i].imag() == 0.0);
        std::tuple<bool, std::string, int> root = getNumberAsString<20>(roots.first[i].real());
        if (i == 0)
        {
            BOOST_TEST(std::get<0>(root)); 
            BOOST_TEST((std::get<1>(root).rfind("6.7656414525108398788", 0) == 0 || std::get<1>(root).rfind("6.7656414525108398787", 0) == 0));
            BOOST_TEST(std::get<2>(root) == -1);
        }
        else if (i == 1)
        {
            BOOST_TEST(std::get<0>(root));
            BOOST_TEST((std::get<1>(root).rfind("3.1295583037285078489", 0) == 0 || std::get<1>(root).rfind("3.1295583037285078488", 0) == 0));
            BOOST_TEST(std::get<2>(root) == 0);
        }
    }

    // Get roots with precision 30
    // Solutions: 0.676564145251083987876879502160 and 3.12955830372850784885781437539
    std::pair<std::vector<mpc_30>, bool> roots_30 = p.roots<30>(max_iter, tol, tol, rng, dist);
    std::sort(roots_30.first.begin(), roots_30.first.end(), [](mpc_30 a, mpc_30 b){ return (a.real() < b.real()); });
    for (int i = 0; i < 2; ++i)
    {
        BOOST_TEST(roots_30.first[i].imag() == 0.0);
        std::tuple<bool, std::string, int> root = getNumberAsString<30>(roots_30.first[i].real());
        if (i == 0)
        {
            BOOST_TEST(std::get<0>(root));
            BOOST_TEST((std::get<1>(root).rfind("6.76564145251083987876879502160", 0) == 0 || std::get<1>(root).rfind("6.76564145251083987876879502159", 0) == 0));
            BOOST_TEST(std::get<2>(root) == -1);
        }
        else if (i == 1)
        {
            BOOST_TEST(std::get<0>(root));
            BOOST_TEST((std::get<1>(root).rfind("3.12955830372850784885781437539", 0) == 0 || std::get<1>(root).rfind("3.12955830372850784885781437538", 0) == 0));
            BOOST_TEST(std::get<2>(root) == 0);
        }
    }
}

BOOST_AUTO_TEST_CASE(testSolveCubic)
{
    /*
     * Get the roots of a cubic using the cubic formula. 
     */
    typedef number<mpfr_float_backend<30> >  mpfr_30;
    typedef number<mpc_complex_backend<30> > mpc_30;
    typedef number<mpfr_float_backend<40> >  mpfr_40;
    typedef number<mpc_complex_backend<40> > mpc_40;

    // Define q as before 
    std::vector<mpfr_30> q_coefs;           // 55.9 + 58.5*x - 96.5*x^2 + 83.4*x^3
    q_coefs.push_back(mpfr_30("55.9"));
    q_coefs.push_back(mpfr_30("58.5"));
    q_coefs.push_back(mpfr_30("-96.5"));
    q_coefs.push_back(mpfr_30("83.4"));
    Polynomial<30> q(q_coefs);
    unsigned max_iter = 100; 
    double tol = std::numeric_limits<double>::epsilon(); 
    boost::random::mt19937 rng;
    boost::random::uniform_real_distribution<> dist(0, 1);

    // Get roots with same precision (29, accounting for error due to computations)
    // Solutions: -0.46225599097010167068346290934,
    //            0.80966516574883980416667150263 - 0.89130596059275302403974487956*i,
    //            0.80966516574883980416667150263 + 0.89130596059275302403974487956*i
    std::pair<std::vector<mpc_30>, bool> roots = q.roots(max_iter, tol, tol, rng, dist);
    std::sort(roots.first.begin(), roots.first.end(), [](mpc_30 a, mpc_30 b)
        {
            if (boost::multiprecision::abs(a.real() - b.real()) < std::numeric_limits<mpfr_30>::epsilon())
                return (a.imag() < b.imag());
            else
                return (a.real() < b.real());
        }
    );
    for (int i = 0; i < 3; ++i)
    {
        std::tuple<bool, std::string, int> root_real = getNumberAsString<30>(roots.first[i].real());
        std::tuple<bool, std::string, int> root_imag = getNumberAsString<30>(roots.first[i].imag());
        if (i == 0)
        {
            BOOST_TEST(!std::get<0>(root_real)); 
            BOOST_TEST((std::get<1>(root_real).rfind("4.6225599097010167068346290934", 0) == 0 || std::get<1>(root_real).rfind("4.6225599097010167068346290933", 0) == 0));
            BOOST_TEST(std::get<2>(root_real) == -1);
            // Imaginary part could technically have either sign in this case
            BOOST_TEST((std::get<1>(root_imag).rfind("0.0000000000000000000000000000", 0) == 0 || std::get<1>(root_imag).rfind("-0.0000000000000000000000000000", 0) == 0));
            BOOST_TEST(std::get<2>(root_imag) == 0);
        }
        else if (i == 1)
        {
            BOOST_TEST(std::get<0>(root_real));
            BOOST_TEST((std::get<1>(root_real).rfind("8.0966516574883980416667150263", 0) == 0 || std::get<1>(root_real).rfind("8.0966516574883980416667150262", 0) == 0));
            BOOST_TEST(std::get<2>(root_real) == -1);
            BOOST_TEST(!std::get<0>(root_imag)); 
            BOOST_TEST((std::get<1>(root_imag).rfind("8.9130596059275302403974487956", 0) == 0 || std::get<1>(root_imag).rfind("8.9130596059275302403974487955", 0) == 0));
            BOOST_TEST(std::get<2>(root_imag) == -1);
        }
        else if (i == 2)
        {
            BOOST_TEST(std::get<0>(root_real)); 
            BOOST_TEST((std::get<1>(root_real).rfind("8.0966516574883980416667150263", 0) == 0 || std::get<1>(root_real).rfind("8.0966516574883980416667150262", 0) == 0));
            BOOST_TEST(std::get<2>(root_real) == -1);
            BOOST_TEST(std::get<0>(root_imag)); 
            BOOST_TEST((std::get<1>(root_imag).rfind("8.9130596059275302403974487956", 0) == 0 || std::get<1>(root_imag).rfind("8.9130596059275302403974487955", 0) == 0));
            BOOST_TEST(std::get<2>(root_imag) == -1);
        }
    }

    // Get roots with precision 40 (39, accounting for error due to computations)
    // Solutions: -0.462255990970101670683462909334950520555,
    //            0.809665165748839804166671502629105955721 - 0.891305960592753024039744879557940668782*i,
    //            0.809665165748839804166671502629105955721 + 0.891305960592753024039744879557940668782*i
    std::pair<std::vector<mpc_40>, bool> roots_40 = q.roots<40>(max_iter, tol, tol, rng, dist);
    std::sort(roots_40.first.begin(), roots_40.first.end(), [](mpc_40 a, mpc_40 b)
        {
            if (boost::multiprecision::abs(a.real() - b.real()) < std::numeric_limits<mpfr_40>::epsilon())
                return (a.imag() < b.imag());
            else
                return (a.real() < b.real());
        }
    );
    for (int i = 0; i < 3; ++i)
    {
        std::tuple<bool, std::string, int> root_real = getNumberAsString<40>(roots_40.first[i].real());
        std::tuple<bool, std::string, int> root_imag = getNumberAsString<40>(roots_40.first[i].imag());
        if (i == 0)
        {
            BOOST_TEST(!std::get<0>(root_real));
            BOOST_TEST((std::get<1>(root_real).rfind("4.62255990970101670683462909334950520555", 0) == 0 || std::get<1>(root_real).rfind("4.62255990970101670683462909334950520554", 0) == 0));
            BOOST_TEST(std::get<2>(root_real) == -1);
            // Imaginary part could technically have either sign in this case
            BOOST_TEST((std::get<1>(root_imag).rfind("0.00000000000000000000000000000000000000", 0) == 0 || std::get<1>(root_imag).rfind("-0.00000000000000000000000000000000000000", 0) == 0));
            BOOST_TEST(std::get<2>(root_imag) == 0);
        }
        else if (i == 1)
        {
            BOOST_TEST(std::get<0>(root_real));
            BOOST_TEST((std::get<1>(root_real).rfind("8.09665165748839804166671502629105955721", 0) == 0 || std::get<1>(root_real).rfind("8.09665165748839804166671502629105955720", 0) == 0));
            BOOST_TEST(std::get<2>(root_real) == -1);
            BOOST_TEST(!std::get<0>(root_imag));
            BOOST_TEST((std::get<1>(root_imag).rfind("8.91305960592753024039744879557940668782", 0) == 0 || std::get<1>(root_imag).rfind("8.91305960592753024039744879557940668781", 0) == 0));
            BOOST_TEST(std::get<2>(root_imag) == -1);
        }
        else if (i == 2)
        {
            BOOST_TEST(std::get<0>(root_real));
            BOOST_TEST((std::get<1>(root_real).rfind("8.09665165748839804166671502629105955721", 0) == 0 || std::get<1>(root_real).rfind("8.09665165748839804166671502629105955720", 0) == 0));
            BOOST_TEST(std::get<2>(root_real) == -1);
            BOOST_TEST(std::get<0>(root_imag));
            BOOST_TEST((std::get<1>(root_imag).rfind("8.91305960592753024039744879557940668782", 0) == 0 || std::get<1>(root_imag).rfind("8.91305960592753024039744879557940668781", 0) == 0));
            BOOST_TEST(std::get<2>(root_imag) == -1);
        }
    }
}

BOOST_AUTO_TEST_CASE(testSolveQuintic)
{
    /*
     * Get the roots of a quintic polynomial with Aberth's method. 
     */
    typedef number<mpfr_float_backend<40> >  mpfr_40;
    typedef number<mpc_complex_backend<40> > mpc_40;
    typedef number<mpfr_float_backend<50> >  mpfr_50;
    typedef number<mpc_complex_backend<50> > mpc_50;

    // Define r as before 
    std::vector<mpfr_40> r_coefs;           // 12.9 + 12.1*x + 0.169*x^2 - 5.72*x^3 + 30.0*x^4 + 56.4*x^5
    r_coefs.push_back(mpfr_40("12.9"));
    r_coefs.push_back(mpfr_40("12.1"));
    r_coefs.push_back(mpfr_40("0.169"));
    r_coefs.push_back(mpfr_40("-5.72"));
    r_coefs.push_back(mpfr_40("30.0"));
    r_coefs.push_back(mpfr_40("56.4")); 
    Polynomial<40> r(r_coefs);
    unsigned max_iter = 100; 
    double tol = std::numeric_limits<double>::epsilon(); 
    boost::random::mt19937 rng;
    boost::random::uniform_real_distribution<> dist(0, 1);

    // Get roots with same precision (40) as coefficients
    // Solutions: -0.8011551530834036092615229000431606689577,
    //            -0.4261332773189638878958458634920648066143 - 0.5703930114618632364573561806651850823911*i
    //            -0.4261332773189638878958458634920648066143 + 0.5703930114618632364573561806651850823911*i
    //            0.5607534070521550542287349730881132261995 - 0.4987237601604699794368211275750465313412*i
    //            0.5607534070521550542287349730881132261995 + 0.4987237601604699794368211275750465313412*i
    std::pair<std::vector<mpc_40>, bool> roots = r.roots(max_iter, tol, tol, rng, dist); 
    BOOST_TEST(roots.second);    // Check for convergence first, then sort by real/imaginary parts
    std::sort(roots.first.begin(), roots.first.end(), [](mpc_40 a, mpc_40 b)
        {
            if (boost::multiprecision::abs(a.real() - b.real()) < std::numeric_limits<mpfr_40>::epsilon())
                return (a.imag() < b.imag());
            else
                return (a.real() < b.real());
        }
    );
    for (int i = 0; i < 5; ++i)
    {
        std::tuple<bool, std::string, int> root_real = getNumberAsString<40>(roots.first[i].real());
        std::tuple<bool, std::string, int> root_imag = getNumberAsString<40>(roots.first[i].imag());
        
        if (i == 0)
        {
            BOOST_TEST(!std::get<0>(root_real)); 
            BOOST_TEST((std::get<1>(root_real).rfind("8.011551530834036092615229000431606689577", 0) == 0 || std::get<1>(root_real).rfind("8.011551530834036092615229000431606689576", 0) == 0));
            BOOST_TEST(std::get<2>(root_real) == -1);
            BOOST_TEST(boost::multiprecision::abs(roots.first[i].imag()) < std::numeric_limits<mpfr_40>::epsilon());
        }
        else if (i == 1)
        {
            BOOST_TEST(!std::get<0>(root_real));
            BOOST_TEST((std::get<1>(root_real).rfind("4.261332773189638878958458634920648066143", 0) == 0 || std::get<1>(root_real).rfind("4.261332773189638878958458634920648066142", 0) == 0));
            BOOST_TEST(std::get<2>(root_real) == -1);
            BOOST_TEST(!std::get<0>(root_imag));
            BOOST_TEST((std::get<1>(root_imag).rfind("5.703930114618632364573561806651850823911", 0) == 0 || std::get<1>(root_imag).rfind("5.703930114618632364573561806651850823910", 0) == 0));
            BOOST_TEST(std::get<2>(root_imag) == -1);
        }
        else if (i == 2)
        {
            BOOST_TEST(!std::get<0>(root_real));
            BOOST_TEST((std::get<1>(root_real).rfind("4.261332773189638878958458634920648066143", 0) == 0 || std::get<1>(root_real).rfind("4.261332773189638878958458634920648066142", 0) == 0));
            BOOST_TEST(std::get<2>(root_real) == -1);
            BOOST_TEST(std::get<0>(root_imag));
            BOOST_TEST((std::get<1>(root_imag).rfind("5.703930114618632364573561806651850823911", 0) == 0 || std::get<1>(root_imag).rfind("5.703930114618632364573561806651850823910", 0) == 0));
            BOOST_TEST(std::get<2>(root_imag) == -1);
        }
        else if (i == 3)
        {
            BOOST_TEST(std::get<0>(root_real));
            BOOST_TEST((std::get<1>(root_real).rfind("5.607534070521550542287349730881132261995", 0) == 0 || std::get<1>(root_real).rfind("5.607534070521550542287349730881132261994", 0) == 0));
            BOOST_TEST(std::get<2>(root_real) == -1);
            BOOST_TEST(!std::get<0>(root_imag));
            BOOST_TEST((std::get<1>(root_imag).rfind("4.987237601604699794368211275750465313412", 0) == 0 || std::get<1>(root_imag).rfind("4.987237601604699794368211275750465313411", 0) == 0));
            BOOST_TEST(std::get<2>(root_imag) == -1);
        }
        else if (i == 4)
        {
            BOOST_TEST(std::get<0>(root_real));
            BOOST_TEST((std::get<1>(root_real).rfind("5.607534070521550542287349730881132261995", 0) == 0 || std::get<1>(root_real).rfind("5.607534070521550542287349730881132261994", 0) == 0));
            BOOST_TEST(std::get<2>(root_real) == -1);
            BOOST_TEST(std::get<0>(root_imag));
            BOOST_TEST((std::get<1>(root_imag).rfind("4.987237601604699794368211275750465313412", 0) == 0 || std::get<1>(root_imag).rfind("4.987237601604699794368211275750465313411", 0) == 0));
            BOOST_TEST(std::get<2>(root_imag) == -1);
        }
    }

    // Get roots with precision 50
    // Solutions: -0.80115515308340360926152290004316066895772474042464,
    //            -0.42613327731896388789584586349206480661427590651186 - 0.57039301146186323645735618066518508239107126817095*i
    //            -0.42613327731896388789584586349206480661427590651186 + 0.57039301146186323645735618066518508239107126817095*i
    //            0.56075340705215505422873497308811322619952125544759 - 0.49872376016046997943682112757504653134116713321119*i
    //            0.56075340705215505422873497308811322619952125544759 + 0.49872376016046997943682112757504653134116713321119*i
    std::pair<std::vector<mpc_50>, bool> roots_50 = r.roots<50>(max_iter, tol, tol, rng, dist);
    BOOST_TEST(roots_50.second);    // Check for convergence first, then sort by real/imaginary parts
    std::sort(roots_50.first.begin(), roots_50.first.end(), [](mpc_50 a, mpc_50 b)
        {
            if (boost::multiprecision::abs(a.real() - b.real()) < std::numeric_limits<mpfr_50>::epsilon())
                return (a.imag() < b.imag());
            else
                return (a.real() < b.real());
        }
    );
    for (int i = 0; i < 5; ++i)
    {
        std::tuple<bool, std::string, int> root_real = getNumberAsString<50>(roots_50.first[i].real());
        std::tuple<bool, std::string, int> root_imag = getNumberAsString<50>(roots_50.first[i].imag());

        // Check all 50 digits
        if (i == 0)
        {
            BOOST_TEST(!std::get<0>(root_real));
            BOOST_TEST((std::get<1>(root_real).rfind("8.0115515308340360926152290004316066895772474042464", 0) == 0 || std::get<1>(root_real).rfind("8.0115515308340360926152290004316066895772474042463", 0) == 0));
            BOOST_TEST(std::get<2>(root_real) == -1);
            BOOST_TEST(boost::multiprecision::abs(roots.first[i].imag()) < std::numeric_limits<mpfr_50>::epsilon());
        }
        else if (i == 1)
        {
            BOOST_TEST(!std::get<0>(root_real));
            BOOST_TEST((std::get<1>(root_real).rfind("4.2613327731896388789584586349206480661427590651186", 0) == 0 || std::get<1>(root_real).rfind("4.2613327731896388789584586349206480661427590651185", 0) == 0));
            BOOST_TEST(std::get<2>(root_real) == -1);
            BOOST_TEST(!std::get<0>(root_imag));
            BOOST_TEST((std::get<1>(root_imag).rfind("5.7039301146186323645735618066518508239107126817095", 0) == 0 || std::get<1>(root_imag).rfind("5.7039301146186323645735618066518508239107126817094", 0) == 0));
            BOOST_TEST(std::get<2>(root_imag) == -1);
        }
        else if (i == 2)
        {
            BOOST_TEST(!std::get<0>(root_real));
            BOOST_TEST((std::get<1>(root_real).rfind("4.2613327731896388789584586349206480661427590651186", 0) == 0 || std::get<1>(root_real).rfind("4.2613327731896388789584586349206480661427590651185", 0) == 0));
            BOOST_TEST(std::get<2>(root_real) == -1);
            BOOST_TEST(std::get<0>(root_imag));
            BOOST_TEST((std::get<1>(root_imag).rfind("5.7039301146186323645735618066518508239107126817095", 0) == 0 || std::get<1>(root_imag).rfind("5.7039301146186323645735618066518508239107126817094", 0) == 0));
            BOOST_TEST(std::get<2>(root_imag) == -1);
        }
        else if (i == 3)
        {
            BOOST_TEST(std::get<0>(root_real));
            BOOST_TEST((std::get<1>(root_real).rfind("5.6075340705215505422873497308811322619952125544759", 0) == 0 || std::get<1>(root_real).rfind("5.6075340705215505422873497308811322619952125544758", 0) == 0));
            BOOST_TEST(std::get<2>(root_real) == -1);
            BOOST_TEST(!std::get<0>(root_imag));
            BOOST_TEST((std::get<1>(root_imag).rfind("4.9872376016046997943682112757504653134116713321119", 0) == 0 || std::get<1>(root_imag).rfind("4.9872376016046997943682112757504653134116713321118", 0) == 0));
            BOOST_TEST(std::get<2>(root_imag) == -1);
        }
        else if (i == 4)
        {
            BOOST_TEST(std::get<0>(root_real));
            BOOST_TEST((std::get<1>(root_real).rfind("5.6075340705215505422873497308811322619952125544759", 0) == 0 || std::get<1>(root_real).rfind("5.6075340705215505422873497308811322619952125544758", 0) == 0));
            BOOST_TEST(std::get<2>(root_real) == -1);
            BOOST_TEST(std::get<0>(root_imag));
            BOOST_TEST((std::get<1>(root_imag).rfind("4.9872376016046997943682112757504653134116713321119", 0) == 0 || std::get<1>(root_imag).rfind("4.9872376016046997943682112757504653134116713321118", 0) == 0));
            BOOST_TEST(std::get<2>(root_imag) == -1);
        }
    }
}
