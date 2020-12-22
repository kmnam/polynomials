#define BOOST_TEST_MODULE testPolynomial
#define BOOST_TEST_DYN_LINK
#include <iostream>
#include <string>
#include <sstream>
#include <limits>
#include <utility>
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
 *     12/22/2020
 */
using boost::multiprecision::number; 
using boost::multiprecision::mpfr_float_backend; 
using boost::multiprecision::mpc_complex_backend; 

template <unsigned N>
std::pair<std::string, int> getNumberAsString(number<mpfr_float_backend<N> > num)
{
    /*
     * Return a maximally precise string representation of the given number. 
     */
    std::stringstream ss;

    // Feed number into stringstream with maximum precision
    ss << std::setprecision(std::numeric_limits<number<mpfr_float_backend<N> > >::max_digits10)
       << std::scientific;
    ss << num;

    // Get output string and parse mantissa and exponent 
    std::string s = ss.str();
    std::string mantissa = s.substr(0, s.find("e"));
    s.erase(0, s.find("e") + 1);
    int exponent = std::stoi(s);

    return std::make_pair(mantissa, exponent);
}

template <unsigned N>
std::vector<number<mpfr_float_backend<N> > > logRandomCoefs(const unsigned deg,
                                                            const double min_exp,
                                                            const double max_exp,
                                                            boost::random::mt19937 rng)
{
    /*
     * Return a vector of randomly sampled coefficients for a polynomial 
     * of given degree, from a log-uniform distribution with the given
     * min/max exponents.
     *
     * The numbers are sampled first as doubles, then cast to the given 
     * precision (N).  
     */
    boost::random::uniform_real_distribution<> dist(min_exp, max_exp);
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
    std::pair<std::string, int> coef = getNumberAsString<20>(g.coefficients()[0]);
    BOOST_TEST(coef.first.size() > 20);
    BOOST_TEST(coef.first.rfind("5.0000000000000000000", 0) == 0);
    BOOST_TEST(coef.second == 0);

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
        BOOST_TEST(coef.first.size() > 20);
        if      (i == 0) BOOST_TEST(coef.first.rfind("1.0000000000000000000", 0) == 0);
        else if (i == 1) BOOST_TEST(coef.first.rfind("2.0000000000000000000", 0) == 0);
        else             BOOST_TEST(coef.first.rfind("3.0000000000000000000", 0) == 0);
        BOOST_TEST(coef.second == 0);
    }

    // Testing constructor with double scalar argument
    double value = 17;
    Polynomial<20> k(value);
    BOOST_TEST(k.degree() == 0);
    BOOST_TEST(k.coefficients().size() == 1);
    BOOST_TEST(k.coefficients()[0] == 17.0);

    // Testing precision of stored coefficient 
    coef = getNumberAsString<20>(k.coefficients()[0]);
    BOOST_TEST(coef.first.size() > 20);
    BOOST_TEST(coef.first.rfind("1.7000000000000000000", 0) == 0);
    BOOST_TEST(coef.second == 1);

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
        BOOST_TEST(coef.first.size() > 20);
        if      (i == 0) BOOST_TEST(coef.first.rfind("7.1000000000000000000", 0) == 0);
        else if (i == 1) BOOST_TEST(coef.first.rfind("4.2000000000000000000", 0) == 0);
        else             BOOST_TEST(coef.first.rfind("9.9000000000000000000", 0) == 0);
        BOOST_TEST(coef.second == 1);
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
    std::pair<std::string, int> zc0 = getNumberAsString<20>(z0.real());
    mpc_20 z1 = f.eval<20, 20>(z);
    std::pair<std::string, int> zc1 = getNumberAsString<20>(z1.real());
    BOOST_TEST(z0 == 86.0);
    BOOST_TEST(zc0.first.size() > 20);
    BOOST_TEST(zc0.first.rfind("8.6000000000000000000", 0) == 0);
    BOOST_TEST(zc0.second == 1);
    BOOST_TEST(z1 == 86.0);
    BOOST_TEST(zc1.first.size() > 20);
    BOOST_TEST(zc1.first.rfind("8.6000000000000000000", 0) == 0);
    BOOST_TEST(zc1.second == 1);

    // Evaluate with constant precision (20) and std::complex<double> input
    std::complex<double> zd(5.0, 0.0);
    mpc_20 zd0 = f.eval(zd);
    std::pair<std::string, int> zdc0 = getNumberAsString<20>(zd0.real());
    mpc_20 zd1 = f.eval<20>(zd);
    std::pair<std::string, int> zdc1 = getNumberAsString<20>(zd1.real());
    BOOST_TEST(zd0 == 86.0);
    BOOST_TEST(zdc0.first.size() > 20);
    BOOST_TEST(zdc0.first.rfind("8.6000000000000000000", 0) == 0);
    BOOST_TEST(zdc0.second == 1);
    BOOST_TEST(zd1 == 86.0);
    BOOST_TEST(zdc1.first.size() > 20);
    BOOST_TEST(zdc1.first.rfind("8.6000000000000000000", 0) == 0);
    BOOST_TEST(zdc1.second == 1);

    // Evaluate with constant precision (20) and real input and real output 
    mpfr_20 x = 5.0;
    mpfr_20 x0 = f.eval(x);
    std::pair<std::string, int> xc0 = getNumberAsString<20>(x0);
    mpfr_20 x1 = f.eval<20, 20>(x);
    std::pair<std::string, int> xc1 = getNumberAsString<20>(x1);
    BOOST_TEST(x0 == 86.0);
    BOOST_TEST(xc0.first.size() > 20);
    BOOST_TEST(xc0.first.rfind("8.6000000000000000000", 0) == 0);
    BOOST_TEST(xc0.second == 1);
    BOOST_TEST(x1 == 86.0);
    BOOST_TEST(xc1.first.size() > 20);
    BOOST_TEST(xc1.first.rfind("8.6000000000000000000", 0) == 0);
    BOOST_TEST(xc1.second == 1);

    // Evaluate with constant precision (20) and double input
    double xd = 5.0;
    mpfr_20 xd0 = f.eval(xd);
    std::pair<std::string, int> xdc0 = getNumberAsString<20>(xd0.real());
    mpfr_20 xd1 = f.eval<20>(xd);
    std::pair<std::string, int> xdc1 = getNumberAsString<20>(xd1.real());
    BOOST_TEST(xd0 == 86.0);
    BOOST_TEST(xdc0.first.size() > 20);
    BOOST_TEST(xdc0.first.rfind("8.6000000000000000000", 0) == 0);
    BOOST_TEST(xdc0.second == 1);
    BOOST_TEST(xd1 == 86.0);
    BOOST_TEST(xdc1.first.size() > 20);
    BOOST_TEST(xdc1.first.rfind("8.6000000000000000000", 0) == 0);
    BOOST_TEST(xdc1.second == 1);

    // Evaluate with input precision = 30 and output precision = 20
    mpc_30 a(5.0, 0.0);
    mpc_20 a0 = f.eval<30, 20>(a);
    std::pair<std::string, int> ac0 = getNumberAsString<20>(a0.real());
    mpfr_30 y = 5.0;
    mpfr_20 y0 = f.eval<30, 20>(y);
    std::pair<std::string, int> yc0 = getNumberAsString<20>(y0);
    BOOST_TEST(a0 == 86.0);
    BOOST_TEST(ac0.first.size() > 20);
    BOOST_TEST(ac0.first.rfind("8.6000000000000000000", 0) == 0);
    BOOST_TEST(ac0.second == 1);
    BOOST_TEST(y0 == 86.0);
    BOOST_TEST(yc0.first.size() > 20);
    BOOST_TEST(yc0.first.rfind("8.6000000000000000000", 0) == 0);
    BOOST_TEST(yc0.second == 1);

    // Evaluate with input precision = 30 and output precision = 40
    mpc_40 a2 = f.eval<30, 40>(a);
    std::pair<std::string, int> ac2 = getNumberAsString<40>(a2.real());
    mpfr_40 y2 = f.eval<30, 40>(y);
    std::pair<std::string, int> yc2 = getNumberAsString<40>(y2);
    BOOST_TEST(a2 == 86.0);
    BOOST_TEST(ac2.first.size() > 40);
    BOOST_TEST(ac2.first.rfind("8.600000000000000000000000000000000000000", 0) == 0);
    BOOST_TEST(ac2.second == 1);
    BOOST_TEST(y2 == 86.0);
    BOOST_TEST(yc2.first.size() > 40);
    BOOST_TEST(yc2.first.rfind("8.600000000000000000000000000000000000000", 0) == 0);
    BOOST_TEST(yc2.second == 1);
}

BOOST_AUTO_TEST_CASE(testOperators)
{
    /*
     * Test the addition, subtraction, multiplication, and division operators.  
     */
    typedef number<mpfr_float_backend<20> > mpfr_20;
    typedef number<mpfr_float_backend<30> > mpfr_30;
    typedef number<mpfr_float_backend<40> > mpfr_40;

    // Define three polynomials of varying precisions and degrees 
    std::vector<mpfr_20> p_coefs;           // 41.5 - 74.6*x + 19.6*x^2
    p_coefs.push_back(mpfr_20("41.5"));
    p_coefs.push_back(mpfr_20("-74.6"));
    p_coefs.push_back(mpfr_20("19.6"));
    Polynomial<20> p(p_coefs);
    for (int i = 0; i < 3; ++i)
    {
        std::pair<std::string, int> coef = getNumberAsString<20>(p.coefficients()[i]);
        if (i == 0)      BOOST_TEST((coef.first.rfind("4.1500000000000000000", 0) == 0  || coef.first.rfind("4.1499999999999999999", 0) == 0));
        else if (i == 1) BOOST_TEST((coef.first.rfind("-7.4600000000000000000", 0) == 0 || coef.first.rfind("-7.4599999999999999999", 0) == 0));
        else if (i == 2) BOOST_TEST((coef.first.rfind("1.9600000000000000000", 0) == 0  || coef.first.rfind("1.9599999999999999999", 0) == 0));
        BOOST_TEST(coef.second == 1); 
    }
    std::vector<mpfr_30> q_coefs;           // 55.9 + 58.5*x - 96.5*x^2 + 83.4*x^3
    q_coefs.push_back(mpfr_30("55.9"));
    q_coefs.push_back(mpfr_30("58.5"));
    q_coefs.push_back(mpfr_30("-96.5"));
    q_coefs.push_back(mpfr_30("83.4"));
    Polynomial<30> q(q_coefs);
    for (int i = 0; i < 4; ++i)
    {
        std::pair<std::string, int> coef = getNumberAsString<30>(q.coefficients()[i]);
        if (i == 0)      BOOST_TEST((coef.first.rfind("5.59000000000000000000000000000", 0) == 0  || coef.first.rfind("5.589999999999999999999999999999", 0) == 0));
        else if (i == 1) BOOST_TEST((coef.first.rfind("5.85000000000000000000000000000", 0) == 0  || coef.first.rfind("5.849999999999999999999999999999", 0) == 0));
        else if (i == 2) BOOST_TEST((coef.first.rfind("-9.65000000000000000000000000000", 0) == 0 || coef.first.rfind("-9.649999999999999999999999999999", 0) == 0));
        else if (i == 3) BOOST_TEST((coef.first.rfind("8.34000000000000000000000000000", 0) == 0  || coef.first.rfind("8.339999999999999999999999999999", 0) == 0));
        BOOST_TEST(coef.second == 1);
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
        std::pair<std::string, int> coef = getNumberAsString<40>(r.coefficients()[i]);
        if (i == 0)
            BOOST_TEST((coef.first.rfind("1.290000000000000000000000000000000000000", 0) == 0 || coef.first.rfind("1.289999999999999999999999999999999999999", 0) == 0));
        else if (i == 1)
            BOOST_TEST((coef.first.rfind("1.210000000000000000000000000000000000000", 0) == 0 || coef.first.rfind("1.209999999999999999999999999999999999999", 0) == 0));
        else if (i == 2)
            BOOST_TEST((coef.first.rfind("1.690000000000000000000000000000000000000", 0) == 0 || coef.first.rfind("1.689999999999999999999999999999999999999", 0) == 0));
        else if (i == 3)
            BOOST_TEST((coef.first.rfind("-5.720000000000000000000000000000000000000", 0) == 0 || coef.first.rfind("-5.719999999999999999999999999999999999999", 0) == 0));
        else if (i == 4)
            BOOST_TEST((coef.first.rfind("3.000000000000000000000000000000000000000", 0) == 0 || coef.first.rfind("2.999999999999999999999999999999999999999", 0) == 0));
        else if (i == 5)
            BOOST_TEST((coef.first.rfind("5.640000000000000000000000000000000000000", 0) == 0 || coef.first.rfind("5.639999999999999999999999999999999999999", 0) == 0));
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
        std::pair<std::string, int> coef = getNumberAsString<40>(s.coefficients()[i]);
        if (i == 0)
            BOOST_TEST((coef.first.rfind("5.250000000000000000000000000000000000000", 0) == 0  || coef.first.rfind("5.249999999999999999999999999999999999999", 0) == 0));
        else if (i == 1)
            BOOST_TEST((coef.first.rfind("5.880000000000000000000000000000000000000", 0) == 0  || coef.first.rfind("5.879999999999999999999999999999999999999", 0) == 0));
        else if (i == 2)
            BOOST_TEST((coef.first.rfind("-8.660000000000000000000000000000000000000", 0) == 0 || coef.first.rfind("-8.659999999999999999999999999999999999999", 0) == 0));
        else if (i == 3)
            BOOST_TEST((coef.first.rfind("7.470000000000000000000000000000000000000", 0) == 0  || coef.first.rfind("7.469999999999999999999999999999999999999", 0) == 0));
        else if (i == 4)
            BOOST_TEST((coef.first.rfind("-2.250000000000000000000000000000000000000", 0) == 0 || coef.first.rfind("-2.249999999999999999999999999999999999999", 0) == 0));
    }

    // Add p and q to get a result of precision 40
    Polynomial<40> p_plus_q = p.operator+<30, 40>(q);
    std::vector<mpfr_40> sum_coefs = p_plus_q.coefficients();
    BOOST_TEST(p_plus_q.degree() == 3); 
    BOOST_TEST(sum_coefs.size() == 4);
    for (int i = 0; i < 4; ++i)
    {
        std::pair<std::string, int> coef = getNumberAsString<40>(sum_coefs[i]);
        if (i == 0)      
            BOOST_TEST((coef.first.rfind("9.740000000000000000000000000000000000000", 0) == 0  || coef.first.rfind("9.739999999999999999999999999999999999999", 0) == 0));
        else if (i == 1)
            BOOST_TEST((coef.first.rfind("-1.610000000000000000000000000000000000000", 0) == 0 || coef.first.rfind("-1.609999999999999999999999999999999999999", 0) == 0));
        else if (i == 2)
            BOOST_TEST((coef.first.rfind("-7.690000000000000000000000000000000000000", 0) == 0 || coef.first.rfind("-7.689999999999999999999999999999999999999", 0) == 0));
        else if (i == 3)
            BOOST_TEST((coef.first.rfind("8.3400000000000000000000000000000000000000", 0) == 0 || coef.first.rfind("8.339999999999999999999999999999999999999", 0) == 0));
        BOOST_TEST(coef.second == 1);
    }

    // Add r and s (both precision 40, output precision 40)
    Polynomial<40> r_plus_s = r + s;
    sum_coefs = r_plus_s.coefficients();
    BOOST_TEST(r_plus_s.degree() == 5);
    BOOST_TEST(sum_coefs.size() == 6);
    for (int i = 0; i < 6; ++i)
    {
        std::pair<std::string, int> coef = getNumberAsString<40>(sum_coefs[i]);
        if (i == 0)
        {      
            BOOST_TEST((coef.first.rfind("1.815000000000000000000000000000000000000", 0) == 0  || coef.first.rfind("1.814999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(coef.second == 1);
        }
        else if (i == 1)
        {
            BOOST_TEST((coef.first.rfind("1.798000000000000000000000000000000000000", 0) == 0  || coef.first.rfind("1.797999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(coef.second == 1);
        }
        else if (i == 2)
        {
            BOOST_TEST((coef.first.rfind("-8.491000000000000000000000000000000000000", 0) == 0 || coef.first.rfind("-8.490999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(coef.second == 0);
        }
        else if (i == 3)
        {
            BOOST_TEST((coef.first.rfind("6.898000000000000000000000000000000000000", 0) == 0  || coef.first.rfind("6.897999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(coef.second == 1);
        }
    } 

    // Subtract q from p to get a result of precision 40
    Polynomial<40> p_minus_q = p.operator-<30, 40>(q);
    std::vector<mpfr_40> diff_coefs = p_minus_q.coefficients();
    BOOST_TEST(p_minus_q.degree() == 3); 
    BOOST_TEST(diff_coefs.size() == 4);
    for (int i = 0; i < 4; ++i)
    {
        std::pair<std::string, int> coef = getNumberAsString<40>(diff_coefs[i]);
        if (i == 0)
        {      
            BOOST_TEST((coef.first.rfind("-1.440000000000000000000000000000000000000", 0) == 0 || coef.first.rfind("-1.439999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(coef.second == 1);
        }
        else if (i == 1)
        {
            BOOST_TEST((coef.first.rfind("-1.331000000000000000000000000000000000000", 0) == 0 || coef.first.rfind("-1.330999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(coef.second == 2);
        }
        else if (i == 2)
        {
            BOOST_TEST((coef.first.rfind("1.161000000000000000000000000000000000000", 0) == 0  || coef.first.rfind("1.160999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(coef.second == 2);
        }
        else if (i == 3)
        {
            BOOST_TEST((coef.first.rfind("-8.340000000000000000000000000000000000000", 0) == 0 || coef.first.rfind("-8.339999999999999999999999999999999999999", 0) == 0));
            BOOST_TEST(coef.second == 1);
        }
    }
}
