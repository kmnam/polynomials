#define BOOST_TEST_MODULE testComplexNumber
#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>
#include "../../include/complex/complex.hpp"
#include "../../include/complex/eigen.hpp"

/*
 * Test module for the Duals::DualNumber class.
 *
 * Authors:
 *     Kee-Myoung Nam, Department of Systems Biology, Harvard Medical School
 * Last updated:
 *     11/14/2019
 */
typedef ComplexNumber<double> cdouble;

BOOST_AUTO_TEST_CASE(testInitialize)
{
    cdouble z;
    BOOST_TEST(z.real() == 0.0);
    BOOST_TEST(z.imag() == 0.0);
}

BOOST_AUTO_TEST_CASE(testInitializeFromDouble)
{
    cdouble z(10.0);
    BOOST_TEST(z.real() == 10.0);
    BOOST_TEST(z.imag() == 0.0);

    cdouble w = 10.0;
    BOOST_TEST(z.real() == 10.0);
    BOOST_TEST(z.imag() == 0.0);
}

BOOST_AUTO_TEST_CASE(testSum)
{
    cdouble z(3.14159, 2.71828);
    cdouble w(1.41421, 2.23606);
    cdouble sum = z + w;
    BOOST_TEST(sum.real() == 3.14159 + 1.41421);
    BOOST_TEST(sum.imag() == 2.71828 + 2.23606);

    z += w;
    BOOST_TEST(z.real() == 3.14159 + 1.41421);
    BOOST_TEST(z.imag() == 2.71828 + 2.23606);

    double x = 17.76;
    BOOST_TEST((z + x).real() == 3.14159 + 1.41421 + 17.76);
    BOOST_TEST((x + z).real() == 3.14159 + 1.41421 + 17.76);
}

