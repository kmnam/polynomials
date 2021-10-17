#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <utility>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpc.hpp>
#include "../include/polynomial/polynomial.hpp"

/*
 * A simple wrapper script that uses the Polynomial class to calculate the
 * zeros of univariate polynomials.
 *
 * Coefficients are stored to 100-digit precision with Boost/MPFR.  
 *
 * Authors:
 *     Kee-Myoung Nam, Department of Systems Biology, Harvard Medical School
 * Last updated:
 *     10/17/2021
 */
using boost::multiprecision::number;
using boost::multiprecision::mpfr_float_backend;

template <unsigned N>
std::vector<Polynomial<N> > parsePolynomials(std::string filename)
{
    /*
     * Parse a polynomial from a comma-separated input file.
     */
    typedef number<mpfr_float_backend<N> > RTN;

    std::vector<Polynomial<N> > polynomials; 
    std::ifstream infile(filename);
    std::string line;
    while (std::getline(infile, line, '\n'))
    {
        // Each line contains a list of coefficients
        std::vector<RTN> coefs; 

        std::istringstream iss(line);
        std::string coef;
        while (std::getline(iss, coef, ','))
        {
            RTN x(coef);
            coefs.push_back(x); 
        }

        // Instantiate polynomial 
        polynomials.emplace_back(Polynomial<N>(coefs)); 
    }

    return polynomials; 
} 

int main(int argc, char* argv[])
{
    /*
     * Parse the given input file of polynomial coefficients and solve each 
     * polynomial. 
     */
    // Set up seeded Boost random number generator 
    boost::random::mt19937 rng(1234567890);
    boost::random::uniform_real_distribution<double> dist;

    const unsigned max_iter = 1000;                       // Max of 1000 iterations 
    const number<mpfr_float_backend<100> > tol = 1e-80;   // Tolerance = 1e-80
    const bool verbose = true;                            // Print all intermediate output 

    // Parse polynomials from input file, storing coefficients up to 100 digits 
    std::string filename(argv[1]);
    std::vector<Polynomial<100> > polynomials = parsePolynomials<100>(filename);

    // Solve each polynomial ...
    for (auto&& poly : polynomials)
    {
        // ... by applying Aberth's method
        std::pair<std::vector<number<mpc_complex_backend<100> > >, bool> roots
            = poly.roots(max_iter, tol, tol, rng, dist, verbose);

        if (roots.second)
        {
            for (auto&& root : roots.first)
                std::cout << root.real() << " " << root.imag() << std::endl;
        }
        else 
        {
            std::cout << "Did not converge" << std::endl;
        }
    }
}
