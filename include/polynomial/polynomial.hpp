#ifndef UNIVAR_POLYNOMIAL_HPP
#define UNIVAR_POLYNOMIAL_HPP

#include <utility>
#include <ostream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <complex>
#include <limits>
#include <boost/random.hpp>
#include <Eigen/Dense>
#include <boost/multiprecision/mpc.hpp>

/*
 * A Polynomial class template with arbitrary real coefficients.
 *
 * The coefficients used in this class should be those that support usage
 * with std::complex.  
 *
 * Authors:
 *     Kee-Myoung Nam, Department of Systems Biology, Harvard Medical School
 * Last updated:
 *     2/18/2020
 */
using namespace Eigen;
using boost::multiprecision::number;
using boost::multiprecision::mpc_complex_backend;
typedef number<mpc_complex_backend<30> > mpc_30;
typedef number<mpc_complex_backend<60> > mpc_60;
typedef number<mpc_complex_backend<100> > mpc_100;
typedef number<mpc_complex_backend<200> > mpc_200;

// ------------------------------------------------------------------ //
//                      VARIOUS HELPER FUNCTIONS                      //
// ------------------------------------------------------------------ //
enum SolveMethod 
{
    Weierstrass,
    Aberth
};

template <typename T>
Matrix<T, Dynamic, 1> appendZeros(const Ref<const Matrix<T, Dynamic, 1> >& v, unsigned j)
{
    /*
     * Pad the given vector, v, with j zeros. 
     */
    Matrix<T, Dynamic, 1> w(v);
    w.conservativeResize(w.size() + j);
    for (unsigned i = v.size(); i < w.size(); ++i) w(i) = 0;
    return w;
}

template <typename T>
Matrix<T, Dynamic, 1> removeTrailingZeros(const Ref<const Matrix<T, Dynamic, 1> >& v)
{
    /*
     * Remove trailing zeros from v.
     */
    Matrix<T, Dynamic, 1> w(v);
    int size = v.size();
    for (int i = v.size() - 1; i >= 0; --i)
    {
        if (v(i) == 0)
        {
            size--;
            w = w.head(size).eval();
        }
        else break;
    }
    return w;
}

template <typename T>
std::vector<std::complex<T> > quadratic(const T b, const T c)
{
    /*
     * Solve the given monic quadratic. 
     */
    std::vector<std::complex<T> > roots; 
    std::complex<T> disc = std::sqrt(std::complex<T>(b * b - 4 * c, 0.0));
    roots.push_back((-b + disc) / 2.0);
    roots.push_back((-b - disc) / 2.0);
    return roots;
}

template <typename T>
std::complex<T> cuberoot(std::complex<T> z)
{
    /*
     * Return the principal cube root of z.
     */
    if (z.real() < 0)
        return -std::pow(-z, 1.0 / 3.0);
    else
        return std::pow(z, 1.0 / 3.0);
}

template <typename T>
std::vector<std::complex<T> > cubic(const T b, const T c, const T d)
{
    /*
     * Solve the given monic cubic. 
     */
    std::vector<std::complex<T> > roots;
    T q = (3 * c - b * b) / 9;
    T r = (9 * b * c - 27 * d - 2 * std::pow(b, 3)) / 54;
    T p = std::pow(q, 3) + std::pow(r, 2);
    std::complex<T> s = cuberoot(r + std::sqrt(std::complex<T>(p, 0.0)));
    std::complex<T> t = cuberoot(r - std::sqrt(std::complex<T>(p, 0.0)));
    roots.push_back(-b / 3.0 + s + t);
    roots.push_back(-b / 3.0 - (s + t) / 2.0 + std::complex<T>(0.0, 1.0) * std::sqrt(3) * (s - t) / 2.0);
    roots.push_back(std::conj(roots[1]));
    return roots;
}

template <typename CT>
CT horner(const std::vector<CT>& coefs, CT x)
{
    /*
     * Perform Horner's method, given the vector of coefficients corresponding
     * to the polynomial and a value at which the polynomial is to be evaluated. 
     */ 
    CT y = coefs[coefs.size() - 1];
    for (int i = coefs.size() - 2; i >= 0; --i)
    {
        y *= x;
        y += coefs[i];
    }
    return y;
}

template <typename CT>
std::pair<std::vector<CT>, std::vector<double> > newton(const std::vector<CT>& coefs,
                                                        const std::vector<CT>& dcoefs,
                                                        const std::vector<CT>& roots)
{
    /*
     * Perform one iteration of Newton's method, given the vector of 
     * coefficients corresponding to the polynomial, its derivative, and 
     * a vector of current root approximations. 
     */
    // Perform the Newton correction for each root
    std::vector<CT> new_roots(roots);
    for (unsigned j = 0; j < roots.size(); ++j)
    {
        new_roots[j] -= (horner<CT>(coefs, roots[j]) / horner<CT>(dcoefs, roots[j]));
    }

    // Compute the absolute difference between each old and new root
    std::vector<double> delta;
    for (unsigned j = 0; j < roots.size(); ++j)
    {
        double dr = static_cast<double>(roots[j].real() - new_roots[j].real());
        double di = static_cast<double>(roots[j].imag() - new_roots[j].imag());
        delta.push_back(std::sqrt(std::pow(dr, 2.0) + std::pow(di, 2.0)));
    }

    return std::make_pair(new_roots, delta);
}

template <typename CT>
std::pair<std::vector<CT>, std::vector<double> > weierstrass(const std::vector<CT>& coefs,
                                                             const std::vector<CT>& roots)
{
    /*
     * Perform one iteration of Weierstrass' method, given the vector 
     * of coefficients corresponding to the polynomial and a vector of 
     * current root approximations.  
     */
    // Perform the Weierstrass correction for each root
    std::vector<CT> new_roots(roots);
    for (unsigned j = 0; j < roots.size(); ++j)
    {
        CT denom = 1.0;
        for (unsigned k = 0; k < roots.size(); ++k)
        {
            if (j != k) denom *= (roots[j] - new_roots[k]);
        }
        CT correction = -horner<CT>(coefs, roots[j]) / denom;
        new_roots[j] += correction;
    }

    // Compute the absolute difference between each old and new root
    std::vector<double> delta;
    for (unsigned j = 0; j < roots.size(); ++j)
    {
        double dr = static_cast<double>(roots[j].real() - new_roots[j].real());
        double di = static_cast<double>(roots[j].imag() - new_roots[j].imag());
        delta.push_back(std::sqrt(std::pow(dr, 2.0) + std::pow(di, 2.0)));
    }

    return std::make_pair(new_roots, delta);
}

template <typename CT>
std::pair<std::vector<CT>, std::vector<double> > aberth(const std::vector<CT>& coefs,
                                                        const std::vector<CT>& dcoefs,
                                                        const std::vector<CT>& roots)
{
    /*
     * Perform one iteration of the Aberth-Ehrlich method, given the vector 
     * of coefficients corresponding to the polynomial and a vector of 
     * current root approximations.  
     */
    // Perform the Aberth-Ehrlich correction for each root
    std::vector<CT> new_roots(roots);
    for (unsigned j = 0; j < roots.size(); ++j)
    {
        CT sum = 0.0;
        for (unsigned k = 0; k < j; ++k)
        {
            sum += (1.0 / (roots[j] - new_roots[k]));
        }
        for (unsigned k = j + 1; k < roots.size(); ++k)
        {
            sum += (1.0 / (roots[j] - roots[k])); 
        }
        CT val = horner<CT>(coefs, roots[j]) / horner<CT>(dcoefs, roots[j]);
        CT denom = 1.0 - (val * sum);
        CT correction = val / denom;
        new_roots[j] -= correction;
    }

    // Compute the absolute difference between each old and new root
    std::vector<double> delta;
    for (unsigned j = 0; j < roots.size(); ++j)
    {
        double dr = static_cast<double>(roots[j].real() - new_roots[j].real());
        double di = static_cast<double>(roots[j].imag() - new_roots[j].imag());
        delta.push_back(std::sqrt(std::pow(dr, 2.0) + std::pow(di, 2.0)));
    }
    
    return std::make_pair(new_roots, delta);
}

bool ccw(std::pair<double, double> x, std::pair<double, double> y, std::pair<double, double> z)
{
    /*
     * Determine if x, y, and z make a counterclockwise turn. 
     */
    double cross = (y.first - x.first) * (z.second - x.second) - (y.second - x.second) * (z.first - x.first);
    return (cross > 0);
}

std::vector<std::complex<double> > biniInitialize(const Ref<const VectorXd>& coefs,
                                                  boost::random::mt19937& rng,
                                                  boost::random::uniform_real_distribution<double>& dist)
{
    /*
     * Use Bini's initialization to yield a set of initial complex roots 
     * for the given vector of polynomial coefficients. 
     */
    using std::abs;
    using std::log2;
    using std::pow;
    using std::sin;
    using std::cos;
    const double two_pi = 2 * std::acos(-1.0); 

    // Find the upper convex hull of the vertices (i, log|f_i|) with 
    // the Andrew scan
    std::vector<std::pair<double, double> > points;
    for (unsigned i = 0; i < coefs.size(); ++i)
    {
        if (abs(coefs(i)) == 0)
            points.push_back(std::make_pair(i, -std::numeric_limits<double>::infinity()));
        else
            points.push_back(std::make_pair(i, log2(abs(coefs(i)))));
    }
    std::vector<std::pair<double, double> > upper_hull;
    for (int i = points.size() - 1; i >= 0; --i)
    {
        while (upper_hull.size() >= 2 && !ccw(upper_hull[upper_hull.size()-2], upper_hull[upper_hull.size()-1], points[i]))
            upper_hull.pop_back();
        upper_hull.push_back(points[i]);
    }

    // Run through the convex hull vertices, grab their x-coordinates, and sort
    std::vector<double> hull_x;
    for (auto&& p : upper_hull)
    {
        if (!std::isinf(p.second))
            hull_x.push_back(static_cast<int>(p.first));
    }
    std::sort(hull_x.begin(), hull_x.end());

    // Compute u_{k_i} for each vertex in the convex hull
    std::vector<double> u;
    for (unsigned i = 0; i < hull_x.size() - 1; ++i)
    {
        int k = hull_x[i];
        int l = hull_x[i+1];
        u.push_back(pow(abs(coefs(k) / coefs(l)), 1.0 / (l - k)));
    }

    // Compute initial root approximations
    // TODO Check that these approximations are correct (how?) on small examples
    int n = coefs.size() - 1;
    std::vector<std::complex<double> > inits; 
    for (unsigned i = 0; i < hull_x.size() - 1; ++i)
    {
        int k = hull_x[i];
        int l = hull_x[i+1];
        for (unsigned j = 0; j < l - k; ++j)
        {
            double theta = (two_pi / (l - k)) * j + (two_pi * i / n) + dist(rng);
            std::complex<double> z(cos(theta), sin(theta));
            inits.push_back(u[i] * z);
        }
    }
    return inits;
}

template <typename T>
class Polynomial
{
    private:
        unsigned deg;                    // Degree of the polynomial
        unsigned prec;                   // Maximum number of digits that can
                                         // be represented by type T 
        Matrix<T, Dynamic, 1> coefs;     // Coefficients of the polynomial
                                         // stored in ascending order of power

        Matrix<std::complex<T>, Dynamic, 1> rootsWeierstrass(unsigned max_iter,
                                                             double atol,
                                                             double rtol,
                                                             std::vector<std::complex<double> > inits = {})
        {
            /*
             * Run Weierstrass' method on the given polynomial, returning the 
             * roots as type std::complex<T>.
             */
            using boost::multiprecision::pow;
            using boost::multiprecision::abs;

            bool converged = false;
            bool quadratic = false;
            bool within_tol = false;
            unsigned iter = 0;

            // Initialize the roots to (0.4 + 0.9 i)^p, for p = 0, ..., d - 1
            // if initializations were not specified 
            if (inits.size() == 0)
            {
                for (unsigned i = 0; i < this->deg; ++i)
                    inits.push_back(std::pow(std::complex<double>(0.4, 0.9), i));
            }
            // Throw exception if improper number of initializations were specified
            else if (inits.size() != this->deg)
            {
                throw std::runtime_error("Invalid number of initial roots specified");
            }

            // Start with type std::complex<double>
            std::vector<std::complex<double> > roots, coefs;
            for (unsigned i = 0; i < this->deg; ++i)
            {
                roots.push_back(inits[i]);
            }
            for (unsigned i = 0; i < this->coefs.size(); ++i)
            {
                std::complex<double> z(this->coefs(i), 0.0);
                coefs.push_back(z);
            }
            // Run the algorithm until all roots exhibit quadratic convergence
            // or the desired number of iterations is reached
            std::vector<double> delta, new_delta;
            for (unsigned i = 0; i < this->deg; ++i) delta.push_back(0.0);
            while (iter < max_iter && !converged)
            {
                auto result = weierstrass<std::complex<double> >(coefs, roots);
                std::vector<std::complex<double> > new_roots = result.first;
                new_delta = result.second;
                iter++;
                // Check that the change is less than the square of the 
                // previous change for each root
                quadratic = true;
                for (unsigned j = 0; j < this->deg; ++j)
                {
                    if (new_delta[j] > 0.01 * (delta[j] * delta[j]))
                    {
                        quadratic = false;
                        break;
                    }
                }
                // Check that the change is within the given tolerances
                within_tol = true;
                for (unsigned j = 0; j < this->deg; ++j)
                {
                    if (new_delta[j] > atol + rtol * std::abs(roots[j]))
                    {
                        within_tol = false;
                        break;
                    }
                }
                if (quadratic || within_tol) converged = true;
                roots = new_roots;
                delta = new_delta;
            }

            // If not converged, try again with a type with greater precision
            unsigned prec = 17;
            while (!converged)
            {
                if (prec < 30)
                {
                    // Re-do the same calculations with mpc_30
                    prec = 30; iter = 0;
                    std::vector<mpc_30> roots2, coefs2;
                    for (unsigned i = 0; i < this->deg; ++i)
                    {
                        roots2.push_back(mpc_30(roots[i].real(), roots[i].imag()));
                        delta[i] = 0.0;
                    }
                    for (unsigned i = 0; i < this->coefs.size(); ++i)
                    {
                        std::stringstream ss;
                        ss << std::setprecision(std::numeric_limits<T>::max_digits10) << this->coefs(i);
                        mpc_30 z(ss.str(), 0.0);
                        coefs2.push_back(z);
                    }
                    while (iter < max_iter && !converged)
                    {
                        auto result = weierstrass<mpc_30>(coefs2, roots2);
                        std::vector<mpc_30> new_roots = result.first;
                        new_delta = result.second;
                        iter++;
                        quadratic = true;
                        for (unsigned j = 0; j < this->deg; ++j)
                        {
                            if (new_delta[j] > 0.01 * (delta[j] * delta[j]))
                            {
                                quadratic = false;
                                break;
                            }
                        }
                        within_tol = true;
                        for (unsigned j = 0; j < this->deg; ++j)
                        {
                            if (new_delta[j] > atol + rtol * abs(roots2[j]))
                            {
                                within_tol = false;
                                break;
                            }
                        }
                        if (quadratic || within_tol) converged = true;
                        roots2 = new_roots;
                        delta = new_delta;
                    }
                    for (unsigned i = 0; i < this->deg; ++i)
                        roots[i] = static_cast<std::complex<double> >(roots2[i]);
                }
                else if (prec < 60)
                {
                    // Re-do the same calculations with mpc_60
                    prec = 60; iter = 0;
                    std::vector<mpc_60> roots2, coefs2;
                    for (unsigned i = 0; i < this->deg; ++i)
                    {
                        roots2.push_back(mpc_60(roots[i].real(), roots[i].imag()));
                        delta[i] = 0.0;
                    }
                    for (unsigned i = 0; i < this->coefs.size(); ++i)
                    {
                        std::stringstream ss;
                        ss << std::setprecision(std::numeric_limits<T>::max_digits10) << this->coefs(i);
                        mpc_60 z(ss.str(), 0.0);
                        coefs2.push_back(z);
                    }
                    while (iter < max_iter && !converged)
                    {
                        auto result = weierstrass<mpc_60>(coefs2, roots2);
                        std::vector<mpc_60> new_roots = result.first;
                        new_delta = result.second;
                        iter++;
                        quadratic = true;
                        for (unsigned j = 0; j < this->deg; ++j)
                        {
                            if (new_delta[j] > 0.01 * (delta[j] * delta[j]))
                            {
                                quadratic = false;
                                break;
                            }
                        }
                        within_tol = true;
                        for (unsigned j = 0; j < this->deg; ++j)
                        {
                            if (new_delta[j] > atol + rtol * abs(roots2[j]))
                            {
                                within_tol = false;
                                break;
                            }
                        }
                        if (quadratic || within_tol) converged = true;
                        roots2 = new_roots;
                        delta = new_delta;
                    }
                    for (unsigned i = 0; i < this->deg; ++i)
                        roots[i] = static_cast<std::complex<double> >(roots2[i]);
                }
                else if (prec < 100)
                {
                    // Re-do the same calculations with mpc_100
                    prec = 100; iter = 0;
                    std::vector<mpc_100> roots2, coefs2;
                    for (unsigned i = 0; i < this->deg; ++i)
                    {
                        roots2.push_back(mpc_100(roots[i].real(), roots[i].imag()));
                        delta[i] = 0.0;
                    }
                    for (unsigned i = 0; i < this->coefs.size(); ++i)
                    {
                        std::stringstream ss;
                        ss << std::setprecision(std::numeric_limits<T>::max_digits10) << this->coefs(i);
                        mpc_100 z(ss.str(), 0.0);
                        coefs2.push_back(z);
                    }
                    while (iter < max_iter && !converged)
                    {
                        auto result = weierstrass<mpc_100>(coefs2, roots2);
                        std::vector<mpc_100> new_roots = result.first;
                        new_delta = result.second;
                        iter++;
                        quadratic = true;
                        for (unsigned j = 0; j < this->deg; ++j)
                        {
                            if (new_delta[j] > 0.01 * (delta[j] * delta[j]))
                            {
                                quadratic = false;
                                break;
                            }
                        }
                        within_tol = true;
                        for (unsigned j = 0; j < this->deg; ++j)
                        {
                            if (new_delta[j] > atol + rtol * abs(roots2[j]))
                            {
                                within_tol = false;
                                break;
                            }
                        }
                        if (quadratic || within_tol) converged = true;
                        roots2 = new_roots;
                        delta = new_delta;
                    }
                    for (unsigned i = 0; i < this->deg; ++i)
                        roots[i] = static_cast<std::complex<double> >(roots2[i]);
                }
                else if (prec < 200)
                {
                    // Re-do the same calculations with mpc_200
                    prec = 200; iter = 0;
                    std::vector<mpc_200> roots2, coefs2;
                    for (unsigned i = 0; i < this->deg; ++i)
                    {
                        roots2.push_back(mpc_200(roots[i].real(), roots[i].imag()));
                        delta[i] = 0.0;
                    }
                    for (unsigned i = 0; i < this->coefs.size(); ++i)
                    {
                        std::stringstream ss;
                        ss << std::setprecision(std::numeric_limits<T>::max_digits10) << this->coefs(i);
                        mpc_200 z(ss.str(), 0.0);
                        coefs2.push_back(z);
                    }
                    while (iter < max_iter && !converged)
                    {
                        auto result = weierstrass<mpc_200>(coefs2, roots2);
                        std::vector<mpc_200> new_roots = result.first;
                        new_delta = result.second;
                        iter++;
                        quadratic = true;
                        for (unsigned j = 0; j < this->deg; ++j)
                        {
                            if (new_delta[j] > 0.01 * (delta[j] * delta[j]))
                            {
                                quadratic = false;
                                break;
                            }
                        }
                        within_tol = true;
                        for (unsigned j = 0; j < this->deg; ++j)
                        {
                            if (new_delta[j] > atol + rtol * abs(roots2[j]))
                            {
                                within_tol = false;
                                break;
                            }
                        }
                        if (quadratic || within_tol) converged = true;
                        roots2 = new_roots;
                        delta = new_delta;
                    }
                    for (unsigned i = 0; i < this->deg; ++i)
                        roots[i] = static_cast<std::complex<double> >(roots2[i]);
                }
                else
                {
                    // Give up at this point (don't throw an exception here, 
                    // as polynomials with singular roots will be encountered)
                    converged = true;
                }
            }

            Matrix<std::complex<T>, Dynamic, 1> final_roots(this->deg);
            for (unsigned i = 0; i < this->deg; ++i)
            {
                std::complex<T> z(static_cast<T>(roots[i].real()), static_cast<T>(roots[i].imag()));
                final_roots(i) = z;
            }
            return final_roots;
        }

        Matrix<std::complex<T>, Dynamic, 1> rootsAberth(unsigned max_iter,
                                                        double atol,
                                                        double rtol,
                                                        unsigned sharpen_iter = 20,
                                                        std::vector<std::complex<double> > inits = {})
        {
            /*
             * Run the Aberth-Ehrlich method on the given polynomial, returning the 
             * roots as type std::complex<T>.
             */
            using boost::multiprecision::pow;
            using boost::multiprecision::abs;

            bool converged = false;
            bool quadratic = false;
            bool within_tol = false;
            unsigned iter = 0;

            // Initialize the roots to (0.4 + 0.9 i)^p, for p = 0, ..., d - 1
            // if initializations were not specified 
            if (inits.size() == 0)
            {
                for (unsigned i = 0; i < this->deg; ++i)
                    inits.push_back(std::pow(std::complex<double>(0.4, 0.9), i));
            }
            // Throw exception if improper number of initializations were specified
            else if (inits.size() != this->deg)
            {
                throw std::runtime_error("Invalid number of initial roots specified");
            }

            // Start with type std::complex<double>
            VectorXd dcoefvec = this->diff().coefficients();
            std::vector<std::complex<double> > roots, coefs, dcoefs;
            for (unsigned i = 0; i < this->deg; ++i)
            {
                roots.push_back(inits[i]);
            }
            for (unsigned i = 0; i < this->coefs.size(); ++i)
            {
                std::complex<double> z(this->coefs(i), 0.0);
                coefs.push_back(z);
            }
            for (unsigned i = 0; i < dcoefvec.size(); ++i)
            {
                std::complex<double> z(dcoefvec(i), 0.0);
                dcoefs.push_back(z);
            }
            // Run the algorithm until all roots exhibit quadratic convergence
            // or the desired number of iterations is reached
            std::vector<double> delta, new_delta;
            for (unsigned i = 0; i < this->deg; ++i) delta.push_back(0.0);
            while (iter < max_iter && !converged)
            {
                auto result = aberth<std::complex<double> >(coefs, dcoefs, roots);
                std::vector<std::complex<double> > new_roots = result.first;
                new_delta = result.second;
                iter++;
                // Check that the change is less than the square of the 
                // previous change for each root
                for (unsigned j = 0; j < this->deg; ++j)
                {
                    quadratic = true;
                    if (new_delta[j] > 0.01 * (delta[j] * delta[j]))
                    {
                        quadratic = false;
                        break;
                    }
                }
                // Check that the change is within the given tolerances
                within_tol = true;
                for (unsigned j = 0; j < this->deg; ++j)
                {
                    if (new_delta[j] > atol + rtol * abs(roots[j]))
                    {
                        within_tol = false;
                        break;
                    }
                }
                if (quadratic || within_tol) converged = true;
                roots = new_roots;
                delta = new_delta;
            }
            // Sharpen each root for the given number of iterations 
            for (unsigned i = 0; i < sharpen_iter; ++i)
            {
                auto result = newton<std::complex<double> >(coefs, dcoefs, roots);
                roots = result.first;
            }

            // If not converged, try again with a type with greater precision
            unsigned prec = 17;
            while (!converged)
            {
                if (prec < 30)
                {
                    // Re-do the same calculations with mpc_30
                    prec = 30; iter = 0;
                    std::vector<mpc_30> roots2, coefs2, dcoefs2;
                    for (unsigned i = 0; i < this->deg; ++i)
                    {
                        roots2.push_back(mpc_30(roots[i].real(), roots[i].imag()));
                        delta[i] = 0.0;
                    }
                    for (unsigned i = 0; i < this->coefs.size(); ++i)
                    {
                        std::stringstream ss;
                        ss << std::setprecision(std::numeric_limits<T>::max_digits10) << this->coefs(i);
                        mpc_30 z(ss.str(), "0.0");
                        coefs2.push_back(z);
                    }
                    for (unsigned i = 0; i < dcoefvec.size(); ++i)
                    {
                        std::stringstream ss;
                        ss << std::setprecision(std::numeric_limits<T>::max_digits10) << dcoefvec(i);
                        mpc_30 z(ss.str(), "0.0");
                        dcoefs2.push_back(z);
                    }
                    while (iter < max_iter && !converged)
                    {
                        auto result = aberth<mpc_30>(coefs2, dcoefs2, roots2);
                        std::vector<mpc_30> new_roots = result.first;
                        new_delta = result.second;
                        iter++;
                        for (unsigned j = 0; j < this->deg; ++j)
                        {
                            quadratic = true;
                            if (new_delta[j] > 0.01 * (delta[j] * delta[j]))
                            {
                                quadratic = false;
                                break;
                            }
                        }
                        within_tol = true;
                        for (unsigned j = 0; j < this->deg; ++j)
                        {
                            if (new_delta[j] > atol + rtol * abs(roots2[j]))
                            {
                                within_tol = false;
                                break;
                            }
                        }
                        if (quadratic || within_tol) converged = true;
                        roots2 = new_roots;
                        delta = new_delta;
                    }
                    for (unsigned i = 0; i < sharpen_iter; ++i)
                    {
                        auto result = newton<mpc_30>(coefs2, dcoefs2, roots2);
                        roots2 = result.first;
                    }
                    for (unsigned i = 0; i < this->deg; ++i)
                        roots[i] = static_cast<std::complex<double> >(roots2[i]);
                }
                else if (prec < 60)
                {
                    // Re-do the same calculations with mpc_60
                    prec = 60; iter = 0;
                    std::vector<mpc_60> roots2, coefs2, dcoefs2;
                    for (unsigned i = 0; i < this->deg; ++i)
                    {
                        roots2.push_back(mpc_60(roots[i].real(), roots[i].imag()));
                        delta[i] = 0.0;
                    }
                    for (unsigned i = 0; i < this->coefs.size(); ++i)
                    {
                        std::stringstream ss;
                        ss << std::setprecision(std::numeric_limits<T>::max_digits10) << this->coefs(i);
                        mpc_60 z(ss.str(), 0.0);
                        coefs2.push_back(z);
                    }
                    for (unsigned i = 0; i < dcoefvec.size(); ++i)
                    {
                        std::stringstream ss;
                        ss << std::setprecision(std::numeric_limits<T>::max_digits10) << dcoefvec(i);
                        mpc_60 z(ss.str(), "0.0");
                        dcoefs2.push_back(z);
                    }
                    while (iter < max_iter && !converged)
                    {
                        auto result = aberth<mpc_60>(coefs2, dcoefs2, roots2);
                        std::vector<mpc_60> new_roots = result.first;
                        new_delta = result.second;
                        iter++;
                        for (unsigned j = 0; j < this->deg; ++j)
                        {
                            quadratic = true;
                            if (new_delta[j] > 0.01 * (delta[j] * delta[j]))
                            {
                                quadratic = false;
                                break;
                            }
                        }
                        within_tol = true;
                        for (unsigned j = 0; j < this->deg; ++j)
                        {
                            if (new_delta[j] > atol + rtol * abs(roots2[j]))
                            {
                                within_tol = false;
                                break;
                            }
                        }
                        if (quadratic || within_tol) converged = true;
                        roots2 = new_roots;
                        delta = new_delta;
                    }
                    for (unsigned i = 0; i < sharpen_iter; ++i)
                    {
                        auto result = newton<mpc_60>(coefs2, dcoefs2, roots2);
                        roots2 = result.first;
                    }
                    for (unsigned i = 0; i < this->deg; ++i)
                        roots[i] = static_cast<std::complex<double> >(roots2[i]);
                }
                else if (prec < 100)
                {
                    // Re-do the same calculations with mpc_100
                    prec = 100; iter = 0;
                    std::vector<mpc_100> roots2, coefs2, dcoefs2;
                    for (unsigned i = 0; i < this->deg; ++i)
                    {
                        roots2.push_back(mpc_100(roots[i].real(), roots[i].imag()));
                        delta[i] = 0.0;
                    }
                    for (unsigned i = 0; i < this->coefs.size(); ++i)
                    {
                        std::stringstream ss;
                        ss << std::setprecision(std::numeric_limits<T>::max_digits10) << this->coefs(i);
                        mpc_100 z(ss.str(), 0.0);
                        coefs2.push_back(z);
                    }
                    for (unsigned i = 0; i < dcoefvec.size(); ++i)
                    {
                        std::stringstream ss;
                        ss << std::setprecision(std::numeric_limits<T>::max_digits10) << dcoefvec(i);
                        mpc_100 z(ss.str(), "0.0");
                        dcoefs2.push_back(z);
                    }
                    while (iter < max_iter && !converged)
                    {
                        auto result = aberth<mpc_100>(coefs2, dcoefs2, roots2);
                        std::vector<mpc_100> new_roots = result.first;
                        new_delta = result.second;
                        iter++;
                        for (unsigned j = 0; j < this->deg; ++j)
                        {
                            quadratic = true;
                            if (new_delta[j] > 0.01 * (delta[j] * delta[j]))
                            {
                                quadratic = false;
                                break;
                            }
                        }
                        within_tol = true;
                        for (unsigned j = 0; j < this->deg; ++j)
                        {
                            if (new_delta[j] > atol + rtol * abs(roots2[j]))
                            {
                                within_tol = false;
                                break;
                            }
                        }
                        if (quadratic || within_tol) converged = true;
                        roots2 = new_roots;
                        delta = new_delta;
                    }
                    for (unsigned i = 0; i < sharpen_iter; ++i)
                    {
                        auto result = newton<mpc_100>(coefs2, dcoefs2, roots2);
                        roots2 = result.first;
                    }
                    for (unsigned i = 0; i < this->deg; ++i)
                        roots[i] = static_cast<std::complex<double> >(roots2[i]);
                }
                else if (prec < 200)
                {
                    // Re-do the same calculations with mpc_200
                    prec = 200; iter = 0;
                    std::vector<mpc_200> roots2, coefs2, dcoefs2;
                    for (unsigned i = 0; i < this->deg; ++i)
                    {
                        roots2.push_back(mpc_200(roots[i].real(), roots[i].imag()));
                        delta[i] = 0.0;
                    }
                    for (unsigned i = 0; i < this->coefs.size(); ++i)
                    {
                        std::stringstream ss;
                        ss << std::setprecision(std::numeric_limits<T>::max_digits10) << this->coefs(i);
                        mpc_200 z(ss.str(), 0.0);
                        coefs2.push_back(z);
                    }
                    for (unsigned i = 0; i < dcoefvec.size(); ++i)
                    {
                        std::stringstream ss;
                        ss << std::setprecision(std::numeric_limits<T>::max_digits10) << dcoefvec(i);
                        mpc_200 z(ss.str(), 0.0);
                        dcoefs2.push_back(z);
                    }
                    while (iter < max_iter && !converged)
                    {
                        auto result = aberth<mpc_200>(coefs2, dcoefs2, roots2);
                        std::vector<mpc_200> new_roots = result.first;
                        new_delta = result.second;
                        iter++;
                        for (unsigned j = 0; j < this->deg; ++j)
                        {
                            quadratic = true;
                            if (new_delta[j] > 0.01 * (delta[j] * delta[j]))
                            {
                                quadratic = false;
                                break;
                            }
                        }
                        within_tol = true;
                        for (unsigned j = 0; j < this->deg; ++j)
                        {
                            if (new_delta[j] > atol + rtol * abs(roots2[j]))
                            {
                                within_tol = false;
                                break;
                            }
                        }
                        if (quadratic || within_tol) converged = true;
                        roots2 = new_roots;
                        delta = new_delta;
                    }
                    for (unsigned i = 0; i < sharpen_iter; ++i)
                    {
                        auto result = newton<mpc_200>(coefs2, dcoefs2, roots2);
                        roots2 = result.first;
                    }
                    for (unsigned i = 0; i < this->deg; ++i)
                        roots[i] = static_cast<std::complex<double> >(roots2[i]);
                }
                else
                {
                    // Give up at this point (don't throw an exception here, 
                    // as polynomials with singular roots will be encountered)
                    converged = true;
                }
            }

            Matrix<std::complex<T>, Dynamic, 1> final_roots(this->deg);
            for (unsigned i = 0; i < this->deg; ++i)
            {
                std::complex<T> z(static_cast<T>(roots[i].real()), static_cast<T>(roots[i].imag()));
                final_roots(i) = z;
            }
            return final_roots;
        }


    public:
        Polynomial()
        {
            /*
             * Empty constructor for the zero polynomial.
             */
            this->deg = 0;
            this->coefs = Matrix<T, Dynamic, 1>::Zero(1);
            this->prec = std::numeric_limits<T>::max_digits10;
        }

        Polynomial(const T coef)
        {
            /*
             * Constructor with user-specified constant.
             */
            this->deg = 0;
            this->coefs = Matrix<T, Dynamic, 1>::Constant(1, 1, coef);
            this->prec = std::numeric_limits<T>::max_digits10;
        }

        Polynomial(const Ref<const Matrix<T, Dynamic, 1> >& coefs)
        {
            /*
             * Constructor with user-specified coefficients.
             */
            this->coefs = removeTrailingZeros<T>(coefs);
            this->deg = this->coefs.size() - 1;
            this->prec = std::numeric_limits<T>::max_digits10;
        }

        ~Polynomial()
        {
            /*
             * Empty destructor.
             */
        }

        unsigned degree() const
        {
            /*
             * Return the degree of this polynomial.
             */
            return this->deg;
        }

        Matrix<T, Dynamic, 1> coefficients() const
        {
            /*
             * Return the coefficients of this polynomial.
             */
            return this->coefs;
        }

        T eval(T x) const
        {
            /*
             * Return the value obtained by evaluating the polynomial at x.
             *
             * Horner's method: starting with the last (highest-degree)
             * coefficient, recursively do the following:
             * 1) Multiply by x
             * 2) Add by the previous coefficient
             */
            T y = this->coefs(this->deg);
            for (int i = this->deg - 1; i >= 0; --i)
            {
                y *= x;
                y += this->coefs(i);
            }
            return y;
        }

        std::complex<T> eval(std::complex<T> x) const
        {
            /*
             * Return the value obtained by evaluating the polynomial at x
             * with Horner's method.
             */
            std::complex<T> y(this->coefs(this->deg), 0.0);
            for (int i = this->deg - 1; i >= 0; --i)
            {
                y *= x;
                y += this->coefs(i);
            }
            return y;
        }

        Matrix<T, Dynamic, 1> eval(const Ref<const Matrix<T, Dynamic, 1> >& x) const
        {
            /*
             * Return the vector of values obtained by evaluating the polynomial
             * at each coordinate of the vector x.
             */
            Array<T, Dynamic, 1> y = Array<T, Dynamic, 1>::Constant(x.size(), 1, this->coefs(this->deg));
            for (int i = this->deg - 1; i >= 0; --i)
            {
                y *= x.array();
                y += this->coefs(i);
            }
            return y.matrix();
        }

        Polynomial operator+(const Polynomial<T>& q) const 
        {
            /*
             * Return the result of adding by q. 
             */
            // Copy over the coefficients into new vectors and resize 
            // as necessary 
            Matrix<T, Dynamic, 1> p_coefs(this->coefs);
            Matrix<T, Dynamic, 1> q_coefs(q.coefficients());
            if (this->deg > q.degree())
                q_coefs = appendZeros<T>(q_coefs, this->deg - q.degree());
            else if (this->deg < q.degree())
                p_coefs = appendZeros<T>(p_coefs, q.degree() - this->deg);

            // Instantiate the sum polynomial
            return Polynomial(p_coefs + q_coefs);
        }

        Polynomial operator+(const T s) const
        {
            /*
             * Return the result of adding by scalar s.
             */
            Matrix<T, Dynamic, 1> p_coefs(this->coefs);
            p_coefs(0) += s;
            return Polynomial(p_coefs);
        }

        Polynomial& operator+=(const Polynomial<T>& q)
        {
            /*
             * In-place addition by q.
             */
            // Copy over the coefficients into new vectors and resize 
            // as necessary 
            Matrix<T, Dynamic, 1> p_coefs(this->coefs);
            Matrix<T, Dynamic, 1> q_coefs(q.coefficients());
            if (this->deg > q.degree())
                q_coefs = appendZeros<T>(q_coefs, this->deg - q.degree());
            else if (this->deg < q.degree())
                p_coefs = appendZeros<T>(p_coefs, q.degree() - this->deg);

            // Update polynomial coefficients
            this->coefs = removeTrailingZeros<T>(p_coefs + q_coefs);
            this->deg = this->coefs.size() - 1;
            return *this;
        }

        Polynomial& operator+=(const T s)
        {
            /*
             * In-place addition by scalar s.
             */
            this->coefs(0) += s;
            return *this;
        }

        Polynomial operator-(const Polynomial<T>& q) const
        {
            /*
             * Return the result of subtracting by q. 
             */
            // Copy over the coefficients into new vectors and resize 
            // as necessary 
            Matrix<T, Dynamic, 1> p_coefs(this->coefs);
            Matrix<T, Dynamic, 1> q_coefs(q.coefficients());
            if (this->deg > q.degree())
                q_coefs = appendZeros<T>(q_coefs, this->deg - q.degree());
            else if (this->deg < q.degree())
                p_coefs = appendZeros<T>(p_coefs, q.degree() - this->deg);

            // Instantiate the difference polynomial
            return Polynomial(p_coefs - q_coefs);
        }

        Polynomial operator-(const T s) const
        {
            /*
             * Return the result of subtracting by scalar s.
             */
            Matrix<T, Dynamic, 1> p_coefs(this->coefs);
            p_coefs(0) -= s;
            return Polynomial(p_coefs);
        }

        Polynomial& operator-=(const Polynomial<T>& q) 
        {
            /*
             * In-place subtraction by q.
             */
            // Copy over the coefficients into new vectors and resize 
            // as necessary 
            Matrix<T, Dynamic, 1> p_coefs(this->coefs);
            Matrix<T, Dynamic, 1> q_coefs(q.coefficients());
            if (this->deg > q.degree())
                q_coefs = appendZeros<T>(q_coefs, this->deg - q.degree());
            else if (this->deg < q.degree())
                p_coefs = appendZeros<T>(p_coefs, q.degree() - this->deg);

            // Update polynomial coefficients
            this->coefs = removeTrailingZeros<T>(p_coefs - q_coefs);
            this->deg = this->coefs.size() - 1;
            return *this;
        }

        Polynomial& operator-=(const T s)
        {
            /*
             * In-place subtraction by scalar s.
             */
            this->coefs(0) -= s;
            return *this;
        }

        Polynomial operator*(const Polynomial<T>& q) const
        {
            /*
             * Return the result of multiplying by q.
             *
             * TODO Replace with a method based on FFT. 
             */
            unsigned p_deg = this->deg;
            unsigned q_deg = q.degree();
            Matrix<T, Dynamic, 1> q_coefs = q.coefficients();
            Matrix<T, Dynamic, 1> pq_coefs = Matrix<T, Dynamic, 1>::Zero(p_deg + q_deg + 1);
            for (unsigned i = 0; i <= p_deg; ++i)
            {
                for (unsigned j = 0; j <= q_deg; ++j)
                {
                    pq_coefs(i + j) += this->coefs(i) * q_coefs(j);
                }
            }
            return Polynomial(pq_coefs);
        }

        Polynomial operator*(const T s) const
        {
            /*
             * Return the result of multiplying by scalar s.
             */
            Matrix<T, Dynamic, 1> p_coefs(this->coefs);
            p_coefs *= s;
            return Polynomial(p_coefs);
        }

        Polynomial& operator*=(const Polynomial<T>& q) 
        {
            /*
             * In-place multiplication by q.
             *
             * TODO Replace with a method based on FFT. 
             */
            unsigned p_deg = this->deg;
            unsigned q_deg = q.degree();
            Matrix<T, Dynamic, 1> q_coefs = q.coefficients();
            Matrix<T, Dynamic, 1> pq_coefs = Matrix<T, Dynamic, 1>::Zero(p_deg + q_deg + 1);
            for (unsigned i = 0; i <= p_deg; ++i)
            {
                for (unsigned j = 0; j <= q_deg; ++j)
                {
                    pq_coefs(i + j) += this->coefs(i) * q_coefs(j);
                }
            }

            // Update polynomial coefficients
            this->coefs = removeTrailingZeros<T>(pq_coefs);
            this->deg = pq_coefs.size() - 1;
            return *this;
        }

        Polynomial& operator*=(const T s)
        {
            /*
             * In-place multiplication by scalar s.
             */
            this->coefs *= s;
            return *this;
        }

        Polynomial operator/(const T s) const
        {
            /*
             * Return the result of dividing by scalar s.
             */
            Matrix<T, Dynamic, 1> p_coefs(this->coefs);
            p_coefs /= s;
            return Polynomial(p_coefs);
        }

        Polynomial& operator/=(const T s)
        {
            /*
             * In-place division by scalar s.
             */
            this->coefs /= s;
            return *this;
        }

        Polynomial diff() const 
        {
            /*
             * Return the derivative of this polynomial.
             */
            if (this->deg == 0)
                return Polynomial();

            Matrix<T, Dynamic, 1> dcoefs = (
                Array<T, Dynamic, 1>::LinSpaced(this->deg, 1, this->deg)
                * this->coefs.tail(this->deg).array()
            ).matrix();
            return Polynomial(dcoefs);
        }

        T leadingCoef() const
        {
            /*
             * Return the leading coefficient of this polynomial.
             */
            return this->coefs[this->deg];
        }

        Polynomial<T> reduce() const
        {
            /*
             * Factor out all powers of x from this polynomial. 
             */
            Matrix<T, Dynamic, 1> reduced(this->coefs);
            while (reduced(0) == 0.0)
                reduced = reduced.tail(reduced.size() - 1).eval();
            return Polynomial<T>(reduced);
        }

        Polynomial<T> monic() const
        {
            /*
             * Factor out the leading coefficient. 
             */
            return Polynomial<T>(this->coefs / this->leadingCoef());
        }

        Matrix<std::complex<T>, Dynamic, 1> roots(SolveMethod method = Aberth,
                                                  unsigned max_iter = 1000,
                                                  double atol = 1e-15,
                                                  double rtol = 1e-15,
                                                  unsigned sharpen_iter = 20,
                                                  std::vector<std::complex<double> > inits = {})
        {
            /*
             * Return all complex roots of the polynomial.  
             */
            if (this->deg == 2)
            {
                Matrix<std::complex<T>, Dynamic, 1> roots(2);
                T b = this->coefs(1) / this->coefs(2);
                T c = this->coefs(0) / this->coefs(2);
                std::vector<std::complex<T> > r = quadratic<T>(b, c);
                roots(0) = r[0];
                roots(1) = r[1];
                return roots;
            }
            //else if (this->deg == 3)
            //{
            //    Matrix<std::complex<T>, Dynamic, 1> roots(3);
            //    T b = this->coefs(2) / this->coefs(3);
            //    T c = this->coefs(1) / this->coefs(3);
            //    T d = this->coefs(0) / this->coefs(3);
            //    std::vector<std::complex<T> > r = cubic<T>(b, c, d);
            //    roots(0) = r[0];
            //    roots(1) = r[1];
            //    roots(2) = r[2];
            //    return roots;
            //}
            else if (method == Aberth)
            {
                return this->rootsAberth(max_iter, atol, rtol, sharpen_iter, inits);
            }
            else
            {
                return this->rootsWeierstrass(max_iter, atol, rtol, inits);
            }
        }
};

#endif
