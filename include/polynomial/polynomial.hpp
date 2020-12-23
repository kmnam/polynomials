#ifndef UNIVAR_POLYNOMIAL_MULTIPREC_HPP
#define UNIVAR_POLYNOMIAL_MULTIPREC_HPP

#include <utility>
#include <ostream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <boost/random.hpp>
#include <boost/multiprecision/mpc.hpp>

/*
 * A Polynomial class template with arbitrary real coefficients.
 *
 * The internal functions defined here use boost::multiprecision::number
 * types with float/complex backends with user-defined precision.  
 *
 * Authors:
 *     Kee-Myoung Nam, Department of Systems Biology, Harvard Medical School
 * Last updated:
 *     12/22/2020
 */
using boost::multiprecision::number;
using boost::multiprecision::mpfr_float_backend;
using boost::multiprecision::mpc_complex_backend;

// ------------------------------------------------------------------ //
//                      VARIOUS HELPER FUNCTIONS                      //
// ------------------------------------------------------------------ //
enum SolveMethod 
{
    Weierstrass,
    Aberth
};

template <typename T>
std::vector<T> appendZeros(const std::vector<T>& v, int j)
{
    /*
     * Return a new copy of v with j appended zeros.
     */
    std::vector<T> w(v);
    for (int i = 0; i < j; ++i) w.push_back(0);
    return w;
}

template <typename T>
std::vector<T> removeTrailingZeros(const std::vector<T>& v)
{
    /*
     * Return a new copy of v with all trailing zeros removed. 
     */
    std::vector<T> w(v);
    for (int i = v.size() - 1; i >= 0; --i)
    {
        if (v[i] == 0) w.pop_back();
        else           break;
    }
    return w;
}

template <typename T, unsigned N>
number<mpfr_float_backend<N> > floatToNumber(const T x)
{
    /*
     * Convert a floating-point number to a number of precision N via 
     * a stringstream. 
     */
    std::stringstream ss; 

    // Write the string to precision set by std::numeric_limits<T>::max_digits10
    ss << std::setprecision(std::numeric_limits<T>::max_digits10) << x;

    return number<mpfr_float_backend<N> >(ss.str());
}

template <unsigned N, unsigned M>
number<mpfr_float_backend<M> > convertPrecision(const number<mpfr_float_backend<N> > x)
{
    /*
     * Convert a number of precision N to a number of precision M via 
     * a stringstream. 
     */
    std::stringstream ss;

    // Write the string to precision N (which is fewer digits than that 
    // set by std::numeric_limits<T>::max_digits10) 
    ss << std::setprecision(N) << x;

    return number<mpfr_float_backend<M> >(ss.str());
}

template <typename T, unsigned N>
number<mpc_complex_backend<N> > floatToNumber(const std::complex<T> x)
{
    /*
     * Convert a floating-point std::complex number to a number of precision N
     * via a stringstream. 
     */
    std::stringstream ss_real, ss_imag; 

    // Write the real and imaginary parts of x to precision set by
    // std::numeric_limits<T>::max_digits10
    ss_real << std::setprecision(std::numeric_limits<T>::max_digits10) << x.real();
    ss_imag << std::setprecision(std::numeric_limits<T>::max_digits10) << x.imag();

    return number<mpc_complex_backend<N> >(ss_real.str(), ss_imag.str());
}

template <unsigned N, unsigned M>
number<mpc_complex_backend<M> > convertPrecision(const number<mpc_complex_backend<N> > x)
{
    /*
     * Convert a number of precision N to a number of precision M via 
     * a stringstream. 
     */
    std::stringstream ss_real, ss_imag;

    // Write the parts to precision N (which is fewer digits than that 
    // set by std::numeric_limits<T>::max_digits10) 
    ss_real << std::setprecision(N) << x.real();
    ss_imag << std::setprecision(N) << x.imag();

    return number<mpc_complex_backend<M> >(ss_real.str(), ss_imag.str());
}

template <unsigned N>
std::vector<number<mpc_complex_backend<N> > > quadratic(const number<mpfr_float_backend<N> > b,
                                                        const number<mpfr_float_backend<N> > c)
{
    /*
     * Solve the given monic quadratic, x^2 + b*x + c = 0, with coefficients of 
     * type number<mpfr_float_backend>.
     */
    typedef number<mpc_complex_backend<N> > CTN;
    std::vector<CTN> roots; 
    CTN disc = boost::multiprecision::sqrt(CTN(b * b - 4 * c, 0.0));
    roots.push_back((-b + disc) * 0.5);
    roots.push_back((-b - disc) * 0.5);
    return roots;
}

template <unsigned N>
number<mpc_complex_backend<N> > cuberoot(number<mpc_complex_backend<N> > z)
{
    /*
     * Return the principal cube root of the complex number z.
     */
    const number<mpfr_float_backend<N> > one_third = number<mpfr_float_backend<N> >("1.0") / number<mpfr_float_backend<N> >("3.0");
    if (z.real() < 0)
        return -boost::multiprecision::pow(-z, one_third);
    else
        return boost::multiprecision::pow(z, one_third);
}

template <unsigned N>
std::vector<number<mpc_complex_backend<N> > > cubic(const number<mpfr_float_backend<N> > b,
                                                    const number<mpfr_float_backend<N> > c,
                                                    const number<mpfr_float_backend<N> > d)
{
    /*
     * Solve the given monic cubic, x^3 + b*x^2 + c*x + d = 0. 
     */
    typedef number<mpfr_float_backend<N> >  RTN;
    typedef number<mpc_complex_backend<N> > CTN;
    std::vector<CTN> roots;
    const RTN sqrt3 = boost::multiprecision::sqrt(RTN("3.0"));
    RTN q = (3 * c - b * b) / 9;
    RTN r = (9 * b * c - 27 * d - 2 * boost::multiprecision::pow(b, 3)) / 54;
    RTN p = boost::multiprecision::pow(q, 3) + boost::multiprecision::pow(r, 2);
    CTN s = cuberoot<N>(r + boost::multiprecision::sqrt(CTN(p, 0.0)));
    CTN t = cuberoot<N>(r - boost::multiprecision::sqrt(CTN(p, 0.0)));
    roots.push_back(-b / 3.0 + s + t);
    roots.push_back(-b / 3.0 - (s + t) / 2.0 + CTN(0.0, 1.0) * sqrt3 * (s - t) / 2.0);
    roots.push_back(boost::multiprecision::conj(roots[1]));
    return roots;
}

template <unsigned N>
number<mpc_complex_backend<N> > horner(const std::vector<number<mpfr_float_backend<N> > >& coefs,
                                       const number<mpc_complex_backend<N> > x)
{
    /*
     * Perform Horner's method for evaluating a (real) polynomial at a
     * (complex) value, given the vector of its coefficients and the value
     * at which the polynomial is to be evaluated.
     *
     * This function assumes that the input precision equals the output 
     * precision, and thus does not call convertPrecision().  
     */
    typedef number<mpc_complex_backend<N> > CT;

    CT y(coefs[coefs.size() - 1], 0.0);
    for (int i = coefs.size() - 2; i >= 0; --i)
    {
        y *= x;
        CT coef(coefs[i], 0.0);
        y += coef;
    }
    return y;
}

template <unsigned N>
number<mpfr_float_backend<N> > horner(const std::vector<number<mpfr_float_backend<N> > >& coefs,
                                      const number<mpfr_float_backend<N> > x)
{
    /*
     * Perform Horner's method for evaluating a (real) polynomial at a 
     * (real) value, given the vector of its coefficients and the value 
     * at which the polynomial is to be evaluated.
     *
     * This function assumes that the input precision equals the output
     * precision, and thus does not call convertPrecision(). 
     */
    typedef number<mpfr_float_backend<N> > RT; 

    RT y = coefs[coefs.size() - 1];
    for (int i = coefs.size() - 2; i >= 0; --i)
    {
        y *= x;
        y += coefs[i]; 
    }
    return y; 
}

template <unsigned N, unsigned M>
number<mpc_complex_backend<M> > hornerWiden(const std::vector<number<mpfr_float_backend<N> > >& coefs,
                                            const number<mpc_complex_backend<N> > x)
{
    /*
     * Perform Horner's method for evaluating a (real) polynomial at a
     * (complex) value, given the vector of its coefficients and the value
     * at which the polynomial is to be evaluated. 
     */
    typedef number<mpfr_float_backend<M> >  RTM;
    typedef number<mpc_complex_backend<M> > CTM;

    CTM value = convertPrecision<N, M>(x);
    RTM yr = convertPrecision<N, M>(coefs[coefs.size() - 1]);
    CTM y(yr, 0.0);
    for (int i = coefs.size() - 2; i >= 0; --i)
    {
        y *= value;
        RTM cr = convertPrecision<N, M>(coefs[i]); 
        CTM coef(cr, 0.0);
        y += coef;
    }
    return y;
}

template <unsigned N, unsigned M>
number<mpfr_float_backend<M> > hornerWiden(const std::vector<number<mpfr_float_backend<N> > >& coefs,
                                           const number<mpfr_float_backend<N> > x)
{
    /*
     * Perform Horner's method for evaluating a (real) polynomial at a
     * (real) value, given the vector of its coefficients and the value
     * at which the polynomial is to be evaluated. 
     */
    typedef number<mpfr_float_backend<M> >  RTM;

    RTM value = convertPrecision<N, M>(x);
    RTM y = convertPrecision<N, M>(coefs[coefs.size() - 1]);
    for (int i = coefs.size() - 2; i >= 0; --i)
    {
        y *= value;
        y += convertPrecision<N, M>(coefs[i]);
    }
    return y; 
}

template <unsigned N, unsigned M>
std::pair<std::vector<number<mpc_complex_backend<M> > >, std::vector<double> >
    newton(const std::vector<number<mpfr_float_backend<N> > >& coefs,
           const std::vector<number<mpfr_float_backend<N> > >& dcoefs,
           const std::vector<number<mpc_complex_backend<N> > >& roots)
{
    /*
     * Perform one iteration of Newton's method, given the vector of 
     * coefficients corresponding to the polynomial, its derivative, and 
     * a vector of current root approximations. 
     */
    typedef number<mpc_complex_backend<M> > CTM;

    // Perform the Newton correction for each root
    std::vector<CTM> new_roots;
    for (int j = 0; j < roots.size(); ++j)
    {
        CTM new_root = convertPrecision<N, M>(roots[j]);
        new_roots.push_back(new_root - (hornerWiden<N, M>(coefs, roots[j]) / hornerWiden<N, M>(dcoefs, roots[j])));
    }

    // Compute the absolute difference between each old and new root
    std::vector<double> delta;
    for (int j = 0; j < roots.size(); ++j)
    {
        double dr = static_cast<double>(roots[j].real()) - static_cast<double>(new_roots[j].real());
        double di = static_cast<double>(roots[j].imag()) - static_cast<double>(new_roots[j].imag());
        delta.push_back(std::sqrt(std::pow(dr, 2.0) + std::pow(di, 2.0)));
    }

    return std::make_pair(new_roots, delta);
}

template <unsigned N, unsigned M>
std::pair<std::vector<number<mpc_complex_backend<M> > >, std::vector<double> >
    weierstrass(const std::vector<number<mpfr_float_backend<N> > >& coefs,
                const std::vector<number<mpc_complex_backend<N> > >& roots)
{
    /*
     * Perform one iteration of Weierstrass' method, given the vector 
     * of coefficients corresponding to the polynomial and a vector of 
     * current root approximations.  
     */
    typedef number<mpc_complex_backend<M> > CTM;

    // Perform the Weierstrass correction for each root
    std::vector<CTM> new_roots;
    for (auto&& root : roots)
    {
        CTM new_root = convertPrecision<N, M>(root);
        new_roots.push_back(new_root);
    }
    for (int j = 0; j < roots.size(); ++j)
    {
        CTM root = convertPrecision<N, M>(roots[j]);
        CTM denom = 1.0;
        for (int k = 0; k < roots.size(); ++k)
        {
            if (j != k) denom *= (root - new_roots[k]);
        }
        CTM correction = -hornerWiden<N, M>(coefs, roots[j]) / denom;
        new_roots[j] += correction;
    }

    // Compute the absolute difference between each old and new root
    std::vector<double> delta;
    for (int j = 0; j < roots.size(); ++j)
    {
        double dr = static_cast<double>(roots[j].real()) - static_cast<double>(new_roots[j].real());
        double di = static_cast<double>(roots[j].imag()) - static_cast<double>(new_roots[j].imag());
        delta.push_back(std::sqrt(std::pow(dr, 2.0) + std::pow(di, 2.0)));
    }

    return std::make_pair(new_roots, delta);
}

template <unsigned N, unsigned M>
std::pair<std::vector<number<mpc_complex_backend<M> > >, std::vector<double> >
    aberth(const std::vector<number<mpfr_float_backend<N> > >& coefs,
           const std::vector<number<mpfr_float_backend<N> > >& dcoefs,
           const std::vector<number<mpc_complex_backend<N> > >& roots)
{
    /*
     * Perform one iteration of the Aberth-Ehrlich method, given the vector 
     * of coefficients corresponding to the polynomial and a vector of 
     * current root approximations.  
     */
    typedef number<mpc_complex_backend<M> > CTM;

    // Perform the Aberth-Ehrlich correction for each root
    std::vector<CTM> new_roots;
    for (auto&& root : roots)
    {
        CTM new_root = convertPrecision<N, M>(root);
        new_roots.push_back(new_root);
    }
    for (int j = 0; j < roots.size(); ++j)
    {
        CTM root_j = convertPrecision<N, M>(roots[j]);
        CTM sum = 0.0;
        for (int k = 0; k < j; ++k)
        {
            sum += (1.0 / (root_j - new_roots[k]));
        }
        for (int k = j + 1; k < roots.size(); ++k)
        {
            CTM root_k(roots[k]);
            sum += (1.0 / (root_j - root_k));
        }
        CTM value = hornerWiden<N, M>(coefs, roots[j]) / hornerWiden<N, M>(dcoefs, roots[j]);
        CTM denom = 1.0 - (value * sum);
        CTM correction = value / denom;
        new_roots[j] -= correction;
    }

    // Compute the absolute difference between each old and new root
    std::vector<double> delta;
    for (int j = 0; j < roots.size(); ++j)
    {
        double dr = static_cast<double>(roots[j].real()) - static_cast<double>(new_roots[j].real());
        double di = static_cast<double>(roots[j].imag()) - static_cast<double>(new_roots[j].imag());
        delta.push_back(std::sqrt(std::pow(dr, 2.0) + std::pow(di, 2.0)));
    }
    
    return std::make_pair(new_roots, delta);
}

bool ccw(std::pair<double, double> x, std::pair<double, double> y, std::pair<double, double> z)
{
    /*
     * Determine if x, y, and z make a counterclockwise turn in 2-D. 
     */
    double cross = (y.first - x.first) * (z.second - x.second) - (y.second - x.second) * (z.first - x.first);
    return (cross > 0);
}

template <unsigned N, unsigned M>
std::vector<number<mpc_complex_backend<M> > > bini(const std::vector<number<mpfr_float_backend<N> > > coefs,
                                                   boost::random::mt19937& rng,
                                                   boost::random::uniform_real_distribution<double>& dist)
{
    /*
     * Use Bini's initialization to yield a set of initial complex roots 
     * for the given vector of polynomial coefficients. 
     */
    typedef number<mpfr_float_backend<N> >  RTM;
    typedef number<mpc_complex_backend<M> > CTM;

    // Find the upper convex hull of the vertices (i, log|f_i|) with 
    // the Andrew scan (use doubles for this computation)
    std::vector<std::pair<double, double> > points;
    for (int i = 0; i < coefs.size(); ++i)
    {
        if (boost::multiprecision::abs(coefs[i]) == 0)
            points.push_back(std::make_pair(i, -std::numeric_limits<double>::infinity()));
        else
            points.push_back(std::make_pair(i, static_cast<double>(boost::multiprecision::log2(boost::multiprecision::abs(coefs[i])))));
    }
    std::vector<std::pair<double, double> > upper_hull;
    for (int i = points.size() - 1; i >= 0; --i)
    {
        while (upper_hull.size() >= 2 && !ccw(upper_hull[upper_hull.size()-2], upper_hull[upper_hull.size()-1], points[i]))
            upper_hull.pop_back();
        upper_hull.push_back(points[i]);
    }

    // Run through the convex hull vertices, grab their x-coordinates, and sort
    std::vector<int> hull_x;
    for (auto&& p : upper_hull)
    {
        if (!std::isinf(p.second))
            hull_x.push_back(static_cast<int>(p.first));
    }
    std::sort(hull_x.begin(), hull_x.end());

    // Compute u_{k_i} for each vertex in the convex hull
    std::vector<double> u;
    for (int i = 0; i < hull_x.size() - 1; ++i)
    {
        int k = hull_x[i];
        int l = hull_x[i+1];
        u.push_back(std::pow(static_cast<double>(boost::multiprecision::abs(coefs[k] / coefs[l])), 1.0 / (l - k)));
    }

    // Compute initial root approximations
    // TODO Check that these approximations are correct (how?) on small examples
    int n = coefs.size() - 1;
    const CTM two_pi = floatToNumber<double, M>(std::complex<double>(2 * std::acos(-1.0), 0)); 
    std::vector<CTM> inits; 
    for (int i = 0; i < hull_x.size() - 1; ++i)
    {
        int k = hull_x[i];
        int l = hull_x[i+1];
        for (int j = 0; j < l - k; ++j)
        {
            RTM theta = floatToNumber<double>((two_pi / (l - k)) * j + (two_pi * i / n) + dist(rng));
            CTM z = floatToNumber<double, M>(std::complex<double>(std::cos(theta), std::sin(theta)));
            inits.push_back(floatToNumber<double, M>(u[i]) * z);
        }
    }
    return inits;
}

template <unsigned N>
class Polynomial
{
    private:
        int deg;                    // Degree of the polynomial
        
        // Coefficients of the polynomial stored in ascending order of power
        std::vector<number<mpfr_float_backend<N> > > coefs;

        template <unsigned M = N>
        std::pair<std::vector<number<mpc_complex_backend<M> > >, bool> rootsWeierstrass(int max_iter,
                                                                                        double atol,
                                                                                        double rtol)
        {
            /*
             * Run Weierstrass' method on the given polynomial. 
             */
            typedef number<mpfr_float_backend<M> >  RTM;
            typedef number<mpc_complex_backend<M> > CTM;

            bool converged = false;
            bool quadratic = false;
            bool within_tol = false;
            int iter = 0;

            // Initialize the roots to (0.4 + 0.9 i)^p, for p = 0, ..., d - 1
            // TODO Fix this to work with Bini's initialization?
            std::vector<CTM> inits;
            for (int i = 0; i < this->deg; ++i)
            {
                CTM init = boost::multiprecision::pow(floatToNumber<double, M>(std::complex<double>(0.4, 0.9)), i);
                inits.push_back(init);
            }
           
            // Set up vectors of coefficients and roots at given precision (M) 
            std::vector<RTM> coefs;
            for (auto&& coef : this->coefs) coefs.push_back(convertPrecision<N, M>(coef));
            std::vector<CTM> roots(inits);
            std::vector<CTM> new_roots;

            // Run the algorithm until all roots exhibit quadratic convergence
            // or the desired number of iterations is reached
            std::vector<double> delta, new_delta;
            for (int i = 0; i < this->deg; ++i) delta.push_back(0.0);
            while (iter < max_iter && !converged)
            {
                std::pair<std::vector<CTM>, std::vector<double> > result = weierstrass<M, M>(coefs, roots);
                new_roots = result.first;
                new_delta = result.second;
                iter++;

                // Check that the change is less than the square of the 
                // previous change for each root
                quadratic = true;
                for (int j = 0; j < this->deg; ++j)
                {
                    if (new_delta[j] > 0.01 * (delta[j] * delta[j]))
                    {
                        quadratic = false;
                        break;
                    }
                }
                // Check that the change is within the given tolerances
                within_tol = true;
                for (int j = 0; j < this->deg; ++j)
                {
                    if (new_delta[j] > atol + rtol * boost::multiprecision::abs(roots[j]))
                    {
                        within_tol = false;
                        break;
                    }
                }
                if (quadratic || within_tol) converged = true; // TODO Want both?
                roots = new_roots;
                delta = new_delta;
            }

            return std::make_pair(roots, converged);
        }

        template <unsigned M = N>
        std::pair<std::vector<number<mpc_complex_backend<M> > >, bool> rootsAberth(int max_iter,
                                                                                   double atol,
                                                                                   double rtol)
        {
            /*
             * Run the Aberth-Ehrlich method on the given polynomial.
             */
            typedef number<mpfr_float_backend<M> >  RTM;
            typedef number<mpc_complex_backend<M> > CTM;

            bool converged = false;
            bool quadratic = false;
            bool within_tol = false;
            int iter = 0;

            // Initialize the roots to (0.4 + 0.9 i)^p, for p = 0, ..., d - 1
            // TODO Fix this to work with Bini's initialization?
            std::vector<CTM> inits;
            for (int i = 0; i < this->deg; ++i)
            {
                CTM init = boost::multiprecision::pow(floatToNumber<double, M>(std::complex<double>(0.4, 0.9)), i);
                inits.push_back(init);
            }
           
            // Set up vectors of coefficients and roots at given precision (M)
            std::vector<RTM> coefs, dcoefs;
            for (int i = 0; i < this->coefs.size(); ++i)
                coefs.push_back(convertPrecision<N, M>(coefs[i]));
            for (int i = 0; i < this->coefs.size() - 1; ++i)
                dcoefs.push_back((i + 1) * coefs[i+1]);
            std::vector<CTM> roots(inits);
            std::vector<CTM> new_roots;

            // Run the algorithm until all roots exhibit quadratic convergence
            // or the desired number of iterations is reached
            std::vector<double> delta, new_delta;
            for (int i = 0; i < this->deg; ++i) delta.push_back(0.0);
            while (iter < max_iter && !converged)
            {
                std::pair<std::vector<CTM>, std::vector<double> > result = aberth<M, M>(coefs, dcoefs, roots);
                new_roots = result.first;
                new_delta = result.second;
                iter++;

                // Check that the change is less than the square of the 
                // previous change for each root
                for (int j = 0; j < this->deg; ++j)
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
                for (int j = 0; j < this->deg; ++j)
                {
                    if (new_delta[j] > atol + rtol * boost::multiprecision::abs(roots[j]))
                    {
                        within_tol = false;
                        break;
                    }
                }
                if (quadratic || within_tol) converged = true;  // TODO Want both?
                roots = new_roots;
                delta = new_delta;
            }
            // Sharpen each root for the given number of iterations 
            //for (unsigned i = 0; i < sharpen_iter; ++i)
            //{
            //    auto result = newton<std::complex<double> >(coefs, dcoefs, roots);
            //    roots = result.first;
            //}

            return std::make_pair(roots, converged);
        }


    public:
        Polynomial()
        {
            /*
             * Empty constructor for the zero polynomial.
             */
            this->deg = 0;
            this->coefs.push_back(0.0);
        }

        Polynomial(const double coef)
        {
            /*
             * Constructor with user-specified double constant.
             */
            this->deg = 0;
            this->coefs.push_back(floatToNumber<double, N>(coef));
        }

        Polynomial(const number<mpfr_float_backend<N> > coef)
        {
            /*
             * Constructor with user-specified constant.
             */
            this->deg = 0;
            this->coefs.push_back(coef);
        }

        Polynomial(const std::vector<double> coefs)
        {
            /*
             * Constructor with user-specified vector of double coefficients.
             */
            typedef number<mpfr_float_backend<N> > RTN;

            for (auto&& coef : coefs) this->coefs.push_back(floatToNumber<double, N>(coef));
            this->coefs = removeTrailingZeros<RTN>(this->coefs);
            this->deg = this->coefs.size() - 1;
        }

        Polynomial(const std::vector<number<mpfr_float_backend<N> > > coefs)
        {
            /*
             * Constructor with user-specified vector of coefficients.
             */
            this->coefs = removeTrailingZeros<number<mpfr_float_backend<N> > >(coefs);
            this->deg = this->coefs.size() - 1;
        }

        ~Polynomial()
        {
            /*
             * Empty destructor.
             */
        }

        int degree() const
        {
            /*
             * Return the degree of this polynomial.
             */
            return this->deg;
        }

        std::vector<number<mpfr_float_backend<N> > > coefficients() const
        {
            /*
             * Return the coefficients of this polynomial.
             */
            return this->coefs;
        }

        number<mpfr_float_backend<N> > eval(const number<mpfr_float_backend<N> > x) const 
        {
            /*
             * Return the value obtained by evaluating the polynomial at the 
             * given (real) value x.
             *
             * Here, precision of the input value matches that of the polynomial.   
             */
            return horner<N>(this->coefs, x);
        }

        template <unsigned M, unsigned P>
        number<mpfr_float_backend<P> > eval(const number<mpfr_float_backend<M> > x) const
        {
            /*
             * Return the value obtained by evaluating the polynomial at the
             * given (real) value x.
             */
            return hornerWiden<N, P>(this->coefs, convertPrecision<M, N>(x));
        }

        template <unsigned M = N>
        number<mpfr_float_backend<M> > eval(const double x) const
        {
            /*
             * Return the value obtained by evaluating the polynomial at the
             * given (real double) value x. 
             */
            if (M == N)
                return horner<N>(this->coefs, floatToNumber<double, N>(x));
            else
                return hornerWiden<N, M>(this->coefs, floatToNumber<double, N>(x));
        }

        number<mpc_complex_backend<N> > eval(const number<mpc_complex_backend<N> > x) const
        {
            /*
             * Return the value obtained by evaluating the polynomial at the 
             * given (complex) value x. 
             *
             * Here, precision of the input value matches that of the polynomial.   
             */
            return horner<N>(this->coefs, x); 
        }

        template <unsigned M, unsigned P>
        number<mpc_complex_backend<P> > eval(const number<mpc_complex_backend<M> > x) const
        {
            /*
             * Return the value obtained by evaluating the polynomial at the
             * given (complex) value x.
             */
            return hornerWiden<N, P>(this->coefs, convertPrecision<M, N>(x));
        }

        template <unsigned M = N>
        number<mpc_complex_backend<M> > eval(const std::complex<double> x) const 
        {
            /*
             * Return the value obtained by evaluating the polynomial at the 
             * given (std::complex<double>) value x. 
             */
            if (M == N)
                return horner<N>(this->coefs, floatToNumber<double, N>(x));
            else
                return hornerWiden<N, M>(this->coefs, floatToNumber<double, N>(x));
        }

        std::vector<number<mpfr_float_backend<N> > > eval(const std::vector<number<mpfr_float_backend<N> > >& x) const 
        {
            /*
             * Return the vector of values obtained by evaluating the polynomial 
             * at each (real) value in the given vector x.
             *
             * Here, precision of the input values matches that of the polynomial.   
             */
            std::vector<number<mpfr_float_backend<N> > > values;
            for (auto&& y : x)
                values.push_back(this->eval(y));
            
            return values;
        }

        template <unsigned M, unsigned P>
        std::vector<number<mpfr_float_backend<P> > > eval(const std::vector<number<mpfr_float_backend<M> > >& x) const
        {
            /*
             * Return the vector of values obtained by evaluating the polynomial
             * at each (real) value in the given vector x.
             */
            std::vector<number<mpfr_float_backend<P> > > values;
            for (auto&& y : x)
                values.push_back(this->eval<M, P>(y));
            
            return values;
        }

        std::vector<number<mpc_complex_backend<N> > > eval(const std::vector<number<mpc_complex_backend<N> > >& x) const 
        {
            /*
             * Return the vector of values obtained by evaluating the polynomial 
             * at each (complex) value in the given vector x.
             *
             * Here, precision of the input values matches that of the polynomial.   
             */
            std::vector<number<mpc_complex_backend<N> > > values;
            for (auto&& y : x)
                values.push_back(this->eval(y));
            
            return values;
        }

        template <unsigned M, unsigned P>
        std::vector<number<mpc_complex_backend<P> > > eval(const std::vector<number<mpc_complex_backend<M> > >& x) const
        {
            /*
             * Return the vector of values obtained by evaluating the polynomial
             * at each (complex) value in the given vector x.
             */
            std::vector<number<mpc_complex_backend<P> > > values;
            for (auto&& y : x)
                values.push_back(this->eval<M, P>(y));
            
            return values;
        }

        Polynomial<N> operator+(const Polynomial<N>& q) const 
        {
            /*
             * Return the result of adding by the same-precision polynomial q. 
             */
            typedef number<mpfr_float_backend<N> > RTN;

            // Copy over the coefficients into new vectors and resize as necessary 
            std::vector<RTN> p_coefs, q_coefs;
            p_coefs = this->coefs;
            q_coefs = q.coefficients(); 
            if (this->deg > q.degree())
                q_coefs = appendZeros<RTN>(q_coefs, this->deg - q.degree());
            else if (this->deg < q.degree())
                p_coefs = appendZeros<RTN>(p_coefs, q.degree() - this->deg);

            // Instantiate the sum polynomial
            std::vector<RTN> sum_coefs;
            for (int i = 0; i < p_coefs.size(); ++i)
                sum_coefs.push_back(p_coefs[i] + q_coefs[i]);

            return Polynomial<N>(sum_coefs);
        }

        template <unsigned M, unsigned P>
        Polynomial<P> operator+(const Polynomial<M>& q) const 
        {
            /*
             * Return the result of adding by q. 
             */
            typedef number<mpfr_float_backend<N> > RTN;
            typedef number<mpfr_float_backend<M> > RTM; 
            typedef number<mpfr_float_backend<P> > RTP;

            // Copy over the coefficients into new vectors with the output
            // precision (P) and resize as necessary 
            std::vector<RTP> p_coefs, q_coefs;
            for (const RTN coef : this->coefs)      p_coefs.push_back(convertPrecision<N, P>(coef));
            for (const RTM coef : q.coefficients()) q_coefs.push_back(convertPrecision<M, P>(coef));
            if (this->deg > q.degree())
                q_coefs = appendZeros<RTP>(q_coefs, this->deg - q.degree());
            else if (this->deg < q.degree())
                p_coefs = appendZeros<RTP>(p_coefs, q.degree() - this->deg);

            // Instantiate the sum polynomial
            std::vector<RTP> sum_coefs;
            for (int i = 0; i < p_coefs.size(); ++i)
                sum_coefs.push_back(p_coefs[i] + q_coefs[i]);

            return Polynomial<P>(sum_coefs);
        }

        Polynomial<N> operator+(const number<mpfr_float_backend<N> > s) const
        {
            /*
             * Return the result of adding by the same-precision scalar s.
             */
            std::vector<number<mpfr_float_backend<N> > > sum_coefs;
            sum_coefs.push_back(this->coefs[0] + s);
            for (int i = 1; i < this->coefs.size(); ++i)
                sum_coefs.push_back(this->coefs[i]);

            return Polynomial<N>(sum_coefs);
        }

        template <unsigned M, unsigned P>
        Polynomial<P> operator+(const number<mpfr_float_backend<M> > s) const
        {
            /*
             * Return the result of adding by scalar s.
             */
            typedef number<mpfr_float_backend<P> > RTP;

            std::vector<RTP> sum_coefs;
            sum_coefs.push_back(convertPrecision<N, P>(this->coefs[0]) + convertPrecision<M, P>(s));
            for (int i = 1; i < this->coefs.size(); ++i)
                sum_coefs.push_back(convertPrecision<N, P>(this->coefs[i]));

            return Polynomial<P>(sum_coefs);
        }

        Polynomial& operator+=(const Polynomial<N>& q)
        {
            /*
             * In-place addition by the same-precision polynomial q.
             */
            typedef number<mpfr_float_backend<N> > RTN;

            // Copy over the coefficients into new vectors and resize 
            // as necessary 
            std::vector<RTN> p_coefs(this->coefs);
            std::vector<RTN> q_coefs(q.coefficients());
            if (this->deg > q.degree())
                q_coefs = appendZeros<RTN>(q_coefs, this->deg - q.degree());
            else if (this->deg < q.degree())
                p_coefs = appendZeros<RTN>(p_coefs, q.degree() - this->deg);

            // Update polynomial coefficients
            for (int i = 0; i < p_coefs.size(); ++i)
                p_coefs[i] += q_coefs[i];
            this->coefs = removeTrailingZeros<RTN>(p_coefs);
            this->deg = this->coefs.size() - 1;
            return *this;
        }

        Polynomial& operator+=(const number<mpfr_float_backend<N> > s)
        {
            /*
             * In-place addition by the same-precision scalar s.
             *
             * Note that the precision of the input scalar must match that 
             * of the current polynomial. 
             */
            this->coefs[0] += s;
            return *this;
        }

        Polynomial<N> operator-(const Polynomial<N>& q) const 
        {
            /*
             * Return the result of subtracting by the same-precision polynomial q. 
             */
            typedef number<mpfr_float_backend<N> > RTN;

            // Copy over the coefficients into new vectors and resize as necessary 
            std::vector<RTN> p_coefs, q_coefs;
            p_coefs = this->coefs;
            q_coefs = q.coefficients(); 
            if (this->deg > q.degree())
                q_coefs = appendZeros<RTN>(q_coefs, this->deg - q.degree());
            else if (this->deg < q.degree())
                p_coefs = appendZeros<RTN>(p_coefs, q.degree() - this->deg);

            // Instantiate the difference polynomial
            std::vector<RTN> diff_coefs;
            for (int i = 0; i < p_coefs.size(); ++i)
                diff_coefs.push_back(p_coefs[i] - q_coefs[i]);

            return Polynomial<N>(diff_coefs);
        }

        template <unsigned M, unsigned P>
        Polynomial<P> operator-(const Polynomial<M>& q) const
        {
            /*
             * Return the result of subtracting by q. 
             */
            typedef number<mpfr_float_backend<P> > RTP;

            // Copy over the coefficients into new vectors with the output
            // precision (P) and resize as necessary 
            std::vector<RTP> p_coefs, q_coefs;
            for (auto&& coef : this->coefs)      p_coefs.push_back(convertPrecision<N, P>(coef));
            for (auto&& coef : q.coefficients()) q_coefs.push_back(convertPrecision<M, P>(coef));
            if (this->deg > q.degree())
                q_coefs = appendZeros<RTP>(q_coefs, this->deg - q.degree());
            else if (this->deg < q.degree())
                p_coefs = appendZeros<RTP>(p_coefs, q.degree() - this->deg);

            // Instantiate the difference polynomial
            std::vector<RTP> diff_coefs;
            for (int i = 0; i < p_coefs.size(); ++i)
                diff_coefs.push_back(p_coefs[i] - q_coefs[i]);

            return Polynomial<P>(diff_coefs);
        }

        Polynomial<N> operator-(const number<mpfr_float_backend<N> > s) const
        {
            /*
             * Return the result of subtracting by the same-precision scalar s.
             */
            std::vector<number<mpfr_float_backend<N> > > diff_coefs;
            diff_coefs.push_back(this->coefs[0] - s);
            for (int i = 1; i < this->coefs.size(); ++i)
                diff_coefs.push_back(this->coefs[i]);

            return Polynomial<N>(diff_coefs);
        }

        template <unsigned M, unsigned P>
        Polynomial<P> operator-(const number<mpfr_float_backend<M> > s) const
        {
            /*
             * Return the result of subtracting by scalar s.
             */
            typedef number<mpfr_float_backend<P> > RTP;

            std::vector<RTP> diff_coefs;
            diff_coefs.push_back(convertPrecision<N, P>(this->coefs[0]) - convertPrecision<M, P>(s));
            for (int i = 1; i < this->coefs.size(); ++i)
                diff_coefs.push_back(convertPrecision<N, P>(this->coefs[i]));

            return Polynomial<P>(diff_coefs);
        }

        Polynomial& operator-=(const Polynomial<N>& q) 
        {
            /*
             * In-place subtraction by q.
             *
             * Note that the precision of the input polynomial must match 
             * that of the current polynomial. 
             */
            typedef number<mpfr_float_backend<N> > RTN;

            // Copy over the coefficients into new vectors and resize 
            // as necessary 
            std::vector<RTN> p_coefs(this->coefs);
            std::vector<RTN> q_coefs(q.coefficients());
            if (this->deg > q.degree())
                q_coefs = appendZeros<RTN>(q_coefs, this->deg - q.degree());
            else if (this->deg < q.degree())
                p_coefs = appendZeros<RTN>(p_coefs, q.degree() - this->deg);

            // Update polynomial coefficients
            for (int i = 0; i < p_coefs.size(); ++i)
                p_coefs[i] -= q_coefs[i];
            this->coefs = removeTrailingZeros<RTN>(p_coefs);
            this->deg = this->coefs.size() - 1;
            return *this;
        }

        Polynomial& operator-=(const number<mpfr_float_backend<N> > s)
        {
            /*
             * In-place subtraction by scalar s.
             *
             * Note that the precision of the input scalar must match that 
             * of the current polynomial. 
             */
            this->coefs[0] -= s;
            return *this;
        }

        Polynomial<N> operator*(const Polynomial<N>& q) const
        {
            /*
             * Return the result of multiplying by the same-precision polynomial q.
             *
             * TODO Replace with a method based on FFT. 
             */
            typedef number<mpfr_float_backend<N> > RTN;

            int p_deg = this->deg;
            int q_deg = q.degree();
            std::vector<RTN> p_coefs, q_coefs, pq_coefs;
            for (auto&& coef : this->coefs)      p_coefs.push_back(coef);
            for (auto&& coef : q.coefficients()) q_coefs.push_back(coef);

            // Compute product coefficients 
            for (int i = 0; i < p_deg + q_deg + 1; ++i) pq_coefs.push_back(0.0);
            for (int i = 0; i <= p_deg; ++i)
            {
                for (int j = 0; j <= q_deg; ++j)
                {
                    pq_coefs[i + j] += p_coefs[i] * q_coefs[j];
                }
            }
            return Polynomial<N>(pq_coefs);
        }

        template <unsigned M, unsigned P>
        Polynomial<P> operator*(const Polynomial<M>& q) const
        {
            /*
             * Return the result of multiplying by q.
             *
             * TODO Replace with a method based on FFT. 
             */
            typedef number<mpfr_float_backend<M> > RTM;
            typedef number<mpfr_float_backend<P> > RTP;

            int p_deg = this->deg;
            int q_deg = q.degree();
            std::vector<RTP> p_coefs, q_coefs, pq_coefs;
            for (auto&& coef : this->coefs)      p_coefs.push_back(convertPrecision<N, P>(coef));
            for (auto&& coef : q.coefficients()) q_coefs.push_back(convertPrecision<M, P>(coef));

            // Compute product coefficients 
            for (int i = 0; i < p_deg + q_deg + 1; ++i) pq_coefs.push_back(0.0);
            for (int i = 0; i <= p_deg; ++i)
            {
                for (int j = 0; j <= q_deg; ++j)
                {
                    pq_coefs[i + j] += p_coefs[i] * q_coefs[j];
                }
            }
            return Polynomial<P>(pq_coefs);
        }

        Polynomial<N> operator*(const number<mpfr_float_backend<N> > s) const
        {
            /*
             * Return the result of multiplying by the same-precision scalar s.
             */
            typedef number<mpfr_float_backend<N> > RTN;

            std::vector<RTN> p_coefs;
            for (auto&& coef : this->coefs) p_coefs.push_back(coef * s);
            return Polynomial<N>(p_coefs);
        }

        template <unsigned M, unsigned P>
        Polynomial<P> operator*(const number<mpfr_float_backend<M> > s) const
        {
            /*
             * Return the result of multiplying by scalar s.
             */
            typedef number<mpfr_float_backend<P> > RTP;

            RTP t = convertPrecision<M, P>(s);
            std::vector<RTP> p_coefs;
            for (auto&& coef : this->coefs) p_coefs.push_back(convertPrecision<N, P>(coef) * t);
            return Polynomial<P>(p_coefs);
        }

        Polynomial& operator*=(const Polynomial<N>& q) 
        {
            /*
             * In-place multiplication by q.
             *
             * Note that the precision of the input polynomial must match 
             * that of the current polynomial.
             *
             * TODO Replace with a method based on FFT. 
             */
            typedef number<mpfr_float_backend<N> > RTN;

            int p_deg = this->deg;
            int q_deg = q.degree();
            std::vector<RTN> q_coefs = q.coefficients();
            std::vector<RTN> pq_coefs;
            for (int i = 0; i < p_deg + q_deg + 1; ++i) pq_coefs.push_back(0.0);
            for (int i = 0; i <= p_deg; ++i)
            {
                for (int j = 0; j <= q_deg; ++j)
                {
                    pq_coefs[i + j] += this->coefs[i] * q_coefs[j];
                }
            }

            // Update polynomial coefficients
            this->coefs = removeTrailingZeros<RTN>(pq_coefs);
            this->deg = this->coefs.size() - 1;
            return *this;
        }

        Polynomial& operator*=(const number<mpfr_float_backend<N> > s)
        {
            /*
             * In-place multiplication by scalar s.
             *
             * Note that the precision of the input scalar must match that 
             * of the current polynomial. 
             */
            for (int i = 0; i < this->coefs.size(); ++i) this->coefs[i] *= s;
            return *this;
        }

        Polynomial<N> operator/(const number<mpfr_float_backend<N> > s) const
        {
            /*
             * Return the result of dividing by the same-precision scalar s.
             */
            typedef number<mpfr_float_backend<N> > RTN;

            std::vector<RTN> p_coefs;
            for (auto&& coef : this->coefs) p_coefs.push_back(coef / s);
            return Polynomial<N>(p_coefs);
        }

        template <unsigned M, unsigned P>
        Polynomial<P> operator/(const number<mpfr_float_backend<M> > s) const
        {
            /*
             * Return the result of dividing by scalar s.
             */
            typedef number<mpfr_float_backend<P> > RTP;

            RTP t = convertPrecision<M, P>(s); 
            std::vector<RTP> p_coefs;
            for (auto&& coef : this->coefs) p_coefs.push_back(convertPrecision<N, P>(coef) / t);
            return Polynomial<P>(p_coefs);
        }

        Polynomial& operator/=(const number<mpfr_float_backend<N> > s)
        {
            /*
             * In-place division by scalar s.
             *
             * Note that the precision of the input scalar must match that 
             * of the current polynomial. 
             */
            for (auto&& coef : this->coefs) coef /= s;
            return *this;
        }

        Polynomial<N> operator-() const
        {
            /*
             * Return the negative of the current polynomial. 
             */
            std::vector<number<mpfr_float_backend<N> > > coefs(this->coefs); 
            for (auto&& coef : coefs) coef *= (-1);
            return Polynomial<N>(coefs);
        }

        number<mpfr_float_backend<N> > leadingCoef() const
        {
            /*
             * Return the leading coefficient of the current polynomial.
             */
            return this->coefs[this->deg];
        }

        Polynomial<N> reduce() const
        {
            /*
             * Factor out all powers of x from the current polynomial. 
             */
            typedef number<mpfr_float_backend<N> > RTN;

            // Skip over all low-degree zero coefficients 
            std::vector<RTN> reduced;
            int i = 0;
            while (i < this->coefs.size() && this->coefs[i] == 0.0) i++;
            for (int j = i; j < this->coefs.size(); ++j)
                reduced.push_back(this->coefs[j]);

            return Polynomial<N>(reduced);
        }

        Polynomial<N> monic() const
        {
            /*
             * Factor out the leading coefficient. 
             */
            typedef number<mpfr_float_backend<N> > RTN;

            std::vector<RTN> divided;
            RTN leading_coef = this->leadingCoef();
            for (auto&& coef : this->coefs)
                divided.push_back(coef / leading_coef);

            return Polynomial<N>(divided);
        }

        std::pair<std::vector<number<mpc_complex_backend<N> > >, bool> roots(SolveMethod method = Aberth,
                                                                             int max_iter = 1000,
                                                                             double atol = 1e-15,
                                                                             double rtol = 1e-15)
        {
            /*
             * Return all complex roots of the polynomial. 
             */
            typedef number<mpfr_float_backend<N> > RTN;

            if (this->deg == 2)
            {
                return std::make_pair(quadratic<N>(this->coefs[1] / this->coefs[2], this->coefs[0] / this->coefs[2]), true);
            }
            else if (this->deg == 3)
            {
                RTN a = this->coefs[3];
                return std::make_pair(cubic<N>(this->coefs[2] / a, this->coefs[1] / a, this->coefs[0] / a), true);
            }
            else if (method == Aberth)
            {
                return this->rootsAberth<N>(max_iter, atol, rtol);
            }
            else
            {
                return this->rootsWeierstrass<N>(max_iter, atol, rtol);
            }
        }

        template <unsigned M>
        std::pair<std::vector<number<mpc_complex_backend<M> > >, bool> roots(SolveMethod method = Aberth,
                                                                             int max_iter = 1000,
                                                                             double atol = 1e-15,
                                                                             double rtol = 1e-15)
        {
            /*
             * Return all complex roots of the polynomial.  
             */
            typedef number<mpfr_float_backend<M> > RTM;

            if (this->deg == 2)
            {
                RTM a = convertPrecision<N, M>(this->coefs[2]);
                RTM b_over_a = convertPrecision<N, M>(this->coefs[1]) / a;
                RTM c_over_a = convertPrecision<N, M>(this->coefs[0]) / a;
                return std::make_pair(quadratic<M>(b_over_a, c_over_a), true);
            }
            else if (this->deg == 3)
            {
                RTM a = convertPrecision<N, M>(this->coefs[3]);
                RTM b_over_a = convertPrecision<N, M>(this->coefs[2]) / a;
                RTM c_over_a = convertPrecision<N, M>(this->coefs[1]) / a;
                RTM d_over_a = convertPrecision<N, M>(this->coefs[0]) / a;
                return std::make_pair(cubic<M>(b_over_a, c_over_a, d_over_a), true);
            }
            else if (method == Aberth)
            {
                return this->rootsAberth<M>(max_iter, atol, rtol);
            }
            else
            {
                return this->rootsWeierstrass<M>(max_iter, atol, rtol);
            }
        }
};

template <unsigned N>
Polynomial<N> operator+(const number<mpfr_float_backend<N> >& s, Polynomial<N> rhs)
{
    /*
     * Left-hand addition by the same-precision scalar s.
     */
    return rhs + s; 
}

template <unsigned N>
Polynomial<N> operator-(const number<mpfr_float_backend<N> >& s, Polynomial<N> rhs)
{
    /*
     * Left-hand subtraction by the same-precision scalar s. 
     */
    return (-rhs) + s;  
}

template <unsigned N>
Polynomial<N> operator*(const number<mpfr_float_backend<N> >& s, Polynomial<N> rhs)
{
    /*
     * Left-hand multiplication by the same-precision scalar s. 
     */
    return rhs * s;
}

#endif
