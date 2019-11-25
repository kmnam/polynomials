#ifndef UNIVAR_POLYNOMIAL_HPP
#define UNIVAR_POLYNOMIAL_HPP

#include <utility>
#include <ostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <complex>
#include <limits>
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
 *     11/24/2019
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
            w = w.head(size);
        }
        else break;
    }
    return w;
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
    // Define a function for evaluating the given polynomial with 
    // Horner's method
    std::function<CT(CT)> eval = [coefs](CT x)
    {
        CT y = coefs[coefs.size() - 1];
        for (int i = coefs.size() - 2; i >= 0; --i)
        {
            y *= x;
            y += coefs[i];
        }
        return y;
    };

    // Perform the Weierstrass correction for each root
    std::vector<CT> new_roots(roots);
    for (unsigned j = 0; j < roots.size(); ++j)
    {
        CT denom = 1.0;
        for (unsigned k = 0; k < roots.size(); ++k)
        {
            if (j != k) denom *= (new_roots[j] - new_roots[k]);
        }
        CT correction = -eval(new_roots[j]) / denom;
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


template <typename T>
class Polynomial
{
    private:
        unsigned deg;                    // Degree of the polynomial
        unsigned prec;                   // Maximum number of digits that can
                                         // be represented by type T 
        Matrix<T, Dynamic, 1> coefs;     // Coefficients of the polynomial
                                         // stored in ascending order of power

        Polynomial multiplyFFT(const Polynomial<T>& q) const
        {
            /*
             * Return the result of multiplying by q.
             *
             * Compute the product with q using the FFT algorithm:
             * 1) Evaluate the two polynomials at the (2d + 2)-th roots of 
             *    unity, where d = max(this->deg, q.degree).
             * 2) Multiply the pairs of values obtained above.
             * 3) Interpolate the pairs of values to obtain the product
             *    polynomial. 
             */
            const double two_pi = 2.0 * std::acos(-1);
            unsigned p_deg = this->deg;
            unsigned q_deg = q.degree();
            Matrix<T, Dynamic, 1> p_coefs;
            Matrix<T, Dynamic, 1> q_coefs;
            if (p_deg > q_deg)
            {
                p_coefs = this->coefs;
                q_coefs = appendZeros<T>(q.coefficients(), p_deg - q_deg);
            }
            else if (p_deg < q_deg)
            {
                q_coefs = q.coefficients();
                p_coefs = appendZeros<T>(this->coefs(), q_deg - p_deg);
            }

            // Evaluate the two polynomials at the (2d)-th roots of unity
            unsigned pq_deg = 2 * std::max(p_deg, q_deg);
            Array<std::complex<T>, Dynamic, 1> values_p(pq_deg);
            Array<std::complex<T>, Dynamic, 1> values_q(pq_deg);
            for (unsigned i = 0; i < pq_deg; ++i)
            {
                T a = std::cos(two_pi * i / pq_deg);
                T b = std::sin(two_pi * i / pq_deg);
                std::complex<T> z(a, b);
                std::cout << i << " " << z << std::endl;
                std::cout << this->eval(z) << std::endl;
                std::cout << q.eval(z) << std::endl;
                values_p(i) = this->eval(z);
                values_q(i) = q.eval(z);
            }
            std::cout << values_p << "\n\n";
            std::cout << values_q << "\n\n";

            // Multiply the two vectors of values pointwise
            auto values_pq = (values_p * values_q).matrix();

            // Interpolate the product values by the inverse DFT
            Matrix<std::complex<T>, Dynamic, Dynamic> inv_dft(pq_deg, pq_deg);
            for (unsigned i = 0; i < pq_deg; ++i)
            {
                inv_dft(0, i) = 1.0;
                inv_dft(i, 0) = 1.0;
            }
            for (unsigned i = 1; i < pq_deg; ++i)
            {
                for (unsigned j = 1; j < pq_deg; ++j)
                {
                    T a = std::cos(two_pi * i * j / pq_deg);
                    T b = -std::sin(two_pi * i * j / pq_deg);
                    std::complex<T> z(a, b);
                    inv_dft(i, j) = z; 
                }
            }
            std::cout << values_pq << "\n\n";
            std::cout << inv_dft << "\n\n";
            std::cout << inv_dft * values_pq << "\n\n";
            Matrix<T, Dynamic, 1> prod_coefs = (inv_dft * values_pq).array().real().matrix();
            return Polynomial(prod_coefs);
        }

        Matrix<std::complex<T>, Dynamic, 1> rootsWeierstrass(unsigned max_iter = 20)
        {
            /*
             * Run Weierstrass' method on the given polynomial, returning the 
             * roots as type std::complex<T>.
             */
            bool converged = false;
            bool quadratic = false;
            unsigned iter = 0;

            // Start with std::complex<T>
            std::vector<std::complex<T> > roots;
            for (unsigned i = 0; i < this->deg; ++i)
            {   // Initialize the roots to (0.4 + 0.9 i)^p, for p = 0, ..., d - 1
                roots.push_back(std::pow(std::complex<T>(0.4, 0.9), i)); 
            }

            // Write coefficients to a std::complex<T> vector
            std::vector<std::complex<T> > coefs;
            for (unsigned i = 0; i < this->coefs.size(); ++i)
            {
                std::complex<T> z(this->coefs(i), 0.0);
                coefs.push_back(z);
            }

            // Run the algorithm until all roots exhibit quadratic convergence
            // or the desired number of iterations is reached
            std::vector<double> delta, new_delta;
            for (unsigned i = 0; i < this->deg; ++i) delta.push_back(0.0);
            while (iter < max_iter && !converged)
            {
                auto result = weierstrass<std::complex<T> >(coefs, roots);
                std::vector<std::complex<T> > new_roots = result.first;
                new_delta = result.second;
                iter++;
   
                // Check that the change is less than the square of the 
                // previous change for each root
                for (unsigned j = 0; j < this->deg; ++j)
                {
                    quadratic = true;
                    if (new_delta[j] >= 0.01 * (delta[j] * delta[j]))
                    {
                        quadratic = false;
                        break;
                    }
                }
                if (quadratic) converged = true;
                roots = new_roots;
                delta = new_delta;
            }

            // If not converged, try again with a boost::multiprecision::mpc type
            unsigned prec = std::numeric_limits<T>::max_digits10;
            while (!converged)
            {
                using boost::multiprecision::pow;
                using boost::multiprecision::abs;
                if (prec < 30)
                {
                    // Re-do the same calculations with mpc_30
                    prec = 30; iter = 0;
                    std::vector<mpc_30> roots2, coefs2;
                    for (unsigned i = 0; i < this->deg; ++i)
                    {
                        roots2.push_back(pow(mpc_30("0.4", "0.9"), i));
                        delta[i] = 0.0;
                    }
                    for (unsigned i = 0; i < this->coefs.size(); ++i)
                    {
                        std::stringstream ss;
                        ss << std::numeric_limits<T>::max_digits10 << this->coefs(i);
                        mpc_30 z(ss.str(), "0.0");
                        coefs2.push_back(z);
                    }
                    while (iter < max_iter && !converged)
                    {
                        auto result = weierstrass<mpc_30>(coefs2, roots2);
                        std::vector<mpc_30> new_roots = result.first;
                        new_delta = result.second;
                        iter++;
                        for (unsigned j = 0; j < this->deg; ++j)
                        {
                            quadratic = true;
                            if (new_delta[j] >= 0.01 * (delta[j] * delta[j]))
                            {
                                quadratic = false;
                                break;
                            }
                        }
                        if (quadratic) converged = true;
                        roots2 = new_roots;
                        delta = new_delta;
                    }
                }
                else if (prec < 60)
                {
                    // Re-do the same calculations with mpc_60
                    prec = 60; iter = 0;
                    std::vector<mpc_60> roots2, coefs2;
                    for (unsigned i = 0; i < this->deg; ++i)
                    {
                        roots2.push_back(pow(mpc_60("0.4", "0.9"), i));
                        delta[i] = 0.0;
                    }
                    for (unsigned i = 0; i < this->coefs.size(); ++i)
                    {
                        std::stringstream ss;
                        ss << std::numeric_limits<T>::max_digits10 << this->coefs(i);
                        mpc_60 z(ss.str(), 0.0);
                        coefs2.push_back(z);
                    }
                    while (iter < max_iter && !converged)
                    {
                        auto result = weierstrass<mpc_60>(coefs2, roots2);
                        std::vector<mpc_60> new_roots = result.first;
                        new_delta = result.second;
                        iter++;
                        for (unsigned j = 0; j < this->deg; ++j)
                        {
                            quadratic = true;
                            if (new_delta[j] >= 0.01 * (delta[j] * delta[j]))
                            {
                                quadratic = false;
                                break;
                            }
                        }
                        if (quadratic) converged = true;
                        roots2 = new_roots;
                        delta = new_delta;
                    }
                }
                else if (prec < 100)
                {
                    // Re-do the same calculations with mpc_100
                    prec = 100; iter = 0;
                    std::vector<mpc_100> roots2, coefs2;
                    for (unsigned i = 0; i < this->deg; ++i)
                    {
                        roots2.push_back(pow(mpc_100("0.4", "0.9"), i));
                        delta[i] = 0.0;
                    }
                    for (unsigned i = 0; i < this->coefs.size(); ++i)
                    {
                        std::stringstream ss;
                        ss << std::numeric_limits<T>::max_digits10 << this->coefs(i);
                        mpc_100 z(ss.str(), 0.0);
                        coefs2.push_back(z);
                    }
                    while (iter < max_iter && !converged)
                    {
                        auto result = weierstrass<mpc_100>(coefs2, roots2);
                        std::vector<mpc_100> new_roots = result.first;
                        new_delta = result.second;
                        iter++;
                        for (unsigned j = 0; j < this->deg; ++j)
                        {
                            quadratic = true;
                            if (new_delta[j] >= 0.01 * (delta[j] * delta[j]))
                            {
                                quadratic = false;
                                break;
                            }
                        }
                        if (quadratic) converged = true;
                        roots2 = new_roots;
                        delta = new_delta;
                    }
                }
                else if (prec < 200)
                {
                    // Re-do the same calculations with mpc_200
                    prec = 200; iter = 0;
                    std::vector<mpc_200> roots2, coefs2;
                    for (unsigned i = 0; i < this->deg; ++i)
                    {
                        roots2.push_back(pow(mpc_200("0.4", "0.9"), i));
                        delta[i] = 0.0;
                    }
                    for (unsigned i = 0; i < this->coefs.size(); ++i)
                    {
                        std::stringstream ss;
                        ss << std::numeric_limits<T>::max_digits10 << this->coefs(i);
                        mpc_200 z(ss.str(), 0.0);
                        coefs2.push_back(z);
                    }
                    while (iter < max_iter && !converged)
                    {
                        auto result = weierstrass<mpc_200>(coefs2, roots2);
                        std::vector<mpc_200> new_roots = result.first;
                        new_delta = result.second;
                        iter++;
                        for (unsigned j = 0; j < this->deg; ++j)
                        {
                            quadratic = true;
                            if (new_delta[j] >= 0.01 * (delta[j] * delta[j]))
                            {
                                quadratic = false;
                                break;
                            }
                        }
                        if (quadratic) converged = true;
                        roots2 = new_roots;
                        delta = new_delta;
                    }
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
                final_roots(i) = roots[i];
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

        Matrix<std::complex<T>, Dynamic, 1> roots()
        {
            /*
             * Return all complex roots of the polynomial.  
             */
            return this->rootsWeierstrass();
        }

        Matrix<T, Dynamic, 1> positiveRoots(double imag_tol = 1e-10)
        {
            /*
             * Return all positive roots with imaginary part less than the
             * given tolerance.
             */
            Matrix<std::complex<T>, Dynamic, 1> roots = this->rootsWeierstrass();
            Matrix<T, Dynamic, 1> pos_roots;
            unsigned i = 0;
            for (unsigned j = 0; j < roots.size(); ++j)
            {
                if (roots(j).real() > 0 && std::abs(roots(j).imag()) < imag_tol)
                {
                    i++;
                    pos_roots.conservativeResize(i);
                    pos_roots(i-1) = roots(j).real();
                }
            }
            return pos_roots;
        }
};

#endif