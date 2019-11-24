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
#include "complex/complex.hpp"

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
            if (this->degree == 0)
                return Polynomial();

            Matrix<T, Dynamic, 1> dcoefs = (
                Array<T, Dynamic, 1>::LinSpaced(this->degree, 1, this->degree)
                * this->coefs.tail(this->degree).array()
            ).matrix();
            return Polynomial(dcoefs);
        }
};

#endif
