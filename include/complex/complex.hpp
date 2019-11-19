#ifndef COMPLEX_NUMBER_HPP
#define COMPLEX_NUMBER_HPP
#include <cmath>

/*
 * Class template implementation of a complex number type with arbitrary
 * real and imaginary part types.
 *
 * Authors:
 *     Kee-Myoung Nam, Department of Systems Biology, Harvard Medical School
 * Last updated:
 *     11/19/2019
 */

template <typename T>
class ComplexNumber
{
    private:
        T a;    // Real part
        T b;    // Imaginary part

    public:
        ComplexNumber()
        {
            /*
             * Empty constructor; initialize to zero.
             */
            this->a = 0.0;
            this->b = 0.0;
        }

        ComplexNumber(T a)
        {
            /*
             * Trivial constructor with specified real part.
             */
            this->a = a;
            this->b = 0.0;
        }

        ComplexNumber(T a, T b)
        {
            /*
             * Trivial constructor with specified values.
             */
            this->a = a;
            this->b = b;
        }

        ~ComplexNumber()
        {
            /*
             * Empty destructor.
             */
        }

        T real() const
        {
            /*
             * Return the real part of this complex number.
             */
            return this->a;
        }

        T imag() const
        {
            /*
             * Return the imaginary part of this complex number.
             */
            return this->b;
        }

        ComplexNumber& operator=(const T x) 
        {
            /*
             * Assignment from real type.
             */
            this->a = x;
            this->b = 0.0;
            return *this;
        }

        ComplexNumber conj() const
        {
            /*
             * Return the conjugate of this complex number.
             */
            return ComplexNumber(this->a, -this->b);
        }

        T abs2() const
        {
            /*
             * Return the squared modulus of this complex number.
             */
            return (this->a * this->a) + (this->b * this->b);
        }

        ComplexNumber operator+(const ComplexNumber<T>& z) const
        {
            /*
             * Return the result of adding by z.
             */
            return ComplexNumber(this->a + z.real(), this->b + z.imag());
        }

        ComplexNumber operator+(const T x) const
        {
            /*
             * Return the result of adding by x.
             */
            return ComplexNumber(this->a + x, this->b);
        }

        ComplexNumber& operator+=(const ComplexNumber<T>& z)
        {
            /*
             * In-place addition by z.
             */
            this->a += z.real();
            this->b += z.imag();
            return *this;
        }

        ComplexNumber& operator+=(const T x)
        {
            /*
             * In-place addition by x.
             */
            this->a += x;
            return *this;
        }

        ComplexNumber operator-(const ComplexNumber<T>& z) const
        {
            /*
             * Return the result of subtracting by z.
             */
            return ComplexNumber(this->a - z.real(), this->b - z.imag());
        }

        ComplexNumber operator-(const T x) const
        {
            /*
             * Return the result of subtracting by x.
             */
            return ComplexNumber(this->a - x, this->b);
        }

        ComplexNumber& operator-=(const ComplexNumber<T>& z)
        {
            /*
             * In-place subtraction by z.
             */
            this->a -= z.real();
            this->b -= z.imag();
            return *this;
        }

        ComplexNumber& operator-=(const T x)
        {
            /*
             * In-place subtraction by x.
             */
            this->a -= x;
            return *this;
        }

        ComplexNumber operator*(const ComplexNumber<T>& z) const
        {
            /*
             * Return the result of multiplying by z.
             */
            return ComplexNumber(
                this->a * z.real() - this->b * z.imag(),
                this->a * z.imag() + this->b * z.real()
            );
        }

        ComplexNumber operator*(const T x) const
        {
            /*
             * Return the result of multiplying by x.
             */
            return ComplexNumber(this->a * x, this->b * x);
        }

        ComplexNumber& operator*=(const ComplexNumber<T>& z)
        {
            /*
             * In-place multiplication by z.
             */
            T a_ = this->a * z.real() - this->b * z.imag();
            T b_ = this->a * z.imag() + this->b * z.real();
            this->a = a_;
            this->b = b_;
            return *this;
        }

        ComplexNumber& operator*=(const T x)
        {
            /*
             * In-place multiplication by x.
             */
            this->a *= x;
            this->b *= x;
            return *this;
        }

        ComplexNumber operator/(const ComplexNumber<T>& z) const
        {
            /*
             * Return the result of dividing by z.
             */
            T modulus = z.abs2();
            return ComplexNumber(
                (this->a * z.real() + this->b * z.imag()) / modulus,
                (this->b * z.real() - this->a * z.imag()) / modulus
            );
        }

        ComplexNumber operator/(const T x) const
        {
            /*
             * Return the result of dividing by x.
             */
            return ComplexNumber(this->a / x, this->b / x);
        }

        ComplexNumber& operator/=(const ComplexNumber<T>& z)
        {
            /*
             * In-place division by z.
             */
            T modulus = z.abs2();
            T a_ = (this->a * z.real() + this->b * z.imag()) / modulus;
            T b_ = (this->b * z.real() - this->a * z.imag()) / modulus;
            this->a = a_;
            this->b = b_;
            return *this;
        }

        ComplexNumber& operator/=(const T x) const
        {
            /*
             * In-place division by x.
             */
            this->a /= x;
            this->b /= x;
            return *this;
        }
};

template <typename T>
ComplexNumber<T> operator+(const T x, const ComplexNumber<T>& z)
{
    /*
     * Return the result of adding real x and complex z.
     */
    return ComplexNumber<T>(x + z.real(), z.imag());
}

template <typename T>
ComplexNumber<T> operator-(const T x, const ComplexNumber<T>& z)
{
    /*
     * Return the result of subtracting complex z from real x.
     */
    return ComplexNumber<T>(x - z.real(), -z.imag());
}

template <typename T>
ComplexNumber<T> operator*(const T x, const ComplexNumber<T>& z)
{
    /*
     * Return the result of multiplying real x and complex z.
     */
    return ComplexNumber<T>(x * z.real(), x * z.imag());
}

template <typename T>
ComplexNumber<T> operator/(const T x, const ComplexNumber<T>& z)
{
    /*
     * Return the result of dividing real x by complex z.
     */
    T modulus = z.abs2();
    return ComplexNumber<T>(x * z.real() / modulus, -x * z.imag() / modulus);
}

#endif
