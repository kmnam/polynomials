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

        T real()
        {
            /*
             * Return the real part of this complex number.
             */
            return this->a;
        }

        T imag()
        {
            /*
             * Return the imaginary part of this complex number.
             */
            return this->b;
        }

        ComplexNumber conj()
        {
            /*
             * Return the conjugate of this complex number.
             */
            return ComplexNumber(this->a, -this->b);
        }

        T abs2()
        {
            /*
             * Return the squared modulus of this complex number.
             */
            return (this->a * this->a) + (this->b * this->b);
        }

        ComplexNumber operator+(const ComplexNumber<T>& z)
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

        ComplexNumber operator-(const ComplexNumber<T>& z)
        {
            /*
             * Return the result of subtracting by z.
             */
            return ComplexNumber(this->a - z.real(), this->b - z.imag());
        }

        ComplexNumber operator*(const ComplexNumber<T>& z)
        {
            /*
             * Return the result of multiplying by z.
             */
            return ComplexNumber(
                this->a * z.real() - this->b * z.imag(),
                this->a * z.imag() + this->b * z.real()
            );
        }

        ComplexNumber operator/(const ComplexNumber<T>& z)
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
};

#endif
