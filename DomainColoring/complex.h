//
//  complex.h
//  Vect
//
//  Created by asd on 02/04/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//

#ifndef complex_h
#define complex_h

template<typename T>
struct Complex
{
    T re=0, im=0;
    
    inline Complex() { re = im = 0; }
    inline Complex(T re, T im) : re(re), im(im) { }
    
    inline T arg() const { return atan2(im, re);  }
    inline T abs() const { return sqrt(sqmag());  }
    inline T sqmag() const {  return re*re + im*im;  }
    
    inline void operator=(const thread Complex<T>&other) {
        re=other.re; im=other.im;
    }
    
    inline Complex<T> operator*(const thread Complex<T>& other) const {
        return Complex(re*other.re - im*other.im,
                       re*other.im + im*other.re);
    }
    inline Complex<T> operator/(const thread Complex<T>&other) const {
        T div=(other.re*other.re) + (other.im*other.im);
        Complex<T>tmp;
        tmp.re=(re*other.re)+(im*other.im);
        tmp.re/=div;
        tmp.im=(im*other.re)-(re*other.im);
        tmp.im/=div;
        return tmp;
    }
    inline Complex<T> operator+(const thread Complex<T>& other)  const {
        return Complex(re + other.re, im + other.im);
    }
    
    inline Complex<T> operator-(const thread  Complex<T>& other)  const {
        return Complex(re - other.re, im - other.im);
    }
    
    inline Complex<T> operator*(const thread T& c) const  {   return Complex(re * c, im * c); }
    inline Complex<T> operator+(const thread  T& c) const {   return Complex(re + c, im);    }
    inline Complex<T> operator-(const thread  T& c) const {   return Complex(re - c, im);    }
};

#endif /* complex_h */
