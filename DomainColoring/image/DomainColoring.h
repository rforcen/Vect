
#ifndef CDomainColoringH
#define CDomainColoringH

#include <complex>
#include "zCompiler.h"

using namespace std;

//typedef float realNum;
#define realNum float

class DomainColoring {
    typedef complex<realNum>Complex;
    typedef Complex(DomainColoring::*MemberFuncPtr)(Complex);
    
    Complex c, v, ci;
    int totCol, icol;
    zCompiler zComp;
    
public:
    
    DomainColoring() {}
    
    bool compileFunc(string expr) {
        zComp.Compile(expr);
        return zComp.err;
    }
    inline realNum pow3(realNum x) { return x*x*x; }
    
    uint32_t*getColor(uint32_t*colors, int w, int h, realNum rmi, realNum rma, realNum imi, realNum ima) {
        realNum PI = M_PI, PI2 = PI * 2;
        realNum E = 2.7182818284590452353602874713527;
        //        realNum limit = PI;
        //        realNum rmi, rma, imi, ima;
        //        rmi = -limit; rma = limit; imi = -limit; ima = limit;
        icol = 0;
        
        try {
            
            for (int j = 0; j < h; j++) {
                realNum im = ima - (ima - imi) * j / (h - 1);
                
                for (int i = 0; i < w; i++) {
                    realNum re = rma - (rma - rmi) * i / (w - 1);
                    
                    v = Complex(re, im); //zComp.Execute(Complex(re, im)); // fun(c); // evaluate here
                    
                    realNum a = arg(v);
                    while (a < 0) a += PI2;
                    a /= PI2;
                    realNum m = abs(v), ranges = 0, rangee = 1;
                    while (m > rangee) {
                        ranges = rangee;
                        rangee *= E;
                    }
                    
                    realNum k = (m - ranges) / (rangee - ranges);
                    realNum kk = (k < 0.5 ? k * 2 : 1 - (k - 0.5) * 2);
                    
                    realNum sat = 0.4 + (1 - pow3(1 - (kk)))     * 0.6;
                    realNum val = 0.6 + (1 - pow3(1 - (1 - kk))) * 0.4;
                    
                    colors[icol++] = HSV2int(a, sat, val);
                }
            }
        }
        catch (...) {
        }
        
        return colors;
    }
    
    int HSV2int(realNum h, realNum s, realNum v) { // convert hsv to int with alpha 0xff00000
        realNum r = 0, g = 0, b = 0;
        if (s == 0)
            r = g = b = v;
        else {
            if (h == 1)
                h = 0;
            realNum z = floor(h * 6);
            int i = (int)(z);
            realNum f = h * 6 - z, p = v * (1 - s), q = v * (1 - s * f), t =
            v * (1 - s * (1 - f));
            
            switch (i) {
                case 0:
                    r = v;
                    g = t;
                    b = p;
                    break;
                case 1:
                    r = q;
                    g = v;
                    b = p;
                    break;
                case 2:
                    r = p;
                    g = v;
                    b = t;
                    break;
                case 3:
                    r = p;
                    g = q;
                    b = v;
                    break;
                case 4:
                    r = t;
                    g = p;
                    b = v;
                    break;
                case 5:
                    r = v;
                    g = p;
                    b = q;
                    break;
            }
        }
        int c, color = 0xff000000;
        // alpha = 0xff
        c = (int)(256 * r) & 0xff;
        color |= c;
        c = (int)(256 * g) & 0xff;
        color |= (c << 8);
        c = (int)(256 * b) & 0xff;
        color |= (c << 16);
        return color;
    }
    
    // randon formula generator
    int rnd(int r) { return rand() % r;    }
    string randomNumber() { return std::to_string((float)rand()/RAND_MAX); }
    string randomConst() {
        switch (rnd(2)) {
            case 0: return "C( " + randomNumber() + ", " + randomNumber() + " )";
            case 1: return randomNumber();
        }
        return "";
    }
    string randomFunc() {      return zComp.toLower(zComp.funcNames()[rnd(12)]);    }
    string randomOperator() {  return " " + string(1, "+-/*^"[rnd(5)]) + " "; }
    string randomFormula(int level=0, int maxLevel=2) {
        if(level<maxLevel)
            for (int i=0; i<rnd(3)+1 ; i++)
                return randomFormula(level+1) + randomOperator() + randomFormula(level+1);
        switch (rnd(4)) {
            case 0: return randomConst();
            case 1: return "z";
            case 2: return randomFunc() + "( " + randomFormula(level+1) + " )";
            case 3: return "( " + randomFormula(level+1) + " )";
        }
        return "";
    }
};


extern DomainColoring *dcCommom;
#endif
