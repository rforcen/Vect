//
//  DomainColoring.metal
//  MetalKernelSum
//
//  Created by asd on 01/04/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//

#include <metal_stdlib>
using namespace metal;
#include "complex.h"

uint HSV2int(float h, const float s, const float v) { // convert hsv to int with alpha 0xff00000
    float r = 0, g = 0, b = 0;
    
    if (s == 0) r = g = b = v;
    else {
        if (h == 1)  h = 0;
        
        float z = floor(h * 6.);
        int i = z;
        
        float   f = h * 6 - z,
        p = v * (1 - s), q = v * (1 - s * f),
        t = v * (1 - s * (1 - f));
        
        switch (i) {
            case 0:          r = v;   g = t;   b = p;             break;
            case 1:          r = q;   g = v;   b = p;             break;
            case 2:          r = p;   g = v;   b = t;             break;
            case 3:          r = p;   g = q;   b = v;             break;
            case 4:          r = t;   g = p;   b = v;             break;
            case 5:          r = v;   g = p;   b = q;             break;
        }
    }
    uint c, color = 0xff000000;
    // alpha = 0xff
    c = (uint)(256. * r) & 0xff;
    color |= c;
    c = (uint)(256. * g) & 0xff;
    color |= (c << 8);
    c = (uint)(256. * b) & 0xff;
    color |= (c << 16);
    return color;
}

// main function evaluator
inline Complex<float> func(const Complex<float>z)  {
    return (z*z+1.)/(z*z-1.); //c*c+c+1.4;
}

inline float pow3(float x) { return x*x*x; }

kernel void generateDomainColor(
                                device uint*colors[[buffer(0)]], // buffer per device
                                const device int &dev[[buffer(1)]], // current device
                                const device int &nDevices[[buffer(2)]],
                                
                                const device uint &w[[buffer(3)]],
                                const device uint &h[[buffer(4)]],
                                
                                const uint tPerTg [[ threads_per_threadgroup ]],
                                const uint t[[thread_position_in_grid]]
                                ) {
    
    // device
    uint szdev=h/nDevices; // chop 'h' in nDevices
    uint fd=dev*szdev, td=(dev==nDevices-1) ? h : (dev+1)*szdev;
    
    uint nThreads=tPerTg*tPerTg; // thread range
    uint segsizeh=h/(nThreads*nDevices), fromh=fd + t * segsizeh, toh = (t==nThreads-1) ? td : fd + (t+1)*segsizeh;
    uint icol=t * segsizeh * w; // inital color index

    /*
    uint nThreads=tPerTg*tPerTg; // thread range
    uint segsizeh=h/nThreads, fromh=t * segsizeh, toh=(t==nThreads-1) ? h : (t+1)*segsizeh;
    uint icol=fromh * w; // inital color idnex
     */

    float E = 2.7182818284590452353602874713527;
    float M_PI = 3.141592653589793238462643383;
    
    float PI = M_PI, PI2 = PI * 2.;
    float limit=PI,  rmi = -limit, rma = limit, imi = -limit, ima = limit;
    
    for (uint j = fromh; j < toh; j++) {
        float im = ima - (ima - imi) * j / (h - 1);
        
        for (uint i = 0; i < w; i++) {
            float re = rma - (rma - rmi) * i / (w - 1);
            
            auto v = func(Complex<float>(re, im)); // evaluate
            
            float ang = v.arg(); // -> hsv
            while (ang < 0) ang += PI2;
            ang /= PI2;
            
            float m = v.abs(), ranges = 0., rangee = 1.;
            while (m > rangee) {
                ranges = rangee;
                rangee *= E;
            }
            
            float k = (m - ranges) / (rangee - ranges);
            float kk = (k < 0.5 ? k * 2. : 1. - (k - 0.5) * 2);
            
            float sat = 0.4 + (1. - pow3(1. - (kk)))     * 0.6;
            float val = 0.6 + (1. - pow3(1. - (1 - kk))) * 0.4;
            
            colors[icol++] = HSV2int(ang, sat, val);
        }
    }
}


