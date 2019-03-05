//
//  Vect.cpp
//  formants
//
//  Created by asd on 03/03/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//

#include "Vect.hpp"
#include <math.h>

void testVect() {
    typedef Vect<float> VF;
    
    VF vinit(5, new float[5]{0,1,2,3,4});
    VF v0, v1(100), v2(100);
    auto v4=v1;
    
    int n=100000; // append test
    float rs=0;
    for (int i=0; i<n; i++) { v0 << i; rs+=i; }
    assert(v0.sum()==rs && "append failure");
    
    // test prefix ++, --
    auto vorg=v0;
    ++v0;
    --v0;
    assert(vorg==v0 && "inc, dec failure");
    
    // seq & Vect[Vect++] vect indexing post increment
    auto vseq=VF::seq(0, 10, 0.1), vs0=vseq;
    for (auto d:vseq) vs0[vseq++]=d+1; // can access & increment iterator index
    assert(vseq.sum() == vs0.sum()-vs0.count());
    
    // aritmethic
    v1.random(); v2.random(); v4.random();
    v0 = (v1+v2/v4)*3;
    
    v0.clear(); // vect append
    
    auto v01=VF::rnd(400), v02=VF::rnd(500);
    v0 << v01;
    v0 << v02;
    
    // for bigger v01, v02 sizes precission float may just not enough
    assert(fabs(v0.sum() - (v01.sum()+v02.sum())) < 1e-3);
    
    for (int i=0; i<v1.length(); i++) // index (v2=v1)
        v2[i] = v1[i];
    // logical ops
    assert(v1==v1     && "== failure");
    assert(v1==v2     && "!= failure");
    assert(!(v1>v2)   && ">  failure");
    
    float sum=0;
    for (auto vi:v1) sum+=vi; // iterator
    assert(sum==v1.sum() && "iterator failure");
    
    // reverse
    auto vrev = VF::seq(10), vw=vrev;
    assert(vrev == vw.reverse().reverse() && "reverse failure");
    
    // subvector (from, to)
    auto vsub = v0.sub(10,20);
    assert(vsub.length()==10 && "sub failure");
    
    // filter
    auto vf = v0.filter([](float x) -> bool { return x>0.5 && x<0.7; }).sort();
    for (auto d:vf) assert((d>0.5 && d<0.7) && "filter failure");
    
    vf=0; // set all values to zero & check scalar
    assert(vf==0    && "== scalar failure");
    assert(vf!=10   && "== scalar failure");
    assert(vf>-10   && "> scalar failure");
    assert(vf<10    && "< scalar failure");
    
    // generate complex wave
    auto vwave = VF::genWave(1, 8000, 3,
                             new float[3]{1,0.4,0.6}, new float[3]{200, 600, 1800}, new float[3]{0,0,0});
    
    VF v3(v0); // misc
    auto v6=v3.scale(0,1);
    v6.sequence(0,5);
    
    // lambda apply func / sort
    v2.apply(sinf).sort().apply([](float x) -> float { return x*x; });
    
    auto str=v3.toString();
    
    
}
