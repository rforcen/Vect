//
//  Vect.cpp
//
//  Created by asd on 03/03/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//

#include "Vect.hpp"
#include <math.h>
#include <iostream>

using std::cout;


void testVect() {
    typedef double real;
    typedef Vect<real> VectReal;
    real sin(real x), cos(real), tan(real); // required to differenciate float/real versions
    
    VectReal vinit(5, new real[5]{0,1,2,3,4});
    VectReal v0, v1(100), v2(100);
    auto v4=v1;
    
    cout << "testing >> ...";
    real f;
    for (auto d:vinit) { // must use iterator to make _inidex=0;
        (void)d; // not used
        vinit >> f;
        cout << f << ",";
    }
    cout << "ok\n";
    
    cout << "testing sequence / lambda constructor...";
    for (int ic=0; ic<100; ic++)    {
        VectReal vs(0, M_PI, 0.001), lsin(0.0, M_PI, 0.01, sin), lcos(0, M_PI_2, 0.001, cos); // from,to,inc, func
        v0 = VectReal(0, 1, 0.0001);
        
        // simulate graph data x-axis, y-lambda
        VectReal vx(0, M_PI, 1e-4), vy(vx, sin), vytest=vx.func(sin);
        assert(vy.norm() == vytest.norm());
        assert(vx.func(tan) == vx.func(sin) / vx.func(cos)); // tan = sin/cos
        
        if (ic==0) cout << vs.sum() << ", " << lsin.sum() << ", " << lcos.sum() << " ok\n";
    }
    // shuffle
    cout << "shuffle test...";
    auto vorg=VectReal(0, 1, 1e-4), vsh=vorg;
    vsh.shuffle();
    assert( vsh != vorg );
    for (auto d:vsh)
        assert(vorg.locate(d)!=-1);
    
    cout << "ok, shuffle sum error=" << abs(vorg.sum() - vsh.sum()) << " due to kahan algo.\n";
    cout << " sums: shuffle sum=" << vsh.sum() << " , original sum=" << vorg.mtSum() << ", kahan sum=" << vsh.sum() << "\n";
    assert( abs(vsh.sum() - vorg.sum()) == 0); //  this is only possible using kahan algo. in other case diff != 0
    
    // multithreading
     cout << "testing multitheading func...";
    for (int ic=0; ic<3; ic++) {
        VectReal vmt(0, 1, 1e-4);
        auto foo = [](real x){ return sin(x);};
        auto vmtf = vmt.mtFunc(foo); // the mt version
        cout << "mt done, now st...";
        auto vff=vmt.func(foo); // st version... see the timing difference?
        assert(vff==vmtf);
        cout << "done\nnow mtSum...";
        auto smt=vff.mtSum();
        cout << "done, now st..."<< abs(vff.sum() - smt);
//        assert(smt == vff.sum()); //  can't assurance due to double sum precission error
        cout << " ok\n";
    }
    
    
    cout << "testing append...";
    for (int ic=0; ic<200; ic++) {
        v0.clear();
        int n=100000; // append test
        real rs=0;
        for (int i=0; i<n; i++) { v0 << i; rs+=i; }
        v0.fit(); // set memory to size
        assert(v0.sum()==rs && "append failure");
    }
    cout << "ok\n";
    
    cout << "testing erase()...";
    for (int i=0; i<80; i++) {
        const int ns=2000, nr=900;
        VectReal v = VectReal::seq(ns), v0=VectReal::rnd(nr);
        
        for (int j=0; j<ns; j++)
            assert( v.locate( v.erase(rand() % v.count()) )==-1 );
        for (int j=0; j<nr; j++)
            v0.locate( v0.erase(rand() % v0.count()) ); // possible dups (random)
        
        assert(v.count()==0 && v0.count()==0);
    }
    cout << "ok\n";
    
    // test prefix ++, --
    cout << "testing ++, -- ...";
    ++v0;
    --v0;
    cout << "ok\n";
    
    // seq & Vect[Vect++] vect indexing post increment
    cout << "testing seq, iterator, sum, count...";
    auto vseq=VectReal::seq(0, 100, 0.1), vs0=vseq;
    for (auto d:vseq) vs0[vseq++]=d+1; // can access & increment iterator index
    assert(vseq.sum() - ( vs0.sum()-vs0.count() ) == 0); //  again kahan algo makes the magic preservin precission
    cout << "ok, error=" << vseq.sum() - ( vs0.sum()-vs0.count() )  << "\n";
    
    // filterIndex, v[lambda]
    cout << "testing filterIndex, [by VectIndex]...";
    for (int i=0; i<100; i++)  {
        int n=10000;
        auto fsel = [](real x) -> bool { return x>0.25 && x<0.85; };
        VectReal v=VectReal::rnd(n);
        
        auto ixs = v.filterIndex(fsel); // v[fsel]
        
        assert(ixs.count() == v.filter(fsel).count());
        assert(v[ixs]==v.filter(fsel));
        
        for (auto ix:ixs)
            assert(fsel(v[ix]));
        auto vsel=v[ ixs ]; // index by VectIndex
        
        assert(vsel.filter(fsel)==vsel);
        assert( v[fsel] == v[v.filterIndex(fsel)] ); // index by lambda
        
    }
    cout << "ok\n";
    
    // aritmethic
    cout << "testing random, algebra...";
    v1.random(); v2.random(); v4.random();
    v0 = (v1+v2/v4)*3;
    cout << "ok\n";
    
    cout << "testing ::rnd, <<, sum, logical ops., iterators, filter, reverse, sub, genWave, apply...";
    v0.clear(); // vect append
    
    auto v01=VectReal::rnd(40000), v02=VectReal::rnd(5000);
    v0 << v01;
    v0 << v02;
    
    // for bigger v01, v02 sizes precission real may just not enough
    cout << "error =" << fabs(v0.sum() - (v01.sum()+v02.sum())) << "\n";
    assert(fabs(v0.sum() - (v01.sum()+v02.sum())) < 1e-8);
    
    
    for (int i=0; i<v1.length(); i++) // test index mutator & accesor. (v2=v1)
        v2[i] = v1[i];
    // logical ops
    assert(v1==v1     && "== failure");
    assert(v1==v2     && "!= failure");
    assert(!(v1>v2)   && ">  failure");
    
    real sum=0;
    for (auto vi:v1) sum+=vi; // iterator
    assert(sum==v1.sum() && "iterator failure");
    
    // reverse
    auto vrev = VectReal::seq(10), vseq1=VectReal::seq(10, 0.02), vw=vrev;
    assert(vrev == vw.reverse().reverse() && "reverse failure");
    vseq = vseq1 * 0.01;
    
    // subvector (from, to)
    auto vsub = v0.sub(10,20);
    assert(vsub.length()==10 && "sub failure");
    
    // filter
    auto fsel = [](real x) -> bool { return x>0.5 && x<0.7; };
    auto vf = v0.filter(fsel).sort();
    for (auto d:vf) assert(fsel(d) && "filter failure");
    
    vf=0; // set all values to zero & check scalar
    assert(vf==0    && "== scalar failure");
    assert(vf!=10   && "== scalar failure");
    assert(vf>-10   && "> scalar failure");
    assert(vf<10    && "< scalar failure");
    
    // generate complex wave
    auto vwave = VectReal::genWave(1, 8000, 3,
                                   new real[3]{1,0.4,0.6},
                                   new real[3]{200, 600, 1800},
                                   new real[3]{0,0,0});
    
    VectReal v3(v0); // misc
    auto v6=v3.norm();
    assert( v6.filter([](real x) -> bool { return x>=1.0 && x<0.0; }).count() == 0 );
    v6.sequence(0,5);
    
    // lambda apply func / sort
    v2.apply(sin).sort().apply([](real x) -> real { return x*x; });
    
    auto str=v3.toString();
    
    cout << "test completed OK\n";
    fflush(stdout);
    
}
