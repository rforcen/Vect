//
//  Vect.cpp
//
//  Created by asd on 03/03/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//

#include "Vect.hpp"
#include <math.h>
#include <iostream>

#include <ctime>
#include <chrono>
#include <ratio>

using std::cout;

class Timer {
    std::chrono::high_resolution_clock::time_point begin, end;
public:
    void start() {
        begin = std::chrono::high_resolution_clock::now();
    }
    long lap() {
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> time_span = end - begin;
        return (long)time_span.count();
    }
};

void testVect() {
    puts("test started\n");
    
   
    
    printf("testing creators...");
    typedef double real;
    typedef Vect<real> VectReal;
    real sin(real x), cos(real), tan(real); // required to differenciate float/real versions
    Timer timer;
    
    VectReal v0, v1(1000), v2(1000), v5(v1), vxx(0, M_PI, 1e-4), vsin(0, M_PI, 1e-4, sin);
    auto v4=v1;
    puts("ok");
    
    
    printf("shuffle test...");
    {
        long ts, tm;
        real inc=1e-7;
        
        auto vorg=VectReal(0, 1, inc), vsh=vorg;
        
        timer.start(); vsh.stshuffle(); ts=timer.lap();
        auto vsht=vsh;
        
        assert(vsh.sum() == vorg.sum());
        
        vorg=VectReal(0, 1, inc); vsh=vorg;
        timer.start(); vsh.shuffle(); tm=timer.lap();
        
        printf("done, st=%ld, mt=%ld, ratio %.1f checking...", ts, tm, 1.*ts/tm);
        
        assert(vsh.sum() == vorg.sum());
        
        assert( vsh != vorg );
        for (auto d:vsh)    assert(vorg.bsearch(d)!=-1);
        
        vsh.sort();
        for (auto d:vorg)   assert(vsh.bsearch(d)!=-1);
        assert(vsh == vorg);
        
        vsht.sort();
        assert(vsht == vorg);
        
    }
    puts("ok");
    
    printf("single vs. multithead...");
    {
        long ts, tm;
        size_t n=1e7;
        
        // low mt gain when using simple funcs as rand, == , sum
        printf("for %ld values:\n", n);
        timer.start();        auto vs=VectReal::strnd(n);      ts=timer.lap();
        timer.start();        auto vm=VectReal::rnd(n);        tm=timer.lap();
        printf("lap for rnd init st:%ld ms, mt:%ld ms, ratio=%.1f\n", ts, tm, 1. * ts/tm);
        
        vs=vs;
        timer.start();        vs.stEQ(vs);        ts=timer.lap();
        timer.start();        void(vm==vm);       tm=timer.lap();
        printf("lap for == worst case init st:%ld ms, mt:%ld ms, ratio %.1f\n", ts, tm, 1.*ts/tm);
        
        timer.start();        auto vss=vs.stadd(vs);     ts=timer.lap();
        timer.start();        vss=vs+vm;                 tm=timer.lap();
        printf("lap for v+v worst case init st:%ld ms, mt:%ld ms, ratio: %.1f\n", ts, tm, 1.*ts/tm);
        
        // this starts better mt performnace
        timer.start();        auto vs1=vs.stfunc(sin);     ts=timer.lap();
        timer.start();        vs1=vs.apply(sin);           tm=timer.lap();
        printf("lap for sin st:%ld ms, mt:%ld ms, ratio %.1f\n", ts, tm, 1.*ts/tm);
        
        // the more complex func to apply the bugger st/mt ratio
        auto cmplxFunc=[](double x) -> double { return sin(x*x)+cos(x*2)+tan(pow(x,3)+log(1+abs(x))+exp(x/3.9)); };
        timer.start();        vs1=vs.stfunc(cmplxFunc);     ts=timer.lap();
        timer.start();        vs1=vs.apply(cmplxFunc);           tm=timer.lap();
        printf("lap for complex trigs. st:%ld ms, mt:%ld ms, ratio %.1f\n", ts, tm, 1.*ts/tm);
        
    }
    puts("ok");
    
    
    printf("testing operators...");
    timer.start();
    { // operators
        double n=1e7;
        VectReal v(n), v1=v+3.0;
        
        assert(v1.sum() == n*3);
        assert((v1-3).sum() == 0);
        assert(v1-3 == v);
        assert(v1*5 == v1/(1./5.));
        assert((v1^2) == v1*v1);
        assert((v1^3) == v1*v1*v1);
        
        v=0; // set all values to zero
        assert(v==0);
        
        v=2; // set all values to 2
        assert(v==2);
        
        VectReal vs(n), vb(n);
        vs=1; vb=2;
        assert(vs<vb);
        assert(vs <= vb);
        assert(vb>vs);
        assert(vb >= vs);
        assert(vb!=vs);
        assert(!(vb==vs));
        
        { // small sized vector aritmethics
            size_t sn=100;
            VectReal vs1(sn), vs2(sn);
            auto vsssum= vs1+1, _vs=vs1+vs2;
            assert(vsssum.sum() == vsssum.length());
            assert(_vs.sum() == 0);
            assert(vsssum == vs1+1);
            assert(_vs == vs1+vs2);
        }
        
        v=0;
        assert(v!=1);
        
        v+=3;
        assert(v == v1 && v == 3);
        v-=3;
        assert(v.sum()==0);
        v=v1;  v1*=2;
        assert(v1 == v*2);
        v=v1;  v1/=2;
        assert(v1 == v/2);
        
        v=v1; v+=v1;
        assert(v == v1*2);
        
        v=v1; v-=v1;
        assert(v.sum() == 0);
        
        v=v1; v*=v1;
        assert(v == v1*v1);
        
        v=v1; v/=v1;
        assert(v == v1/v1);
        
        auto vunion = v1 | v;
        assert((vunion == (v1 << v)) && (vunion.length() == v1.length()));
        auto vinter = v & v1;
        for (auto d:vinter)
            assert(v.bsearch(d) != -1 && v1.bsearch(d) != -1);
        
        printf("done for size %.0f in %ld ms\n", n, timer.lap());
    }
    
    printf("testing vector append / sub / join / cmp / rnd...");
    {
        size_t n=2000, n1=10000;
        auto v=VectReal::rnd(n1), v1=VectReal::rnd(n), v1tmp=v1, v2=(v1<<v);
        
        assert(v2(n, n+n1) == v); // sub vector test
        assert(v2==v1);
        assert(!(v2!=v1));
        assert((v1tmp|v) == v2); // join test |
    }
    puts("ok");
    
    printf("testing sequence / lambda constructor...");
    timer.start();
    for (int ic=0; ic<100; ic++)    {
        VectReal    vs(0, M_PI, 0.001),
        lsin(0.0, M_PI, 0.01, sin),
        lcos(0, M_PI_2, 0.001, cos); // from,to,inc, func
        
        v0 = VectReal(0, 1, 0.0001);
        
        // simulate graph data x-axis, y-lambda
        VectReal vx(0, M_PI, 1e-4), vy(vx, sin), vytest=vx.func(sin);
        
        assert(vy.max() == vytest.stmax());
        assert(vy.min() == vytest.stmin());
        assert(vy == vytest);
        assert(vy.norm() == vytest.norm());
        assert(vx.func(tan) == vx.func(sin) / vx.func(cos)); // tan = sin/cos
        
        if (ic==0) cout << vs.sum() << ", " << lsin.sum() << ", " << lcos.sum() << " time:" << timer.lap() <<  " ms ok\n";
    }
    
    puts("locate, seq generation st/mt"); // sequential generation st/mt, ratio=2.3
    {
        long ls, lm;
        real inc=1e-8;
        
        // seq st/mt ratio = 4 for 1e8 values
        printf("generating seq...");
        timer.start();   VectReal v=VectReal::stSeq(0, 1, inc);
        printf("st in %ld ms, ", ls=timer.lap());
        
        timer.start();  VectReal vm(0, 1, inc);  lm=timer.lap();
        printf("mt in %ld ms, ratio %.1f, ", lm, (double)ls/lm );
        auto diff = abs((v - vm).sum());
        printf("error=%e are eq:%d\n", diff, v==vm);
        
        
        printf("locate st...");
        timer.start();
        for (size_t i=v.count()-5; i<v.count(); i++) assert(v.stlocate(v[i])!=-1);
        printf("done in %ld ms, now in mt...", ls=timer.lap());
        timer.start();
        for (size_t i=v.count()-5; i<v.count(); i++) assert(v.locate(v[i])!=-1);
        lm=timer.lap();
        printf("done in %ld ms, rate st/mt = %.1f ", lm, (double)ls/lm );
        assert(v==vm);
    }
    puts("ok");
    
    puts("test mt vector min/max."); // ratio = 3.1
    {
        int n=1e8;
        Vect<double> v1=VectReal::rnd(n);
        real maxst, minst, maxmt, minmt;
        long lmmmt,lmt, lst;
        
        printf("mt min/max...");
        
        timer.start(); // mt
        auto minmaxmt = v1.minmax();
        printf("done minmax in %ld ms, now min,max mt...", lmmmt=timer.lap());
        
        timer.start(); // mt
        maxmt = v1.max();
        minmt = v1.min();
        printf("done in %ld ms, now st...", lmt=timer.lap());
        
        timer.start(); //  st
        maxst = v1.max(); minst = v1.min();
        lst=timer.lap();
        
        printf("done in %ld ms, ratio (st/mt)= %.1f ", lst, (double)lst/lmt);
        
        assert(maxmt == maxst && minmt == minst && minmaxmt.first == minst && minmaxmt.second == maxst);
    }
    puts("ok");
    
    puts("test mt vector mult."); // ratio = 1 not worth it
    {
        int n=1e8;
        Vect<double> v1=VectReal::rnd(n), v2(v1);
        long lmt, lst;
        
        printf("mt sum v1*v2...");
        
        timer.start(); // mt
        auto vmt=v1 / v2; // v1.MToper(v2, VectReal::opMUL);
        printf("done in %ld ms, now st...", lmt=timer.lap());
        
        timer.start(); //  st
        auto vst=v1 / v2;
        lst=timer.lap();
        
        printf("done in %ld ms, ratio (st/mt)= %.1f ", lst, (double)lst/lmt);
        
        assert(vmt == vst);
    }
    puts("ok");
    
    
    
    puts("testing init vector, iterator...");
    {
        VectReal vinit={0,1,2,3,4,3,2,1,0};
        real f;
        assert(vinit.length() == 9);
        for (auto d:vinit) { // must use iterator to make _inidex=0;
            (void)d; // not used
            vinit >> f;
            cout << f << ",";
        }
    }
    puts("ok");
    
    printf("testing bsearch...");
    {
        double inc=1e-6;
        VectReal v(0,1, inc);
        v.sort();
        
        timer.start();
        int iters=0;
        for (double d=0.2; d<0.4; d+=inc) {
            assert(v.bsearch(d)!=-1);
            iters++;
        }
        printf("%d iterations in %ld ms, ",iters, timer.lap());
    }
    puts("ok");
    
    
    // shuffle
    printf("shuffle test...");
   
    timer.start();   auto vorg=VectReal(0, 1, 1e-4), vsh=vorg;    vsh.shuffle();
    assert( vsh != vorg );
    for (auto d:vsh)
        assert(vorg.locate(d)!=-1);
    
    cout << "ok, shuffle sum error=" << abs(vorg.sum() - vsh.sum()) << " due to kahan algo.\n";
    cout << " sums: shuffle sum=" << vsh.sum() << " , original sum=" << vorg.sum() << ", kahan sum=" << vsh.sum() << "\n";
    assert( abs(vsh.sum() - vorg.sum()) == 0); //  this is only possible using kahan algo. in other case diff != 0
    cout << "duration: " << timer.lap() << "ms\n";
    
    // multithreading
    printf("testing multitheading func...");
    {
        long t1,t2;
        
        VectReal vmt(0, 1, 1e-7);
        auto foo = [](real x){ return sin(x);};
        
        timer.start();
        auto vmtf = vmt.func(foo); // the mt version
        cout << "mt done in " << timer.lap() << "ms, now st...";
        
        timer.start();
        auto vff=vmt.stfunc(foo); // st version... see the timing difference?
        assert(vff==vmtf);
        
        cout << "done in " << timer.lap() << "ms\nnow mtSum...";
        
        timer.start(); // mtSum vs. sum
        auto smt=vff.sum();
        t1=timer.lap();
        printf("done in %ld ms, now st...", t1);
        timer.start();
        vff.stsum();
        t2=timer.lap();
        printf("done in %ld ms, ratio =%.1f",t2, 1.*t2/t1);
        assert(smt == vff.sum()); //  ok due to kahan sum algo.
        puts("ok");
    }
    
    {
        printf("testing append...");
        for (int ic=0; ic<200; ic++) {
            v0.clear();
            int n=100000; // append test
            real rs=0;
            for (int i=0; i<n; i++) { v0 << i; rs+=i; }
            v0.fit(); // set memory to size
            assert(v0.sum()==rs && "append failure");
        }
        puts("ok");
    }
    
    printf("testing erase()...");
    for (int i=0; i<80; i++) {
        const int ns=2000, nr=900;
        VectReal v = VectReal::seq(ns), v0=VectReal::rnd(nr);
        
        for (int j=0; j<ns; j++) // locate is mt and will last longer than st version
            assert( v.stlocate( v.erase(rand() % v.count()) )==-1 );
        for (int j=0; j<nr; j++)
            v0.stlocate( v0.erase(rand() % v0.count()) ); // possible dups (random)
        
        assert(v.count()==0 && v0.count()==0);
    }
    puts("ok");
    
    // test prefix ++, --
    cout << "testing ++, -- ...";
    ++v0;
    --v0;
    cout << "ok\n";
    
    // seq & Vect[Vect++] vect indexing post increment
    cout << "testing seq, iterator, sum, count...";
    auto vseq=VectReal::stSeq(0, 100, 0.1), vs0=vseq;
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
    auto vsub = v0(10,20);
    assert(vsub.length()==10 && "sub failure");
    for (auto d: vsub)
        assert(v0.locate(d)>=10);
    
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
    v6=VectReal::seq(5,1);
    
    // lambda apply func / sort
    v2.apply(sin).sort().apply([](real x) -> real { return x*x; });
    
    // mt vs. st apply -> ratio is 2.5 approx.
    {
        printf("testing apply in mt vs. st modes, now in MT...");
        
        VectReal v=VectReal(0,1,1e-8), v1=v;
        timer.start(); // mt
        v.apply(sin);
        auto tmt=timer.lap();
        printf("done in %ld ms, now in ST...", tmt);
        
        timer.start(); // st
        v1.stapply(sin);
        auto tst=timer.lap();
        printf("done in %ld ms, ratio %.1f\n", tst, 1.0 * tst/tmt);
        assert(v1==v);
    }
    auto str=v3.toString();
    
    puts("\ntest completed OK\n");
    fflush(stdout);
}
