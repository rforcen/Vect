//
//  DynVect.cpp
//  Vect
//
//  Created by asd on 16/03/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//

//#define NDEBUG

#include "DynVect.hpp"
#include "Timer.h"
#include <vector>
#include <math.h>
#include <chrono>



void testDynVect() {
    
    size_t n=1e6;
    double x=10.001, err=1e-3;
    auto v=DynVect<double>::rnd(n);
    char lin[100];
    
//    assert( v.mtsort() == v.sort() );
   
    sprintf(lin, "DynVect test for %ld items\n",n);
    Timer().chrono(lin, [&](){

        Timer().chrono("sort & selection...\n", [&]() { // sort
            DynVect<double>::LambdaBool fbool = [](double x)->bool{return x>0.5 && x<0.6;};
            auto vsel=v(fbool);
            assert(v.sort().isSorted());
            assert(vsel.sort() == vsel.mtsort()); // mtsort experimental!
            assert(vsel.sort() == v(fbool).mtsort()); // no necessarily equal as order can be different
        });
        
    
        Timer().chrono("sum...\n", [&]() {
            double sm, ss, sa;
            auto ts=Timer().chrono("stsum...\n", [&]() {
                sm=v.stsum();
            });
            auto tm=Timer().chrono("sum...\n", [&]() {
                ss=v.sum();
            });
            auto ta=Timer().chrono("accumulate...\n", [&]() {
                sa=std::accumulate(v.begin(), v.end(), 0.0);
            });
            
            printf("(%f - %f = %f) s/m ratio: %.1f  acc:(%f, %ld)\n", sm, ss, abs(sm-ss), 1.*ts/tm, sa, ta);
            assert(abs(sm-ss) < err);
        });
        
        Timer().chrono("shuffle...\n", [&]() {
            auto vorg=v;
            assert(v.shuffle()!=v); // non mutator func
            assert(v.shuffle()!=vorg);
            assert(vorg == vorg.reverse().reverse());
        });
        
        Timer().chrono("find...\n", [&]() {
            auto vv=DynVect<double>::rnd(4000), vvv=vv;
            assert(vv == vvv.reverse().reverse());
            for (auto d:vv)
                assert(vvv.find(d) != -1);
            
            DynVect<double>::LambdaBool fbool = [](double x)->bool{return x>0.5 && x<0.5001;};
            auto vsel=v(fbool);
            for (auto i=0; i<std::min<size_t>(10, vsel.length()); i++)
                assert( v.find(vsel[i]) != -1);
        });
        
        Timer().chrono("iterator...\n", [&]() {
            double s=0., ss=0., vs=v.sum();
            for (auto d:v) s+=d;
            for (auto i=0; i<v.length(); i++) ss+=v[i];
            assert( abs(s-vs) < err && abs(ss-vs) < err && s==ss);
        });
        
      
        
        Timer().chrono("(func)...\n", [&](){
            DynVect<double>::Lambda _sin=[](double x)->double{ return sin(x); };
            
            auto vsel=v((DynVect<double>::LambdaBool)[](double x)->bool{ return x>=0 && x<=1;} );
            assert(vsel.count()>0);
            
            auto vf=v(_sin);
            assert(vf==v(_sin));
            
            DynVect<double>::LambdaBool fbool = [](double x)->bool{return x>0.5 && x<0.6;};
            vsel=v(fbool);
            printf("selected items: %ld\n",vsel.count());
            for (auto d:vsel) assert(fbool(d));
            
            assert(vsel.sort() == v(fbool).sort()); // no necessarily equal as order can be different
        });
        
        Timer().chrono("rnd performance in st/mt...\n", [=](){
            auto ts=Timer().chrono("single thread...",  [n](){ auto v=DynVect<double>::strnd(n); });
            auto tm=Timer().chrono("multi thread...",   [n](){ auto v=DynVect<double>::rnd(n);   });
            printf("s/m ratio: %.1f:\n", 1.*ts/tm);
        });
        
        Timer().chrono("fill...\n", [=](){
            DynVect<double>v1(n);
            auto ts=Timer().chrono("single thread...",  [&v1,n,x](){ v1.stfill(x); });
            auto tm=Timer().chrono("multiple thread...",  [&v1,n,x](){ v1=v1.fill(x+2); });
            
            printf("s/m ratio: %.1f:\n", 1.*ts/tm);
            assert(v1==x+2);
        });
        
        
        Timer().chrono("sum performance in st/mt...\n", [n,x](){
            auto v=DynVect<double>::rnd(n);
            double ss=0, sm=0;
            
            auto ts=Timer().chrono("single thread...", [&v, &ss](){ ss = v.stsum();  });
            auto tm=Timer().chrono("multi thread...",  [&v, &sm](){ sm = v.sum();    });
            
            printf("checking...");
            
            assert(abs(ss-sm)<1e-5);
            assert(v==v);
            assert(v!=v+1);
            
            printf("s/m ratio: %.1f:\n", 1.*ts/tm);
        });
        
        Timer().chrono("apply func...\n", [=](){
            DynVect<double>v1(n);
            auto ts=Timer().chrono("single thread...",  [&](){ v1.stapply([n, x]() -> double { return x; });});
            auto tm=Timer().chrono("multiple thread...", [&](){ v1=v1.apply([n, x]() -> double { return x*x; });});
            
            printf("s/m ratio: %.1f:\n", 1.*ts/tm);
            assert(v1==x*x);
        });
        
        Timer().chrono("iterator vs. index\n", [=](){
            
            auto ti=Timer().chrono("iterator...", [=](){
                DynVect<double>v1(n);
                v1.set(x);
                assert(v1==x && abs(v1.sum()-n*x) < err);
            });
            
            Timer().chrono("index...", [=](){
                DynVect<double>v1(n);
                v1=v1.fill(x);
                assert(v1==x && abs(v1.sum()-n*x)<err);
            });
            
            auto tf=Timer().chrono("fast index...", [=](){
                DynVect<double>v1(n);
                v1.setFast(x);
                assert(v1==x && abs(v1.sum()-n*x)<err);
            });
            
            
            printf("ratio iterator/fast index: %.2f\n",1. *ti/tf);
        });
        
        Timer().chrono("testing operators...", [](){
            size_t n=1e6;
            DynVect<double>dv, dv1(n), dv2(n);
            
            dv = dv1+1;
            assert(dv.sum() == n);
            
            dv2 = dv+1;
            
            assert(dv+1==dv2);
            assert(dv.sum() == dv2.sum()-n);
            
            assert(dv1+1 == dv1+1);
            dv1=dv1-1;
            
            dv.zero(); dv1.zero();
            dv=3; dv1=2; dv1.set(2);
            assert(dv*dv1 == dv*2);
            
            auto v=DynVect<double>::rnd(n), v1=DynVect<double>::rnd(n); // randoms
            assert(v+v != v1+v1);
            assert(v+v == v*2);
            
        });
        
        
        
        Timer().chrono("testing append, assign, ==, !=...", [](){
            DynVect<double>dv;
            
            size_t n=1e6;
            for (auto i=0; i<n; i++) dv << i << i+1; // items
            
            assert(dv.length() == n*2);
            for (auto i=0; i<n; i++) {
                assert(dv[i*2+0]==i);
                assert(dv[i*2+1]==i+1);
            }
            
            DynVect<double>dv1; // dynvect's
            dv1 << dv;
            
            assert(dv1==dv); // ==, !=
            assert(!(dv1!=dv));
            
            dv1 = dv; //  test assign
            auto dv2=dv;
            
            assert(dv1==dv);
            assert(dv2==dv);
            
            dv1.clear();
            dv2.clear();
            for (int i=0; i<100000; i++) dv1 << i; // this with more blocks than other
            for (int i=0; i<100; i++) dv2 << i;
            
            dv1=dv2; // dv1 has more blocks than dv2
            assert(dv1==dv2);
            
            dv=dv;
            assert(dv==dv);
            
            dv2 = dv;
            assert(dv2==dv);
            
        });
        
        Timer().chrono("testing constructor & index access...", [](){
            DynVect<double>dv, dv1(10000), dv2(10);
            
            double s=0, sv=0, sv1=0;
            for (auto i=0; i<dv1.length(); i++) //  test index mutable & accessor
            { dv1[i]=i; s+=i; sv+=dv1[i]; }
            
            for (auto i=0; i<dv1.length(); i++) sv1+=dv1[i];
            assert(s==sv && sv==sv1);
        });
        
    });
}
