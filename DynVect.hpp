//
//  DynVect.hpp
//  Vect
//
//  Created by asd on 16/03/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//
// solves the realoc issue in linear vector using blocks of 2^nBits size
// _foo -> mutable

#ifndef DynVect_hpp
#define DynVect_hpp

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <functional>
#include <string>
#include <thread>
#include <vector>
#include <future>
#include <random>

#include "cpuRandom.hpp"

// #define NDEBUG // optimize performance

template <class T> class DynVect {
public:
    typedef size_t Index;
    typedef std::function<T(T)>const& Lambda;
    typedef std::function<T()>const& LambdaNoArg;
    typedef std::function<void()>const& LambdaVoid;
    typedef std::function<bool(T)>const& LambdaBool;
    typedef std::function<void(size_t, size_t)>const& LambdaRange;
private:
    const int nBits=14; // define 2^nBits block size
    const size_t maskIndex=(__SIZE_MAX__<<nBits) ^ __SIZE_MAX__;
    const size_t blockSize=1<<nBits;
    const size_t blockSizeBytes = blockSize * sizeof(T);
    
    static const int minsizeMT=100000; // min vector size for multitheading
    
    size_t size=0;
    size_t usedBlocks=0;
    
    T**blocks=nullptr; // pointer to blocks each of 2^nBits size
    
    static T cpuRand() {
        static bool hasRandomCPU=CPUrandom::hasSupport();
        return hasRandomCPU ? CPUrandom::randd() : ::rand();
    }
    
    // memory mgr
    T*newBlock() {
        auto p = (T*)calloc(blockSize, sizeof(T));
        assert(p!=nullptr);
        return p;
    }
    void resizeBlocks(size_t n) {
        blocks=(T**)realloc(blocks, n * sizeof(T*));
        assert(blocks!=nullptr);
    }
    
    void copyBlock(T*bdst, T*borg) { memcpy(bdst, borg, blockSizeBytes); }
    
    inline size_t i2b(size_t index) const { return index >> nBits; } // block # of index
    inline size_t i2i(size_t index) const { return index & maskIndex; } // index in block of global index
    
public:
    class Kahan {
        
        T s=0, c=0;
    public:
        void add(T v) {
            auto y = v - c, t = s + y;
            c = (t - s) - y;
            s = t;
        }
        T sum() { return s; }
    };
    
    T kahanSum() const { // sum vector range usign kahan algorithm
        Kahan kh;
        for (auto b=0; b<usedBlocks; b++) // the fast way
            for (auto p=blocks[b]; p<blocks[b] + ((b!=usedBlocks-1) ? blockSize : (size % blockSize)) ; p++)
                kh.add(*p);
        return kh.sum();
    }
    inline T kahanSum(size_t from, size_t to) const { //
        Kahan kh;
        
        // first block
        for (auto pb=blocks[i2b(from)]+i2i(from), p=pb; p<pb+(blockSize - i2i(from)); p++) kh.add(*p);
        for (auto b=i2b(from)+1; b<i2b(to); b++) // rest of complete blocks
            for (auto pb=blocks[b], p=pb; p<pb + blockSize ; p++) kh.add(*p);
        for (auto pb=blocks[i2b(to)], p=pb; p<pb+i2i(to); p++) kh.add(*p); // last block
        return kh.sum();
    }
    
    inline void applyFunc(size_t from, size_t to, Lambda lambda) const { // in 3 stages: 1st block, middle blocks, last block.
        // first block
        for (auto pb=blocks[i2b(from)]+i2i(from), p=pb; p<pb+(blockSize - i2i(from)); p++) *p=lambda(*p);
        for (auto b=i2b(from)+1; b<i2b(to); b++) // middle, rest of complete blocks
            for (auto pb=blocks[b], p=pb; p<pb + blockSize ; p++) *p=lambda(*p);
        for (auto pb=blocks[i2b(to)], p=pb; p<pb+i2i(to); p++) *p=lambda(*p); // last block
    }
    
public:
    
    // constructors
    DynVect() {}
    DynVect(size_t size) : size(size) {
        usedBlocks = size/blockSize + 1;
        // create 'usedBlocks' blocks
        resizeBlocks(usedBlocks);
        for (auto i=0; i<usedBlocks; i++) blocks[i]=newBlock();
    }
    DynVect(const DynVect&v) : size(v.size) { // required for: auto v=v1, v2(v);
        usedBlocks = size/blockSize + 1;
        resizeBlocks(usedBlocks); // create 'usedBlocks' blocks
        for (auto i=0; i<usedBlocks; i++) { // copy 'v' content
            blocks[i]=newBlock();
            copyBlock(blocks[i], v.blocks[i]);
        }
    }
    
    ~DynVect(){
        clear();
    }
    
    static DynVect strnd(size_t size) { // nomralized random
        DynVect v(size);
        for (auto b=0; b<v.usedBlocks; b++)
            for (auto p=v.blocks[b]; p<v.blocks[b] + ((b!=v.usedBlocks-1) ? v.blockSize : (v.size % v.blockSize)) ; p++) *p=(T)cpuRand(); ///RAND_MAX;
        return v;
    }
    static DynVect rnd(size_t size) { // nomralized random
        return (size < minsizeMT) ?
        strnd(size) :
        *Thread(new DynVect(size), (Lambda)([](T x)->T{ return (T)cpuRand(); })).vect();
    }
    
    inline size_t length() const { return size; }
    inline size_t count() const { return size; }
    
    void clear() {
        if(blocks!=nullptr) {
            for (auto i=0; i<usedBlocks; i++) free(blocks[i]);
            free(blocks);
        }
        blocks=nullptr;
        usedBlocks = size = 0;
    }
    void zero() {  fill(0);  }
    void one() { fill(1); }
    void stfill(T c) {
        for (auto b=0; b<usedBlocks; b++)
            for (auto p=blocks[b]; p<blocks[b] + ((b!=usedBlocks-1) ? blockSize : (size % blockSize)) ; p++) *p=c;
    }
    DynVect&_fill(T c) {
        if (size<minsizeMT) stfill(c);
        else                Thread(this, (Lambda)([c](T) -> T { return c; }));
        return *this;
    }
    DynVect fill(T c) {  return DynVect(*this)._fill(c);    }
    
    T stsum() const {
        return kahanSum();
    }
  
    void debug(size_t n=20) {
        for (auto i=0; i<n; i++)
            printf("%.2f, ", (*this)[i]);
    }
    
    class Thread {
    public:
        const int nthreads = std::thread::hardware_concurrency()-1;
        const DynVect*v;
        std::mutex mutex;
        
        Thread(const DynVect*v) : v(v) {}
        Thread(const DynVect*v, LambdaRange lambda) : v(v) {
            auto segsize=v->size / nthreads;
            std::thread ths[nthreads];
            
            for (auto t=0; t<nthreads; t++)
                ths[t] = std::thread ( [this, t, segsize, lambda]()
                                      { lambda(t*segsize, ((t==nthreads-1) ? this->v->size : ((t+1)*segsize))); });
            for (auto t=0; t<nthreads; t++)
                ths[t].join();
        }
        Thread(const DynVect*v, Lambda lambda) : v(v) {
            auto segsize=v->size / nthreads;
            std::thread ths[nthreads];
            
            for (auto t=0; t<nthreads; t++)
                ths[t] = std::thread ( [this, t, segsize, lambda]()
                                      {
                                          auto v=this->v; auto mask=v->maskIndex;
                                          for (auto i=t*segsize; i < ((t==nthreads-1) ? v->size : ((t+1)*segsize)); i++) {
                                              auto &d=v->blocks[i >> v->nBits][i & mask];
                                              d=lambda(d);
                                          }
                                      });
            for (auto t=0; t<nthreads; t++) ths[t].join();
        }
        Thread(const DynVect*v, DynVect&other, LambdaBool lambda) : v(v) {
            auto segsize=v->size / nthreads;
            std::thread ths[nthreads];
            
            for (auto t=0; t<nthreads; t++)
                ths[t] = std::thread ( [this, t, segsize, lambda, &other]()
                      {
                          auto v=this->v; auto mask=v->maskIndex;
                          auto vl=DynVect();
                          for (auto i=t*segsize; i < ((t==nthreads-1) ? v->size : ((t+1)*segsize)); i++) {
                              auto d=v->blocks[i >> v->nBits][i & mask];
                              if(lambda(d)) vl << d;
                          }
                          // lock append th partial result
                          mutex.lock();  other << vl;   mutex.unlock();
                      });
            for (auto t=0; t<nthreads; t++) ths[t].join();
        }
        T ret(T val) { return val; }
        const DynVect *vect() { return v; }
        
        bool eq(const DynVect&other) {
            volatile bool eq=true;
            auto segsize=v->size / nthreads;
            std::thread ths[nthreads];
            
            for (auto t=0; t<nthreads; t++)
                ths[t] = std::thread ( [=, &eq]()  {
                      auto v=this->v;
                      auto mask=v->maskIndex;
                      auto nBits=v->nBits;
                      auto size=v->size;
                      
                      for (auto i=t*segsize; i < ((t==nthreads-1) ? size : ((t+1)*segsize)) && eq; i++) {
                          auto d=v->blocks[i >> nBits] + (i & mask);
                          auto dother=other.blocks[i >> nBits] + (i & mask);
                          if (*d != *dother) { eq=false; break; }
                      }
                  });
            for (auto t=0; t<nthreads; t++) ths[t].join();
            return eq;
        }
        
        size_t find(const T c) {
            volatile auto found=false;
            volatile size_t res=-1;
            
            auto segsize=v->size / nthreads;
            std::thread ths[nthreads];
            
            for (auto t=0; t<nthreads; t++)
                ths[t] = std::thread([this, t, c, segsize, &found, &res]()   {
                    auto v=this->v;
                    auto mask=v->maskIndex;
                    auto nBits=v->nBits;
                    auto size=v->size;
                   
                    
                    for (auto i=t*segsize; i < ((t==nthreads-1) ? size : ((t+1)*segsize)) && !found; i++) {
                        auto d=v->blocks[i >> nBits] + (i & mask);
                        if (*d == c) {
                            res=i; found=true; break;
                        };
                    }
                });

            for (auto t=0; t<nthreads; t++) ths[t].join();
            return res;
        }
    };
    
    T sum() const {
        Kahan kh;
        
        return (size<minsizeMT) ? kahanSum() : Thread(this, [&](size_t from, size_t to) {
            kh.add( kahanSum(from, to) );
        }).ret(kh.sum());
    }
    
   
   
    class Iterator : public std::iterator<std::random_access_iterator_tag, T, ptrdiff_t, T*, T&> {
    public:
        
        DynVect<T>*v=nullptr;
        long index=0;
        
        Iterator() = default;
        Iterator(DynVect<T>*v) : v(v) {}

        Iterator begin() {
            Iterator it(v);
            it.index=0;
            return it;
        }
        Iterator end() {
            Iterator it(v);
            it.index=v->length();
            return it;
        }
        inline T& operator*() {     return (*v)[index];    }
        inline Iterator& operator++()    { ++index;    return *this; } // ++iter
        inline Iterator& operator++(int) { index++;    return *this; } // iter++
        inline Iterator& operator--()    { --index;    return *this; } // --iter
        inline Iterator& operator--(int) { index--;    return *this; } // iter--
        
        inline Iterator& operator+=(long d) { index+=d; return *this;}
        inline Iterator& operator-=(long d) { index-=d; return *this;}
        inline Iterator  operator+(long d) { Iterator it(*this); it.index+=d; return it; }
        inline Iterator  operator-(long d) { Iterator it(*this); it.index-=d; return it; }
        
        inline long operator-(const Iterator&other) { return index-other.index; }
        inline long operator+(const Iterator&other) { return index+other.index; }
        
        inline bool operator!=(const Iterator& other) const {   return other.index!=index;    }
        inline bool operator==(const Iterator& other) const {   return other.index==index;    }
        
        inline bool operator < (const Iterator&other) const {  return index<other.index;    }
        inline bool operator > (const Iterator&other) const {  return index>other.index;    }
        inline bool operator <= (const Iterator&other) const { return index<=other.index;   }
        inline bool operator >= (const Iterator&other) const { return index>=other.index;   }
    };

    Iterator begin()    { return Iterator(this).begin(); }
    Iterator end()      { return Iterator(this).end();   }
    
    T parSum(Iterator beg, Iterator end) { // experimental std::async
        auto len = end - beg;
        if (len < 1024*32)
            return (T) std::accumulate(beg, end, (T)0);
        
        Iterator mid = beg + len/2;
        auto handle = std::async(std::launch::async,[=]()->T{ return parSum(mid, end); });

        return parSum(beg, mid) + handle.get();
    }
    
    int multiThreadSort(Iterator begin, Iterator end) {
        const int szPage=4*1024;
        auto const sz = end - begin;
        if (sz <= 1) return 0;
        
        auto pivot = begin + sz/2;
        auto const pivot_v = *pivot;
        
        std::swap(*pivot, *(end - 1));
        auto p = std::partition(begin, end, [&](const T& a) { return a < pivot_v; } );
        std::swap(*p, *(end - 1));
        
        if (sz > szPage) {
            auto handle = std::async(std::launch::async, [&]() {
                return multiThreadSort(begin, p);
            });
            multiThreadSort(p + 1, end);
        } else {
            multiThreadSort(begin, p);
            multiThreadSort(p + 1, end);
        }
        return 0;
    }
    
    DynVect& _mtsort() { // mutable experimental use w/caution
        multiThreadSort( begin(), end() );
        return *this;
    }
    DynVect mtsort() { // non mutable experimental use w/caution
        return DynVect(*this)._mtsort();
    }
    
    DynVect&_sort() { // mutable _sort
        std::sort( begin(), end() );
        return *this;
    }
    DynVect sort() { // non mutable sort
        return DynVect(*this)._sort();
    }
    
    bool isSorted() {
        return std::is_sorted(begin(), end());
    }
    
    DynVect&_shuffle() {
        std::random_device rd;
        std::mt19937 g(rd());
        
        std::shuffle(begin(), end(), g);
        return *this;
    }
    DynVect shuffle() {
        return DynVect(*this)._shuffle();
    }
    DynVect&_reverse() {
        std::reverse(begin(), end());
        return *this;
    }
    DynVect reverse() {
        return DynVect(*this)._reverse();
    }
    size_t find(const T c) {
        if (size<minsizeMT) {
            auto ri = std::find( begin(), end(), c);
            return ri == end() ? -1 : end()-ri;
        } else return Thread(this).find(c);
    }
    
    DynVect&set(T c){ // iterator usage
        for (auto &d:*this) d=c;
        return *this;
    }
    
    void setFast(T c) { // fastest iteration possible
        for (auto b=0; b<usedBlocks; b++)
            for (auto p=blocks[b]; p<blocks[b] + ((b!=usedBlocks-1) ? blockSize : (size % blockSize)) ; p++) *p=c;
    }
    
    DynVect&stapply(Lambda foo) { // mutable versions
        for (auto &d:*this) d=foo(d);
        return *this;
    }
    DynVect&stapply(LambdaNoArg foo) {
        for (auto &d:*this) d=foo();
        return *this;
    }
    DynVect&_apply(Lambda foo) {
        if (size<minsizeMT) for (auto &d:*this) d=foo(d);
        else Thread(this, foo);
        return *this;
    }
    DynVect&_apply(LambdaNoArg foo) {
        for (auto &d:*this) d=foo();
        return *this;
    }
    // non mutable
    DynVect apply(Lambda foo) { return DynVect(*this)._apply(foo);  }
    DynVect apply(LambdaNoArg foo) { return DynVect(*this)._apply(foo);  }
    
    
    DynVect &operator=(const DynVect &other) { // asignment
        if (&other == this) return *this; // check for self-assignment to avoid delloc problems
        
        if (usedBlocks < other.usedBlocks) { // other > this -> resize
            resizeBlocks(other.usedBlocks);
            for (auto b=0; b<other.usedBlocks; b++) { // add other.size-size more blocks & copy
                if(b>=usedBlocks)  // additional space required?
                    blocks[b]=newBlock();
                copyBlock(blocks[b], other.blocks[b]);
            }
        } else { // this >= other ->  then only copy blocks content & release rest
            for (auto b=0; b<usedBlocks; b++)
                if(b<other.usedBlocks)  copyBlock(blocks[b], other.blocks[b]);
                else free(blocks[b]); // release non used blocks
        }
        size=other.size;
        usedBlocks=other.usedBlocks;
        
        return *this;
    }
    DynVect &operator=(const T c) { // asignment
        fill(c);
        return *this;
    }
    
    inline T& operator[](size_t i) { // index mutator
        assert(i<size && "Vect: index error");
        return blocks[i >> nBits][i & maskIndex];
    }
    inline const T& operator[](size_t i) const { // index accessor
        assert(i<size  && "Vect: index error");
        return blocks[i >> nBits][i & maskIndex];
    };
    
    DynVect operator()(Lambda lambda) {
        DynVect v(*this);
        if(size < minsizeMT) for (auto &d:v) d=lambda(d);
        else Thread(&v, lambda);
        return v;
    }
    DynVect operator()(LambdaNoArg lambda) {
        DynVect v(*this);
        if(size < minsizeMT) for (auto &d:v) d=lambda();
        else Thread(&v, lambda);
        return v;
    }
    DynVect operator()(LambdaBool lambda) {
        DynVect v;
        if(size < minsizeMT) { for (auto &d:v) if(lambda(d)) v<<d; }
        else Thread(this, v, lambda);
        return v;
    }
    
    DynVect& operator<<(T c) { //  append item
        if((size & maskIndex) == 0) { // fits? -> no, add 1 more block
            resizeBlocks(usedBlocks+1); // expand blocks
            blocks[usedBlocks]=newBlock(); // add 1 more
            usedBlocks++;
        }
        (*this)[size++]=c;
        return *this;
    }
    
    DynVect& operator<<(const DynVect& v) { //  append other DynVect
        for (Index i=0; i<v.length(); i++) (*this) << v[i];
        return *this;
    }
    
    // compare DynVect
    bool operator==(const DynVect&other) const { //  all eq
        bool eq=(size==other.size);
        if(eq) {
            if (size<minsizeMT) {
                for (Index i=0; i<size; i++)
                    if((*this)[i]!=other[i]) { eq=false; break; }
            } else
                eq=Thread(this).eq(other);
        }
        return eq;
    }
    bool operator!=(const DynVect&other) const { //  any neq
        bool neq=(size!=other.size);
        if(!neq)
            for (Index i=0; i<size; i++)
                if((*this)[i]!=other[i]) { neq=true; break; }
        return neq;
    }
    // compare 'c'
    bool operator==(const T&c) const { //  all eq
        bool eq=true;
        for (Index i=0; i<size; i++)
            if((*this)[i]!=c) { eq=false; break; }
        return eq;
    }
    bool operator!=(const T&c) const { //  all eq
        return !(*this==c);
    }
    
    // operator aritmethic c
    DynVect operator+(const T c) const { DynVect v(this->size);  for (Index i=0; i<size; i++) v[i] = (*this)[i]+c;    return v; }
    DynVect operator-(const T c) const { DynVect v(this->size);  for (Index i=0; i<size; i++) v[i] = (*this)[i]-c;    return v; }
    DynVect operator*(const T c) const { DynVect v(this->size);  for (Index i=0; i<size; i++) v[i] = (*this)[i]*c;    return v; }
    DynVect operator/(const T c) const { DynVect v(this->size);  for (Index i=0; i<size; i++) v[i] = (*this)[i]/c;    return v; }
    
    DynVect& operator+=(const T c) {  for (Index i=0; i<size; i++) (*this)[i]+=c;    return *this; }
    DynVect& operator-=(const T c) {  for (Index i=0; i<size; i++) (*this)[i]-=c;    return *this; }
    DynVect& operator*=(const T c) {  for (Index i=0; i<size; i++) (*this)[i]*=c;    return *this; }
    DynVect& operator/=(const T c) {  for (Index i=0; i<size; i++) (*this)[i]/=c;    return *this; }

    // operator aritmetic DynVect
    DynVect operator+(const DynVect& o) const {
        assert(size==o.size);
        DynVect v(this->size); for (Index i=0; i<size; i++) v[i] = (*this)[i]+o[i];    return v; }
    DynVect operator-(const DynVect& o) const {
        assert(size==o.size);
        DynVect v(this->size); for (Index i=0; i<size; i++) v[i] = (*this)[i]-o[i];    return v; }
    DynVect operator*(const DynVect& o) const {
        assert(size==o.size);
        DynVect v(this->size); for (Index i=0; i<size; i++) v[i] = (*this)[i]*o[i];    return v; }
    DynVect operator/(const DynVect& o) const {
        assert(size==o.size);
        DynVect v(this->size); for (Index i=0; i<size; i++) v[i] = (*this)[i]/o[i];    return v; }
    
};
#endif /* DynVect_hpp */
