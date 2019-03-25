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


using std::thread, std::vector, std::function, std::generate,
std::mutex, std::atomic, std::iterator, std::accumulate,
std::random_access_iterator_tag, std::launch, std::random_device,
std::mt19937, std::swap;

// #define NDEBUG // optimize performance

template <class T> class DynVect {
public:
    typedef size_t Index;
    typedef function<T(T)>const& Lambda;
    typedef function<T()>const& LambdaNoArg;
    typedef function<void()>const& LambdaVoid;
    typedef function<bool(T)>const& LambdaBool;
    typedef function<void(size_t, size_t)>const& LambdaRange;
    
    typedef enum { opADD, opSUB, opMUL, opDIV, opPOW } Operator;
    
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
        void reset() { s = c = 0; }
        
        inline T add(T v) {
            auto y = v - c, t = s + y;
            c = (t - s) - y;
            s = t;
            return s;
        }
        inline T total() { return s; }
    };
    
    
    
public:
    
    // constructors
    DynVect() {}
    DynVect(size_t size) : size(size) {
        if(size>0) {
            usedBlocks = size/blockSize + 1;
            // create 'usedBlocks' blocks
            resizeBlocks(usedBlocks);
            for (auto i=0; i<usedBlocks; i++) blocks[i]=newBlock();
        }
    }
    DynVect(const DynVect&v) : size(v.size) { // required for: auto v=v1, v2(v);
        if(size>0) {
            usedBlocks = size/blockSize + 1;
            resizeBlocks(usedBlocks); // create 'usedBlocks' blocks
            for (auto i=0; i<usedBlocks; i++) { // copy 'v' content
                blocks[i]=newBlock();
                copyBlock(blocks[i], v.blocks[i]);
            }
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
    
    T kahanSum() const { // sum vector range usign kahan algorithm
        Kahan kh;
        for (auto b=0; b<usedBlocks; b++) // the fast way
            for (auto p=blocks[b]; p<blocks[b] + ((b!=usedBlocks-1) ? blockSize : (size % blockSize)) ; p++)
                kh.add(*p);
        return kh.total();
    }
    inline T kahanSum(size_t from, size_t to) const { //
        Kahan kh;
        
        // first block
        for (auto pb=blocks[i2b(from)]+i2i(from), p=pb; p<pb+(blockSize - i2i(from)); p++) kh.add(*p);
        for (auto b=i2b(from)+1; b<i2b(to); b++) // rest of complete blocks
            for (auto pb=blocks[b], p=pb; p<pb + blockSize ; p++) kh.add(*p);
        for (auto pb=blocks[i2b(to)], p=pb; p<pb+i2i(to); p++) kh.add(*p); // last block
        return kh.total();
    }
    
    
    inline void applyFunc(size_t from, size_t to, Lambda lambda) const { // in 3 stages: 1st block, middle blocks, last block.
        // first block
        for (auto pb=blocks[i2b(from)]+i2i(from), p=pb; p<pb+(blockSize - i2i(from)); p++) *p=lambda(*p);
        for (auto b=i2b(from)+1; b<i2b(to); b++) // middle, rest of complete blocks
            for (auto pb=blocks[b], p=pb; p<pb + blockSize ; p++) *p=lambda(*p);
        for (auto pb=blocks[i2b(to)], p=pb; p<pb+i2i(to); p++) *p=lambda(*p); // last block
    }
    
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
        else                Thread(this).fill(c);
        return *this;
    }
    DynVect fill(T c) {  return DynVect(*this)._fill(c);    } // much slower than _mutable version
    
    DynVect generate(Lambda lambda) {
        DynVect v(*this);
        generate(v.begin(), v.end(), lambda);
        return v;
    }
    
    T stsum() const {
        return kahanSum();
    }
    
    void debug(size_t n=20) {
        for (auto i=0; i<n; i++)
            printf("%.2f, ", (*this)[i]);
    }
    
    class Iterator : public iterator<random_access_iterator_tag, T, ptrdiff_t, T*, T&> {
    public:
        
        DynVect*v=nullptr;
        long index=0;
        
        Iterator() = default;
        Iterator(DynVect*v) : v(v) {}
        
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
    
    
    class Thread {
    public:
        static const int maxThreads=256;
        int nthreads = thread::hardware_concurrency()-1; // std::min<int> ( thread::hardware_concurrency()-1, maxThreads );
        const DynVect*v=nullptr;
        mutex mutex;
        atomic<bool>sw; //=true;
        atomic<size_t> position; // =-1;
        size_t segsize=0;
        thread ths[maxThreads];
        
        // in a from, to scenario cal from, to value based in current thread 't' and segsize=size/nthreads
        inline size_t begin(int t)   const { return t*segsize; }
        inline size_t end(int t)     const { return (t==nthreads-1) ? v->size : (t+1)*segsize; }
        
        
        Thread(const DynVect*v) : v(v), segsize(v->size / nthreads), sw(true), position(-1) {}
        
        Thread(const DynVect*v, LambdaRange lambda) : v(v), segsize(v->size / nthreads), sw(true), position(-1) {
            for (auto t=0; t<nthreads; t++)
                ths[t] = thread ( [this, t, lambda] { lambda(begin(t), end(t)); });
            joinAll();
        }
        Thread(const DynVect*v, Lambda lambda) : v(v), segsize(v->size / nthreads), sw(true), position(-1) {
            for (auto t=0; t<nthreads; t++)
                ths[t] = thread ( [this, v, t, lambda] {
                    auto vv=(DynVect*)v; // override const v
                    for (auto i=begin(t); i < end(t); i++) {
                        T &d=(*vv)[i];
                        d=lambda(d);
                    }
                });
            joinAll();
        }
        Thread(const DynVect*v, LambdaNoArg lambda) : v(v), segsize(v->size / nthreads), sw(true), position(-1) {
            for (auto t=0; t<nthreads; t++)
                ths[t] = thread ( [this, v, t, lambda] {
                    auto vv=(DynVect*)v; // override const v
                    for (auto i=begin(t); i < end(t); i++) {
                        T &d=(*vv)[i];
                        d=lambda();
                    }
                });
            joinAll();
        }
        Thread(const DynVect*v, DynVect&other, LambdaBool lambda) : v(v), segsize(v->size / nthreads), sw(true), position(-1) {
            for (auto t=0; t<nthreads; t++)
                ths[t] = thread ( [this, t, lambda, v, &other] {
                    auto vl=DynVect();
                    for (auto i=begin(t); i < end(t); i++) {
                        T d=(*v)[i];
                        if(lambda(d)) vl << d;
                    }
                    // lock append th partial result
                    mutex.lock();  other << vl;   mutex.unlock();
                });
            joinAll();
        }
        void joinAll() {
            for (auto t=0; t<nthreads; t++) ths[t].join();
        }
        
        inline T kahanSum(int t) const {
            return v->kahanSum(begin(t), end(t));
        }
        T ret(T val) { return val; }
        const DynVect *vect() { return v; }
        
        T sumthread(int t=0) { // elegant & efficient recursive/multithreaded sum
            Kahan k;
            if (t<nthreads) {
                auto th=thread([&](){   k.add( sumthread(t+1) );   });
                k.add( kahanSum(t) );
                th.join();
            }
            return k.sum();
        }
        T sum(int t=0) { // async version, same performance as thread
            Kahan k;
            if (t<nthreads) {
                auto th=async( launch::async, [this, t] { return sum(t+1); } );
                k.add( kahanSum(t) + th.get() );
            }
            return k.total();
        }
        
        
        bool eq(const DynVect&other) {
            sw=true; // eq switch
            
            for (auto t=0; t<nthreads; t++)
                ths[t] = thread ( [this, t, other]()  {
                    for (auto i=begin(t); i < end(t) && sw; i++)
                        if ((*v)[i] != other[i]) { sw=false; break; }
                });
            joinAll();
            return sw;
        }
        
        size_t find(const T c) {
            sw=false; // found switch
            position=-1; // position
            
            for (auto t=0; t<nthreads; t++)
                ths[t] = thread([this, t, c]()   {
                    for (auto i=begin(t); i < end(t) && !sw; i++)
                        if ((*v)[i] == c) { // found -> break this and all threads
                            position=i; sw=true;
                            break;
                        };
                });
            joinAll();
            return position;
        }
        
        void fill(T c, int t=0) {
            if (t<nthreads) {
                auto th=async( launch::async, [this, t, c] { fill(c, t+1); } );
                auto vv=(DynVect*)v;
                for (auto i=begin(t); i<end(t); i++) (*vv)[i] = c;
                th.get();
            }
        }
        
        DynVect evaluate(const DynVect &other, Operator op) { //  this->v op= other
            assert(v->size==other.size);
            
            for (int t=0; t<nthreads; t++)
                ths[t] = std::thread([this, &other, t, op]()   {
                    DynVect&vm=(DynVect&)*v; // create a mutable ref as 'v' is const
                    switch (op) {
                        case opADD: for (auto i=begin(t); i<end(t); i++)  vm[i]+=other[i];   break;
                        case opSUB: for (auto i=begin(t); i<end(t); i++)  vm[i]-=other[i];   break;
                        case opMUL: for (auto i=begin(t); i<end(t); i++)  vm[i]*=other[i];   break;
                        case opDIV: for (auto i=begin(t); i<end(t); i++)  vm[i]/=other[i];   break;
                        case opPOW: for (auto i=begin(t); i<end(t); i++)  vm[i]=pow((*v)[i], other[i]);   break;
                    }
                });
            joinAll();
            return *v;
        }
        DynVect evaluate(T c, Operator op) { //  this->v op= c
            for (int t=0; t<nthreads; t++)
                ths[t] = std::thread([this, c, t, op]()   {
                    DynVect&vm=(DynVect&)*v; // create a mutable ref as 'v' is const
                    switch (op) {
                        case opADD: for (auto i=begin(t); i<end(t); i++)  vm[i]+=c;   break;
                        case opSUB: for (auto i=begin(t); i<end(t); i++)  vm[i]-=c;   break;
                        case opMUL: for (auto i=begin(t); i<end(t); i++)  vm[i]*=c;   break;
                        case opDIV: for (auto i=begin(t); i<end(t); i++)  vm[i]/=c;   break;
                        case opPOW: for (auto i=begin(t); i<end(t); i++)  vm[i]=pow((*v)[i], c);   break;
                    }
                });
            joinAll();
            return *v;
        }
    };
    
    T sum() const {
        return (size<minsizeMT) ? kahanSum() : Thread(this).sum();
    }
    
    T parSum(Iterator beg, Iterator end) { // experimental async
        auto len = end - beg;
        if (len < 1024*32)
            return (T) accumulate(beg, end, (T)0);
        
        Iterator mid = beg + len/2;
        auto handle = async(launch::async,[=]()->T{ return parSum(mid, end); });
        
        return parSum(beg, mid) + handle.get();
    }
    
    int multiThreadSort(Iterator begin, Iterator end) {
        const int szPage=4*1024;
        auto const sz = end - begin;
        if (sz <= 1) return 0;
        
        auto pivot = begin + sz/2;
        auto const pivot_v = *pivot;
        
        swap(*pivot, *(end - 1));
        auto p = partition(begin, end, [&](const T& a) { return a < pivot_v; } );
        swap(*p, *(end - 1));
        
        if (sz > szPage) {
            auto handle = async(launch::async, [&]() {
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
        if(size>1) multiThreadSort( begin(), end() );
        return *this;
    }
    DynVect mtsort() { // non mutable experimental use w/caution
        return DynVect(*this)._mtsort();
    }
    
    DynVect&_sort() { // mutable _sort
        if(size>1) std::sort( begin(), end() );
        return *this;
    }
    DynVect sort() { // non mutable sort
        return DynVect(*this)._sort();
    }
    
    bool isSorted() {
        return is_sorted(begin(), end());
    }
    
    DynVect&_shuffle() {
        random_device rd;
        mt19937 g(rd());
        
        if(size>1)std::shuffle(begin(), end(), g);
        return *this;
    }
    DynVect shuffle() {
        return DynVect(*this)._shuffle();
    }
    DynVect&_reverse() {
        if(size>1)std::reverse(begin(), end());
        return *this;
    }
    DynVect reverse() {
        return DynVect(*this)._reverse();
    }
    size_t find(const T c) {
        if(size>1) {
            if (size<minsizeMT) {
                auto ri = std::find( begin(), end(), c);
                return ri == end() ? -1 : end()-ri;
            } else return Thread(this).find(c);
        } return -1;
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
        if(size>1) {
            if (size<minsizeMT) for (auto &d:*this) d=foo(d);
            else Thread(this, foo);
        }
        return *this;
    }
    DynVect&_apply(LambdaNoArg foo) {
        if(size>1) {
            if (size<minsizeMT) for (auto &d:*this) d=foo();
            else Thread(this, foo);
        }
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
    DynVect operate(const T&c, Operator op) const { // non mutable  res = v <op> c
        DynVect v(*this);
        
        if (size<minsizeMT) {
            switch (op) {
                case opADD: for (Index i=0; i<size; i++) v[i] = (*this)[i]+c; break;
                case opSUB: for (Index i=0; i<size; i++) v[i] = (*this)[i]-c; break;
                case opMUL: for (Index i=0; i<size; i++) v[i] = (*this)[i]*c; break;
                case opDIV: for (Index i=0; i<size; i++) v[i] = (*this)[i]/c; break;
                case opPOW: for (Index i=0; i<size; i++) v[i] = pow((*this)[i],c); break;
            }
        } else Thread(&v).evaluate(c, op); // v=this <op> o
        return v;
    }
    DynVect operator+(const T c) const { return operate(c, opADD); }
    DynVect operator-(const T c) const { return operate(c, opSUB); }
    DynVect operator*(const T c) const { return operate(c, opMUL); }
    DynVect operator/(const T c) const { return operate(c, opDIV); }
    DynVect operator^(const T c) const { return operate(c, opPOW); }
    
    DynVect operateMutable(const T&c, Operator op) const { // mutable  v <op>= c
        if (size<minsizeMT) {
            switch (op) {
                case opADD: for (Index i=0; i<size; i++) (*this)[i]+=c; break;
                case opSUB: for (Index i=0; i<size; i++) (*this)[i]-=c; break;
                case opMUL: for (Index i=0; i<size; i++) (*this)[i]*=c; break;
                case opDIV: for (Index i=0; i<size; i++) (*this)[i]/=c; break;
                case opPOW: for (Index i=0; i<size; i++) (*this)[i]=pow((*this)[i],c); break;
            }
        } else Thread(this).evaluate(c, op); // v=this <op> o
        return *this;
    }
    DynVect& operator+=(const T c) {  for (Index i=0; i<size; i++) (*this)[i]+=c;    return *this; }
    DynVect& operator-=(const T c) {  for (Index i=0; i<size; i++) (*this)[i]-=c;    return *this; }
    DynVect& operator*=(const T c) {  for (Index i=0; i<size; i++) (*this)[i]*=c;    return *this; }
    DynVect& operator/=(const T c) {  for (Index i=0; i<size; i++) (*this)[i]/=c;    return *this; }
    
    // operator aritmetic, non mutable
    DynVect operate(const DynVect&other, Operator op) const { // res <op>= o
        assert(size == other.size);
        DynVect v(*this);
        
        if (size<minsizeMT) {
            switch (op) {
                case opADD: for (Index i=0; i<size; i++) v[i] = (*this)[i]+other[i]; break;
                case opSUB: for (Index i=0; i<size; i++) v[i] = (*this)[i]-other[i]; break;
                case opMUL: for (Index i=0; i<size; i++) v[i] = (*this)[i]*other[i]; break;
                case opDIV: for (Index i=0; i<size; i++) v[i] = (*this)[i]/other[i]; break;
                case opPOW: for (Index i=0; i<size; i++) v[i] = pow((*this)[i],other[i]); break;
            }
        } else Thread(&v).evaluate(other, op); // v=this <op> o
        return v;
    }
    DynVect operator+(const DynVect& other) const { return operate(other, opADD);  }
    DynVect operator-(const DynVect& other) const { return operate(other, opSUB);  }
    DynVect operator*(const DynVect& other) const { return operate(other, opMUL);  }
    DynVect operator/(const DynVect& other) const { return operate(other, opDIV);  }
    DynVect operator^(const DynVect& other) const { return operate(other, opPOW);  }
    
    DynVect& operateMutable(const DynVect&o, Operator op)  { // mutable ops: v += other;
        assert(size == o.size);
        
        if (size<minsizeMT) {
            switch (op) {
                case opADD: for (Index i=0; i<size; i++) (*this)[i]+=o[i]; break;
                case opSUB: for (Index i=0; i<size; i++) (*this)[i]-=o[i]; break;
                case opMUL: for (Index i=0; i<size; i++) (*this)[i]*=o[i]; break;
                case opDIV: for (Index i=0; i<size; i++) (*this)[i]/=o[i]; break;
                case opPOW: for (Index i=0; i<size; i++) (*this)[i]=pow((*this)[i],o[i]); break;
            }
        } else Thread(this).evaluate(o, op); // this=this <op> other
        return *this;
    }
    DynVect& operator+=(const DynVect& other)  { return operateMutable(other, opADD);  }
    DynVect& operator-=(const DynVect& other)  { return operateMutable(other, opSUB);  }
    DynVect& operator*=(const DynVect& other)  { return operateMutable(other, opMUL);  }
    DynVect& operator/=(const DynVect& other)  { return operateMutable(other, opDIV);  }
    DynVect& operator^=(const DynVect& other)  { return operateMutable(other, opPOW);  }
    
};
#endif /* DynVect_hpp */
