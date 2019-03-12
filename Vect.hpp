//
//  Vect.hpp
//  general purpose numerical vect class
//  debug array: p *(float(*)[10])data
// usage: see testVect() func

#ifndef Vect_hpp
#define Vect_hpp

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <string>
#include <functional>
#include <sstream>
#include <thread>


template <class T>
class Vect {
public:
    typedef Vect<size_t> VectIndex;
    typedef std::function<T(T)>const& Lambda;
    typedef enum { opADD, opSUB, opMUL, opDIV, opPOW } Operator;
    
    class Kahan {
    private:
        T s=0, c=0, inc=1, y, t;
    public:
        Kahan(T inc=1)  : inc(inc) {
            s = c = 0;
        }
        T next() {
            y = s+inc /*v[i]*/ - c;
            t = s + y;
            c = (t - s) - y;
            s = t;
            
            return s;
        }
    };
    
    class ThreadedCalc { //  multithreading support
    public:
        std::thread *thds=nullptr;
        int nthreads=0;
        size_t *ranges=nullptr;
        T *thValues=nullptr;
        T *maxv, *minv;
        Vect*vect=nullptr;
        
        ThreadedCalc(Vect&v) : vect(&v) {  init(v.size);  }
        ThreadedCalc(Vect&v, std::function<T()> const& lambda) : vect(&v) {
            init(v.size);
            apply(lambda);
        }
        Vect getVect() { return *vect; }
    private:
        void init(size_t size) {
            nthreads = std::thread::hardware_concurrency();
            thds = new std::thread[nthreads];
            
            auto segSize=size/nthreads;
            ranges=new size_t[nthreads+1];
            for (size_t t=0, rg=0; t<nthreads; t++, rg+=segSize) ranges[t]=rg;
            ranges[nthreads]=size;
            
            thValues=new T[nthreads];
            maxv=new T[nthreads];
            minv=new T[nthreads];
        }
        void joinAll() {
            for (int i=0; i<nthreads; i++)
                thds[i].join();
        }
    public:
        void apply(Lambda lambda) {
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t, lambda]() {
                    for (auto i=ranges[t]; i<ranges[t+1]; i++)
                        vect->data[i] = lambda(vect->data[i]);
                });
            }
            joinAll();
        }
        void apply(T from, T to, T inc, Lambda lambda) {
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t, from, inc, lambda]() {
                    for (auto i=ranges[t]; i<ranges[t+1]; i++)
                        vect->data[i] = lambda(from + inc * i);
                });
            }
            joinAll();
        }
        void apply( std::function<T()> const& lambda) {
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t, lambda]() {
                    for (auto i=ranges[t]; i<ranges[t+1]; i++)
                        vect->data[i] = lambda();
                });
            }
            joinAll();
        }
        
        T sum() {
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t]()   {
                    thValues[t]=kahanSum(vect->data, ranges[t], ranges[t+1]);
                });
            }
            joinAll();
            
            return kahanSum(thValues, 0, nthreads);
        }
        
        Vect evaluate(const Vect &other, Operator op) { //  ratio 1.2 improvement w/8 threads
            assert(vect->size==other.size);
            Vect v(other.size);
            
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, &v, &other, t, op]()   {
                    switch (op) {
                        case opADD: for (auto i=ranges[t]; i<ranges[t+1]; i++)  v[i]=vect->data[i] + other[i];   break;
                        case opSUB: for (auto i=ranges[t]; i<ranges[t+1]; i++)  v[i]=vect->data[i] - other[i];   break;
                        case opMUL: for (auto i=ranges[t]; i<ranges[t+1]; i++)  v[i]=vect->data[i] * other[i];   break;
                        case opDIV: for (auto i=ranges[t]; i<ranges[t+1]; i++)  v[i]=vect->data[i] / other[i];   break;
                        case opPOW: for (auto i=ranges[t]; i<ranges[t+1]; i++)  v[i]=pow(vect->data[i], other[i]);   break;
                    }
                });
            }
            joinAll();
            return v;
        }
        Vect evaluateMutable(const Vect &other, Operator op) { //  ratio 1.2 improvement w/8 threads
            assert(vect->size==other.size);
            
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, &other, t, op]()   {
                    switch (op) {
                        case opADD: for (auto i=ranges[t]; i<ranges[t+1]; i++)  vect->data[i] += other[i];   break;
                        case opSUB: for (auto i=ranges[t]; i<ranges[t+1]; i++)  vect->data[i] -= other[i];   break;
                        case opMUL: for (auto i=ranges[t]; i<ranges[t+1]; i++)  vect->data[i] *= other[i];   break;
                        case opDIV: for (auto i=ranges[t]; i<ranges[t+1]; i++)  vect->data[i] /= other[i];   break;
                        case opPOW: for (auto i=ranges[t]; i<ranges[t+1]; i++)  vect->data[i]=pow(vect->data[i], other[i]);   break;
                    }
                });
            }
            joinAll();
            return *vect;
        }
        Vect evaluate(const T c, Operator op) { //  ratio 1.2 improvement w/8 threads
            
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t, c, op]()   {
                    switch (op) {
                        case opADD: for (auto i=ranges[t]; i<ranges[t+1]; i++)  vect->data[i]+=c;   break;
                        case opSUB: for (auto i=ranges[t]; i<ranges[t+1]; i++)  vect->data[i]-=c;   break;
                        case opMUL: for (auto i=ranges[t]; i<ranges[t+1]; i++)  vect->data[i]*=c;   break;
                        case opDIV: for (auto i=ranges[t]; i<ranges[t+1]; i++)  vect->data[i]/=c;   break;
                        case opPOW: for (auto i=ranges[t]; i<ranges[t+1]; i++)  vect->data[i]=pow(vect->data[i],c);   break;
                    }
                });
            }
            joinAll();
            return *vect;
        }
        
        T min() { return minmax().first;  }
        T max() { return minmax().second; }
        std::pair<T,T>minmax() {
            assert(vect->size>0);
            volatile T resmax=*vect->data, resmin=*vect->data;
            
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t]()   {
                    volatile T max=vect->data[ranges[t]], min=max;
                    for (auto i=ranges[t]; i<ranges[t+1]; i++) {
                        if (vect->data[i] > max) max=vect->data[i];
                        if (vect->data[i] < min) min=vect->data[i];
                    }
                    maxv[t]=max;
                    minv[t]=min;
                });
            }
            joinAll();
            
            resmax = *std::max_element(maxv, maxv+nthreads);
            resmin = *std::min_element(minv, minv+nthreads);
            
            return std::pair<T,T>(resmin, resmax);
        }
        size_t locate(const T c) {
            assert(vect->size>0);
            volatile size_t res=-1;
            volatile bool found=false;
            
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, &res, t, c, &found]()   {
                    for (auto i=ranges[t]; i<ranges[t+1] && !found; i++)
                        if (vect->data[i] == c) { res=i; found=true; break; };
                });
            }
            joinAll();
            return res;
        }
        void seq() {
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t]() {
                    for (auto i=ranges[t]; i<ranges[t+1]; i++)
                        vect->data[i] = (T)i;
                });
            }
            joinAll();
        }
        void seq(T inc) {
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t, inc]() {
                    for (auto i=ranges[t]; i<ranges[t+1]; i++)
                        vect->data[i] = inc * i;
                });
            }
            joinAll();
        }
       
        void seq(T from, T to, T inc) {
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t, from, inc]() {
                    for (auto i=ranges[t]; i<ranges[t+1]; i++)
                        vect->data[i] = from + inc*i;
                });
            }
            joinAll();
        }
        
        ~ThreadedCalc() {
            if(ranges!=nullptr) {
                delete[] thds;
                delete[] ranges;
                delete[] thValues;
                delete[] maxv;
                delete[] minv;
                ranges=nullptr;
            }
        }
    };
    
private:
    size_t realSize=0;
    const int nchunk=128; //  memory increments on each resize
    
    T *data=nullptr;
    size_t size=0;
    size_t _index=0; // used in iterator
    
    static T kahanSum(const T *v, size_t from, size_t to) { // sum vector range usign kahan algorithm
        T s=0, c=0;
        for (auto i=from; i<to; i++) {
            auto y = v[i] - c, t = s + y;
            c = (t - s) - y;
            s = t;
        }
        return s;
    }
    static size_t sizeOfRange(const T from, const T to, const T inc) {
        return (size_t)ceil((to-from)/inc);
    }
    
public:
    
    // constructors
    Vect() { }
    
    Vect(const size_t size) {
        alloc(size);
    }
    Vect(const Vect &other) { // Vect v0(v1), v2=v1;
        if(other.size>0) {
            alloc(other.size);
            memcpy(data, other.data, other.size * sizeof(T));
        }
    }
    Vect(const size_t size, const T v[]){
        alloc(size);
        for (auto &d:*this) d=v[_index++];
    }
    Vect(const T from, const T to, const T inc) {
        alloc(sizeOfRange(from, to, inc));
        ThreadedCalc(*this).seq(from, to, inc);
    }
    Vect(T from, T to, T inc, Lambda lambda) {
        assert(from < to && inc>0);
        alloc( sizeOfRange(from, to, inc) );
        ThreadedCalc(*this).apply(from, to, inc, lambda);
    }
    Vect(const Vect &v, Lambda lambda) { // Vect<double>vx(0, M_PI, 0.001), vy(vx, sin);
        if(v.size>0) {
            *this = v; // make a copy & apply lambda to copy (*this)
            ThreadedCalc(*this).apply(lambda);
        }
    }
    
    ~Vect() { dealloc(); }
    
    static Vect rnd(int size) { // norm random
        assert(size>0);
        Vect v(size);
        ThreadedCalc(v, []() -> T{ return (T)rand()/(T)RAND_MAX; });
        return v;
    }
    static Vect seq(size_t size) { // 0,1...n-1
        assert(size>0);
        Vect v(size);
        ThreadedCalc(v).seq();
        return v;
    }
    static Vect seq(size_t size, T inc) { // 0,inc, 2*inc, ...(n-1)*inc
        assert(size>0);
        Vect v(size);
        ThreadedCalc(v).seq(inc);
        return v;
    }
   
    static Vect stSeq(T from, T to, T inc) {
        assert((from<to && inc>0) && "seq paramaters error");
        Vect v(sizeOfRange(from, to, inc));
        for (size_t i=0; i<v.length(); i++) v.data[i]=(T)i * inc;
        return v;
    }
    
    
    static Vect genWave(T secs, int sampleRate, int nw, T amp[], T hz[], T phase[]) {
        Vect v;
        for (T t=0; t<secs*sampleRate; t++) {
            double res=0, freq2inc = 2.*M_PI/sampleRate * t;
            for (int i=0; i<nw; i++)
                res += amp[i] * sin( hz[i]*freq2inc + phase[i] );
            v << (T)(res/nw);
        }
        return v;
    }
    
private:
    // memory mgr.
    void dealloc() {
        if (data!=nullptr && size!=0) {
            free(data);
            data=nullptr;
            size=realSize=0;
        }
    }
    void alloc(size_t size) {
        assert(((data = (T*)calloc(size, sizeof(T))) != nullptr ) && "calloc error, memory full?");
        realSize = this->size = size;
    }
    void resize(size_t size) {
        assert (((data = (T*)realloc(data, size*sizeof(T))) != nullptr ) && "realloc error, memory full?");
        if(size>this->size)
            memset(this->data, 0, (size-this->size) * sizeof(T));
        realSize = this->size = size;
    }
    
public:
    std::string toString()  {
        std::ostringstream oss;
        for (auto const d:*this) oss << d << ",";
        return oss.str();
    }
    size_t length() {   return size;    }
    size_t count() {    return size;    }
    
    Vect sub(size_t from, size_t to) { // subvector from,to
        assert(!(from>to || from>=size || to>=size) && "sub: from > to or range error");
        auto sz=to-from;
        Vect res(sz);
        for (auto i=from; i<to; i++)
            res.data[i-from]=data[i];
        return res;
    }
    Vect left(size_t n) {
        return sub(0, n);
    }
    Vect right(size_t n) {
        return sub(size-n, size);
    }
    void random() {
        ThreadedCalc(*this, []() -> T {return (T)rand();});
    }
    void rand01() {
        ThreadedCalc(*this, []() -> T { return (T)rand()/(T)RAND_MAX; });
    }
    T seqSum() { // don't use as it generates precission errors -> apply kahan algo.
        T s=0;
        for (auto d:*this) s+=d;
        return s;
    }
    std::pair<T,T>seqMinmax() {
        T mx=*data, mn=*data;
        for (size_t i=0; i<size; i++) {
            if(data[i]>mx) mx=data[i];
            if(data[i]<mn) mn=data[i];
        }
        return std::pair<T,T>(mn,mx);
    }
    std::tuple<T,T,T>seqMinmaxdiff() {
        auto mm=seqMinmax();
        return std::tuple<T,T,T>(mm.first,mm.second,mm.second-mm.first);
    }
    T seqMin() {
        T mn=*data;
        for (size_t i=0; i<size; i++)
            if(data[i]<mn) mn=data[i];
        return mn;
    }
    T seqMax() {
        T mx=*data;
        for (size_t i=0; i<size; i++)
            if(data[i]>mx) mx=data[i];
        return mx;
    }
    
    std::pair<T,T>minmax() {
        return ThreadedCalc(*this).minmax();
    }
    std::tuple<T,T,T>minmaxdiff() {
        auto mm=minmax();
        return std::tuple<T,T,T>(mm.first,mm.second,abs(mm.second-mm.first));
    }
   
    T min() {   return ThreadedCalc(*this).min();    }
    T max() {   return ThreadedCalc(*this).max();    }
    T avg()   {   return sum()/(T)size;  }
    T diff() {
        auto mm=minmax();
        return abs( mm.first - mm.second );
    }
    Vect scale(T from=0, T to=1) {
        auto mmd=minmaxdiff();
        auto mn=std::get<0>(mmd), df=std::get<2>(mmd);
        auto pdf=to-from;
        
        if(df)
            for (auto &d:*this)
                d = from + ((d-mn) / df) * pdf;
        return *this;
    }
    Vect norm() {
        auto mmd=minmaxdiff(); 
        auto mn=std::get<0>(mmd), df=std::get<2>(mmd);
        if(df)
            for (auto &d:*this)
                d = (d-mn) / df;
        return *this;
    }
    Vect fill(T value) {
        for (auto &d:*this) d = value;
        return *this;
    }
    Vect sequence(T from, T to) {
        assert(from<to && "sequence error: from > to");
        if(size) {
            T inc=(T)(to-from)/(T)size;
            for (auto &d:*this) d = inc * (_index++);
        }
        return *this;
    }
    Vect reverse() {
        for (size_t i=0, j=size-1; i<size/2; i++,j--)
            std::swap<T>(data[i], data[j]);
        return *this;
    }
    T polyEval(T x) {
        T r=0, ex=1;
        for (size_t i=0; i<size; i++, ex*=x)
            r+=data[i] * ex;
        return r;
    }
    
    Vect sort(bool increasing=true) {
        if (increasing)
            std::sort(begin(), end(), [](T a, T b){ return a<b; });
        else
            std::sort(begin(), end(), [](T a, T b){ return b<a; });
        return *this;
    }
    Vect sortInv() {
        std::sort(begin(), end(), [](T a, T b){ return b<a; });
        return *this;
    }
    Vect shuffle() {
        if (size > 1) {
            for (int i=0; i<size; i++)
                std::swap<T>(data[i], data[(rand() % (i + 1))]);
        }
        return *this;
    }
    T stSum() { // Single Thread Kahan summation algorithm, sequential sum generates precission errors
        return kahanSum(data, 0, size);
    }
    
    Vect apply(Lambda lambda) { // apply labda func -> mutable, usage: v.apply(sin)
        for (auto &d:*this) d=lambda(d);
        return *this;
    }
    Vect apply(std::function<T()> const& lambda) { // usage: v.apply(funcNoArgs)
        for (auto &d:*this) d=lambda();
        return *this;
    }
    Vect seqFunc(Lambda lambda) { // non mutable -> return a copy applying func
        Vect v(*this);
        for (auto &d:v) d=lambda(d);
        return v;
    }
    Vect func(std::function<T()> const& lambda) { // usage: func(funcNoArgs)
        Vect v(*this);
        for (auto &d:v) d=lambda();
        return v;
    }
    Vect func(Lambda lambda) { // non mutable -> aply on copy
        Vect v(*this);
        ThreadedCalc(v).apply(lambda);
        return v;
    }
   
    T sum() {  return ThreadedCalc(*this).sum();  }
    
    Vect filter(std::function<bool(T)> const& lambda) { // filter by bool lambda conditional expr.
        Vect v;
        for (auto const d:*this) if (lambda(d)) v<<d;
        return v;
    }
    VectIndex filterIndex(std::function<bool(T)> const& lambda) { // indexes of selected items
        VectIndex res;
        for (size_t ix=0; ix<size; ix++)
            if (lambda(data[ix])) res<<ix;
        return res;
    }
    
    T*begin() { _index=0; return data; }
    T*end() { return data+size; }
    size_t incIndex() { return _index++; } // for (auto d:vect) v1[vect.incIndex()]=d;
    
    // append, void return as Vect return -> very bad performnce
    // so v0 << v1 << v2 is v0<<v1; v0<<v2;
    void operator << (const T c) {
        append(c);
    }
    void operator >> (T &c) {
        assert(_index<size);
        c=data[_index++];
    }
    void operator << (const Vect &other) {
        return append(other);
    }
    void append(const T c) {
        if ( size >= realSize ) {
            realSize+=nchunk;
            data = (T*)realloc(data, realSize*sizeof(T));
        }
        data[size++]=c;
    }
    void append(const Vect &other) {
        auto os=other.size;
        if ( size+os >= realSize ) {
            realSize=size+os;
            data = (T*)realloc(data, realSize*sizeof(T));
        }
        for (int i=0; i<os; i++)
            data[i+size]=other.data[i];
        size+=os;
    }
    void clear() {
        dealloc();
        size=realSize=0;
    }
    Vect zero() {
        memset(data, 0, sizeof(T) * size);
        return *this;
    }
    T erase(size_t ix) { // return removed item
        assert(ix>=0 && ix<size);
        T c=data[ix];
        if (ix<=size-1)
            memmove(data+ix, data+ix+1, sizeof(T) * (size-ix-1));
        size--;
        return c;
    }
    void fit() { // fit current memory to size
        if (size)
            resize(size);
    }
    size_t locate(T c) {
        for (size_t i=0; i<size; i++)
            if (c==data[i]) return i;
        return -1;
    }
    size_t mtLocate(T c) {
        return ThreadedCalc(*this).locate(c);
    }
    size_t bsearch(T c) { //  binary search -> must be sorted
        auto ret = std::bsearch(&c, begin(), size, sizeof(T),
                                [](const void*a, const void*b) ->int { return (int)(*((T*)a) - *((T*)b));});
        if(ret==nullptr) return -1;
        return (T*)ret-begin();
    }
    
    
    Vect &operator=(const Vect &other) { // asignment
        // check for self-assignment
        if (&other == this) return *this;
        
        if (size != other.size) resize(other.size);
        std::copy(&other.data[0], &other.data[0] + size, &data[0]);
        return *this;
    }
    Vect &operator=(const T c) { // asignment
        for (auto &d:*this) d=c;
        return *this;
    }
    
    
    // index mutator
    inline T& operator[](size_t index) {
        assert(!(index<0 || index>=size) && "Vect: index error");
        return data[index];
    }
    // index accessor
    inline const T& operator[](size_t index) const {
        assert(!(index<0 || index>=size) && "Vect: index error");
        return data[index];
    };
    inline const T& operator[](const Vect&v) const { // vect[v] is vect[v._index]
        assert(!(v._index<0 || v._index>=size) && "Vect: index error");
        return data[v._index];
    };
    Vect operator[](VectIndex ixs) const { // index vect by vector of indexes
        Vect v;
        for (auto ix:ixs) v << data[ix];
        return v;
    }
    Vect operator[](std::function<bool(T)> const& lambda) { // index vect by vector of indexes
        Vect v;
        for (auto d : *this)
            if (lambda(d)) v << d;
        return v;
    }
    
    // boolean,
    bool operator==(const Vect &other) { // vect==vect
        bool eq = (size==other.size);
        size_t itemNE=0;
        if(eq)
            for (size_t i = 0; i < size; i++)
                if (data[i] != other.data[i]) { eq=false; itemNE=i; break; }
        return eq;
    }
    bool operator==(const T c) { // vect == scalar, true if all ==
        bool eq=true;
        for (auto d:*this) if(d!=c) {eq=false; break;}
        return eq;
    }
    bool operator!=(const Vect &other) {
        bool eq = (size==other.size);
        if(eq) {
            for (size_t i = 0; i < size; i++)
                if (data[i] != other.data[i]) break;
        }
        return eq;
    }
    bool operator!=(const T c) { // vect != scalar
        bool res=true;
        for (auto d:*this) if(d==c) {res=false; break;}
        return res;
    }
    bool operator > (const T c) { // vect > scalar, all > scalar
        bool res=true;
        for (auto d:*this) if(d<=c) {res=false; break;}
        return res;
    }
    bool operator < (const T c) { // vect < scalar, all < scalar
        bool res=true;
        for (auto d:*this) if(d>=c) {res=false; break;}
        return res;
    }
    // compares smallest vectors
    bool operator>(const Vect &other) {
        bool gt=true;
        for (size_t i = 0; i < std::min(size, other.size); i++)
            if (data[i] <= other.data[i]) { gt=false; break; }
        return gt;
    }
    bool operator<(const Vect &other) {
        bool lt=true;
        for (size_t i = 0; i < std::min(size, other.size); i++)
            if (data[i] >= other.data[i]) { lt=false; break; }
        return lt;
    }
    
    
    // Vect op constant
    Vect operator+(const T &c) {
        Vect v(*this);
        return ThreadedCalc(v).evaluate(c, opADD);
    }
    Vect operator-(const T &c) {
        Vect v(*this);
        return ThreadedCalc(v).evaluate(c, opSUB);
    }
    Vect operator*(const T &c) {
        Vect v(*this);
        return ThreadedCalc(v).evaluate(c, opMUL);
    }
    Vect operator/(const T &c) {
        Vect v(*this);
        return ThreadedCalc(v).evaluate(c, opDIV);
    }
    Vect operator^(const T &c) {
        Vect v(*this);
        return ThreadedCalc(v).evaluate(c, opPOW);
    }
    
    // Vect (+-*/)= const<T>
    
    Vect &operator+=(const T &c) {
        if(c!=0) ThreadedCalc(*this).evaluate(c, opADD);
        return *this;
    }
    Vect &operator-=(const T &c) {
        if(c!=0) ThreadedCalc(*this).evaluate(c, opSUB);
        return *this;
    }
    Vect &operator*=(const T &c) {
        if(c!=0) ThreadedCalc(*this).evaluate(c, opMUL);
        return *this;
    }
    Vect &operator/=(const T &c) {
        if(c!=0) ThreadedCalc(*this).evaluate(c, opDIV);
        return *this;
    }
    
    Vect operator|(const Vect &other) { // union
    }
    Vect operator&(const Vect &other) { // intersection
    }
    Vect operator||(const Vect &other) { //
    }
    Vect operator&&(const Vect &other) {
    }
    
    
    // Vect op(+-*/) Vect
    
    Vect operator+(const Vect &other) {
        assert(size==other.size);
        Vect v(*this);
        return ThreadedCalc(v).evaluate(other, opADD);
    }
    Vect operator-(const Vect &other) {
        assert(size==other.size);
        Vect v(*this);
        return ThreadedCalc(v).evaluate(other, opSUB);
    }
    Vect operator*(const Vect &other) {
        assert(size==other.size);
        Vect v(*this);
        return ThreadedCalc(v).evaluate(other, opMUL);
    }
    Vect operator/(const Vect &other) {
        assert(size==other.size);
        Vect v(*this);
        return ThreadedCalc(v).evaluate(other, opDIV);
    }
    
    // Vect (+-*/)= Vect
    
    Vect &operator+=(const Vect &other) {
        for (size_t i = 0; i < std::min(size, other.size); i++)
            data[i] += other.data[i];
        return *this;
    }
    Vect &operator-=(const Vect &other) {
        for (size_t i = 0; i < std::min(size, other.size); i++)
            data[i] -= other.data[i];
        return *this;
    }
    Vect &operator*=(const Vect &other) {
        for (size_t i = 0; i < std::min(size, other.size); i++)
            data[i] *= other.data[i];
        return *this;
    }
    Vect &operator/=(const Vect &other) {
        for (size_t i = 0; i < std::min(size, other.size); i++)
            data[i] /= other.data[i];
        return *this;
    }
    
    // ++, --
    size_t &operator++() {
        return ++_index;
    }
    size_t &operator--() {
        return --_index;
    }
    // postfix: &operator++(int) implements _index++ -> v++
    size_t operator++(int) {
        return _index++;
    }
    size_t operator--(int) {
        return _index--;
    }
};


void testVect();

#endif /* Vect_hpp */
