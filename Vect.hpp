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
#include <typeinfo>


template <class T>
class Vect {
public:
    typedef Vect<size_t> VectIndex;
    typedef std::function<T(T)>const& Lambda;
    typedef enum { opADD, opSUB, opMUL, opDIV, opPOW } Operator;
    
    static const size_t szSingle=10000; // below use single thread, over multithread
    
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
        
        bool found=false;
        size_t res=-1;
        
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
                    for (auto d=vect->data+ranges[t], end=vect->data+ranges[t+1]; d<end; d++)
                        *d = lambda(*d);
                });
            }
            joinAll();
        }
        void apply(T from, T to, T inc, Lambda lambda) {
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t, from, inc, lambda]() {
                    int i=0;
                    for (auto d=vect->data+ranges[t], end=vect->data+ranges[t+1]; d<end; d++)
                        *d = lambda(from + inc * (i++));
                });
            }
            joinAll();
        }
        void apply( std::function<T()> const& lambda) {
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t, lambda]() {
                    for (auto d=vect->data+ranges[t], end=vect->data+ranges[t+1]; d<end; d++)
                        *d = lambda();
                });
            }
            joinAll();
        }
        void scale( T from, T to, T min, T diff) {
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t, from, to, min, diff]() {
                    auto pdf=(to-from);
                    for (auto d=vect->data+ranges[t], end=vect->data+ranges[t+1]; d<end; d++)
                        *d = from + ((*d-min) / diff) * pdf;
                });
            }
            joinAll();
        }
        Vect &norm(T min, T diff) {
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t, min, diff]() {
                    for (auto d=vect->data+ranges[t], end=vect->data+ranges[t+1]; d<end; d++)
                        *d = (*d-min) / diff;
                });
            }
            joinAll();
            return *vect;
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
        Vect& rnd() {
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t]() {
                    unsigned int seedp;
                    for (auto d=vect->data+ranges[t], end=vect->data+ranges[t+1]; d<end; d++)
                        *d = (T)rand_r(&seedp) / (T)RAND_MAX; // use reentrant version of rand
                });
            }
            joinAll();
            return *vect;
        }
        Vect& filter(Vect&vres, std::function<bool(T)> const& lambda) {
            Vect *vs=new Vect[nthreads];
            
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t, lambda, &vs]() {
                    for (auto d=vect->data+ranges[t], end=vect->data+ranges[t+1]; d<end; d++)
                        if(lambda(*d)) vs[t]<<*d;
                });
            }
            joinAll();
            
            // join all vs[] in vect
            for (int t=0; t<nthreads; t++) vres << vs[t];
            delete[]vs;
            
            return vres;
        }
        VectIndex& filterIndex(VectIndex& vres, std::function<bool(T)> const& lambda) {
            VectIndex *vs=new VectIndex[nthreads];
            
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t, lambda, &vs]() {
                    size_t i=ranges[t];
                    for (auto d=vect->data+ranges[t], end=vect->data+ranges[t+1]; d<end; d++, i++)
                        if(lambda(*d)) vs[t] << i;
                });
            }
            joinAll();
            
            // join all vs[] in vect
            for (int t=0; t<nthreads; t++) vres << vs[t];
            delete[]vs;
            
            return vres;
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
        Vect& evaluateMutable(const Vect &other, Operator op) { //  ratio 1.2 improvement w/8 threads
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
        Vect& evaluate(const T c, Operator op) { //  ratio 1.2 improvement w/8 threads
            
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
            found=false;  res=-1;
            
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t, c]()   {
                    for (auto i=ranges[t]; i<ranges[t+1] && !found; i++)
                        if (vect->data[i] == c) { res=i; found=true; break; };
                });
            }
            joinAll();
            return res;
        }
        bool eq(const Vect &other) {
            found = vect->size == other.size; // read found as 'eq'
            
            if(found) {
                for (int t=0; t<nthreads; t++) {
                    thds[t] = std::thread([this, t, &other]()   {
                        for (auto i=ranges[t]; i<ranges[t+1] && found; i++)
                            if (vect->data[i] != other.data[i]) { found=false; break; };
                    });
                }
                joinAll();
            }
            return found; // eq
        }
        bool eq(const T c) {
            found = true; // read found as 'eq'
            
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t, c]()   {
                    for (auto i=ranges[t]; i<ranges[t+1] && found; i++)
                        if (vect->data[i] != c) { found=false; break; };
                });
            }
            joinAll();
            return found; // eq
        }
        Vect& seq() {
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t]() {
                    for (auto i=ranges[t]; i<ranges[t+1]; i++)
                        vect->data[i] = (T)i;
                });
            }
            joinAll();
            return *vect;
        }
        Vect& seq(T inc) {
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t, inc]() {
                    for (auto i=ranges[t]; i<ranges[t+1]; i++)
                        vect->data[i] = inc * i;
                });
            }
            joinAll();
            return *vect;
        }
       
        Vect& seq(T from, T to, T inc) {
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t, from, inc]() {
                    for (auto i=ranges[t]; i<ranges[t+1]; i++)
                        vect->data[i] = from + inc*i;
                });
            }
            joinAll();
            return *vect;
        }
        
        Vect& shuffle() {
            for (int t=0; t<nthreads; t++) {
                thds[t] = std::thread([this, t]() {
                    auto s=0u;
                    for (auto i=ranges[t], from=i, diff=ranges[t+1]-ranges[t]; i<ranges[t+1]; i++)
                        std::swap<T>(vect->data[i], vect->data[from + (rand_r(&s) % diff)]);
                });
            }
            joinAll();
            return *vect;
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
    Vect (const std :: initializer_list<T> inputs) {
        for (auto &i:inputs) *this << i;
    }
    Vect(const T from, const T to, const T inc) {
        assert(from<to && inc!=0);
        auto sz=sizeOfRange(from, to, inc);
        alloc(sz);
        if (sz<szSingle)  for (size_t i=0; i<sz; i++) data[i]=from+(T)i * inc;
        else              ThreadedCalc(*this).seq(from, to, inc);
    }
    Vect(T from, T to, T inc, Lambda lambda) {
        assert(from < to && inc>0);
        auto sz=sizeOfRange(from, to, inc);
        alloc(sz);
        if (sz<szSingle)  for (size_t i=0; i<sz; i++) data[i]=lambda((T)i * inc);
        else              ThreadedCalc(*this).apply(from, to, inc, lambda);
    }
    Vect(const Vect &v, Lambda lambda) { // Vect<double>vx(0, M_PI, 0.001), vy(vx, sin);
        if(v.size>0) {
            *this = v; // make a copy & apply lambda to copy (*this)
            if (size<szSingle) for (auto &d:*this) d=lambda(d);
            else ThreadedCalc(*this).apply(lambda);
        }
    }
    
    ~Vect() { dealloc(); }
    
    static Vect rnd(size_t size) { // norm random
        assert(size>0);
        Vect v(size);
        if (v.size<szSingle)  for (auto &d:v) d=(T)rand()/(T)RAND_MAX;
        else                  ThreadedCalc(v).rnd(); 
        return v;
    }
    static Vect strnd(size_t size) { // norm random
        assert(size>0);
        Vect v(size);
        for (auto &d:v) d=(T)rand()/(T)RAND_MAX;
        return v;
    }
    static Vect seq(size_t size) { // 0,1...n-1
        assert(size>0);
        Vect v(size);
        if (v.size<szSingle)  for (auto &d:v) d=v++; // v++ is v.index++;
        else                  ThreadedCalc(v).seq();
        return v;
    }
    static Vect seq(size_t size, T inc) { // 0,inc, 2*inc, ...(n-1)*inc
        assert(size>0);
        Vect v(size);
        if (v.size<szSingle)  for (auto &d:v) d=inc * (T)v++;
        else                  ThreadedCalc(v).seq(inc);
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
    Vect&csv(std::string s, char delimiter=',') {
        size_t pos_ant=0;
        for (size_t pos=0; pos!=std::string::npos; pos_ant=pos ? pos+1:pos, pos=s.find(delimiter, pos+1))
            if(pos)
                *this << s.substr(pos_ant, pos-pos_ant);
        return *this << s.substr(pos_ant, s.length()-pos_ant);
    }
    std::string toString()  {
        std::ostringstream oss;
        for (auto const d:*this) oss << d << ",";
        return oss.str();
    }
    std::string debug(size_t n=10) {
        std::ostringstream oss;
        for (auto i=0; i<n; i++) oss << data[i] << ",";
        return oss.str();
    }
    size_t length() {   return size;    }
    size_t count() {    return size;    }
    
    Vect sub(size_t from, size_t to) { // subvector from,to
        assert(from<to && from<size && to<=size && "sub: from > to or range error");
        auto sz=to-from;
        Vect res(sz);
        memcpy(res.data, data+from, sz * sizeof(T)); // for (auto i=from; i<to; i++) res.data[i-from]=data[i];
        return res;
    }
    Vect left(size_t n) {
        return sub(0, n);
    }
    Vect right(size_t n) {
        return sub(size-n, size);
    }
    void random() {
        if (size<szSingle)  for (auto &d:*this) d=(T)rand();
        else                ThreadedCalc(*this).rnd(); //, []() -> T {return (T)rand_r(seedp);});
    }
    
    std::pair<T,T>minmax() {
        return ThreadedCalc(*this).minmax();
    }
    std::tuple<T,T,T>minmaxdiff() {
        auto mm=minmax();
        return std::tuple<T,T,T>(mm.first,mm.second,abs(mm.second-mm.first));
    }
   
    T min() { return (size<szSingle) ? stmin() : ThreadedCalc(*this).min();    }
    T max() { return (size<szSingle) ? stmax() : ThreadedCalc(*this).max();    }
    T avg() { return sum()/(T)size;  }
    T diff() {
        auto mm=(size<szSingle) ? stminmax() : minmax();
        return abs( mm.first - mm.second );
    }
    Vect scale(T from=0, T to=1) {
        assert(from<to);
        auto mmd=minmaxdiff();
        auto min=std::get<0>(mmd), diff=std::get<2>(mmd);
        
        return (diff) ? ((size<szSingle) ? stscale(from, to, min,diff) : ThreadedCalc(*this).scale(from, to, min, diff)) : *this;
    }
    Vect norm() {
        auto mmd=minmaxdiff(); 
        auto min=std::get<0>(mmd), diff=std::get<2>(mmd);
        return (diff) ? ((size<szSingle) ? stnorm(min, diff) : ThreadedCalc(*this).norm(min, diff)) : *this;
    }
    Vect& fill(T value) {
        std::fill(begin(), end(), value);
        return *this;
    }
    
    Vect& reverse() {
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
    
    Vect& sort(bool increasing=true) {
        if (increasing)
            std::sort(begin(), end(), [](T a, T b){ return a<b; });
        else
            std::sort(begin(), end(), [](T a, T b){ return b<a; });
        return *this;
    }
    Vect& sortInv() {
        std::sort(begin(), end(), [](T a, T b){ return b<a; });
        return *this;
    }
    Vect& shuffle() {
        if (size > 1) {
            if (size<szSingle) {
                for (int i=0; i<size; i++)
                    std::swap<T>(data[i], data[(rand() % (i + 1))]);
            } else
                ThreadedCalc(*this).shuffle();
        }
        return *this;
    }
    
    // apply mutable, func non mutable
    Vect& apply(Lambda lambda) { // apply labda func -> mutable, usage: v.apply(sin)
        if (size<szSingle) for (auto &d:*this) d=lambda(d);
        else ThreadedCalc(*this).apply(lambda);
        return *this;
    }
    Vect& apply(std::function<T()> const& lambda) { // usage: v.apply(funcNoArgs)
        if (size<szSingle) for (auto &d:*this) d=lambda();
        else ThreadedCalc(*this).apply(lambda);
        return *this;
    }
    
    Vect func(std::function<T()> const& lambda) { // usage: func(funcNoArgs)
        Vect v(*this);
        if (size<szSingle) for (auto &d:v) d=lambda();
        else ThreadedCalc(v).apply(lambda);
        return v;
    }
    Vect func(Lambda lambda) { // non mutable -> aply on copy
        Vect v(*this);
        if (size<szSingle) for (auto &d:v) d=lambda(d);
        else ThreadedCalc(v).apply(lambda);
        return v;
    }
   
    T sum() {  return (size<szSingle) ? stsum() : ThreadedCalc(*this).sum();  }
    
    Vect filter(std::function<bool(T)> const& lambda) { // filter by bool lambda conditional expr.
        Vect v;
        if (size<szSingle) { for (auto const d:*this) if (lambda(d)) v<<d; }
        else ThreadedCalc(*this).filter(v, lambda);
        return v;
    }
   
    VectIndex filterIndex(std::function<bool(T)> const& lambda) { // indexes of selected items
        VectIndex res;
        if (size<szSingle) {  for (size_t ix=0; ix<size; ix++)
            if (lambda(data[ix])) res<<ix;
        } else ThreadedCalc(*this).filterIndex(res, lambda);
        return res;
    }
    
    // iterator support
    T*begin()   { _index=0; return data; }
    T*end()     { return data+size; }
    size_t incIndex() { return _index++; } // for (auto d:vect) v1[vect.incIndex()]=d;
    
    // append
    Vect& operator << (const T c) {
        return append(c);
    }
    Vect& operator << (const std::string sn) {
        try {
            switch (*typeid(T).name()) {
                default:
                case 'd': *this << std::stod( sn ); break;
                case 'f': *this << std::stof( sn ); break;
                case 'i': *this << std::stoi( sn ); break;
                case 'm': *this << std::stoul( sn ); break;
            }
        } catch (...) {}
        return *this;
    }
    void operator >> (T &c) {
        assert(_index<size);
        c=data[_index++];
    }
    Vect& operator << (const Vect &other) {
        return append(other);
    }
    Vect& append(const T c) {
        if ( size >= realSize ) {
            realSize+=nchunk;
            data = (T*)realloc(data, realSize*sizeof(T));
        }
        data[size++]=c;
        return *this;
    }
    Vect& append(const Vect &other) {
        auto os=other.size;
        if ( size+os >= realSize ) {
            realSize=size+os;
            data = (T*)realloc(data, realSize*sizeof(T));
        }
        memcpy(data+size, other.data, os*sizeof(T));
        size+=os;
        return *this;
    }
    Vect& clear() {
        dealloc();
        size=realSize=0;
        return *this;
    }
    Vect& zero() {
        memset(data, 0, sizeof(T) * size);
        return *this;
    }
    Vect& one() {
        fill((T)1);
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
    Vect& fit() { // fit current memory to size
        if (size) resize(size);
        return *this;
    }
    
    size_t locate(T c) { // use locate on big enough vectors
        return (size < szSingle) ? stlocate(c) : ThreadedCalc(*this).locate(c);
    }
    size_t bsearch(T c) { //  binary search -> this must be sorted, return index if found, -1 if not
        auto ret = std::bsearch(&c, begin(), size, sizeof(T),
                                [](const void*a, const void*b) ->int { return (int)(*((T*)a) - *((T*)b));});
        return (ret==nullptr) ? -1 : (T*)ret-begin();
    }
    
    
    Vect &operator=(const Vect &other) { // asignment
        if (&other == this) return *this; // check for self-assignment to avoid delloc problems
        
        if (size != other.size) resize(other.size);
        memcpy(data, other.data, size * sizeof(T));
        return *this;
    }
    Vect &operator=(const T c) { // asignment
        std::fill(begin(), end(), c); //  for (auto &d:*this) d=c;
        return *this;
    }
    
    
    // index mutator
    inline T& operator[](size_t index) {
        assert(index>=0 && index<size && "Vect: index error");
        return data[index];
    }
    // index accessor
    inline const T& operator[](size_t index) const {
        assert(index>=0 && index<size  && "Vect: index error");
        return data[index];
    };
    inline const T& operator[](const Vect&v) const { // vect[v] is vect[v._index]
        assert(v._index>=0 && v._index<size && "Vect: index error");
        return data[v._index];
    };
    Vect operator[](VectIndex ixs) const { // index vect by vector of indexes
        Vect v;
        for (auto ix:ixs) v << data[ix];
        return v;
    }
    Vect operator[](std::function<bool(T)> const& lambda) { // index vect by vector of indexes
        Vect v;
        if(size<szSingle) { for (auto d : *this) if (lambda(d)) v << d; }
        else ThreadedCalc(*this).filter(v, lambda);
        return v;
    }
    Vect operator()(size_t from, size_t to) { // sub vector
        return sub(from, to);
    }
    
    // boolean,
    bool operator==(const Vect &other) { // vect==vect
        return (size<szSingle) ? stEQ(other) : ThreadedCalc(*this).eq(other);
    }
    bool operator==(const T c) { // vect == scalar, true if all ==
        bool eq=true;
        if (size<szSingle) { for (auto d:*this) if(d!=c) {eq=false; break;} }
        else               eq=ThreadedCalc(*this).eq(c);
        return eq;
    }
    bool operator!=(const Vect &other) {
        return ! (*this==other);
    }
    bool operator!=(const T c) { // vect != scalar, any different
        bool res=false;
        for (auto d:*this) if(d!=c) {res=true; break;} // find 1st different
        return res;
    }
    bool operator > (const T c) { // vect > scalar, all > scalar -> any <= c
        bool res=true;
        for (auto d:*this) if(d<=c) {res=false; break;}
        return res;
    }
    bool operator < (const T c) { // vect < scalar, all < scalar -> any >= c
        bool res=true;
        for (auto d:*this) if(d>=c) {res=false; break;}
        return res;
    }
    bool operator >= (const T c) { // vect > scalar, all > scalar -> any <= c
        bool res=true;
        for (auto d:*this) if(d<c) {res=false; break;}
        return res;
    }
    bool operator <= (const T c) { // vect < scalar, all < scalar -> any >= c
        bool res=true;
        for (auto d:*this) if(d>c) {res=false; break;}
        return res;
    }
    // compares smallest vectors
    bool operator>(const Vect &other) {
        assert(size==other.size);
        bool gt=true;
        for (auto d:*this) if(d <= other[_index++]) { gt=false; break; }
        return gt;
    }
    bool operator<(const Vect &other) {
        assert(size==other.size);
        bool lt=true;
        for (auto d:*this) if(d >= other[_index++]) { lt=false; break; }
        return lt;
    }
    bool operator>=(const Vect &other) {
        assert(size==other.size);
        bool gte=true;
        for (auto d:*this) if(d < other[_index++]) { gte=false; break; }
        return gte;
    }
    bool operator<=(const Vect &other) {
        assert(size==other.size);
        bool lte=true;
        for (auto d:*this) if(d > other[_index++]) { lte=false; break; }
        return lte;
    }
    
    
    // Vect op constant
    Vect operator+(const T &c) {
        Vect v(*this);
        return (c!=0) ? ((size<szSingle) ? v.staddc(c) : ThreadedCalc(v).evaluate(c, opADD)) : v;
    }
    Vect operator-(const T &c) {
        Vect v(*this);
        return (c!=0) ? ((size<szSingle) ? v.stsubc(c) : ThreadedCalc(v).evaluate(c, opSUB)) : v;
    }
    Vect operator*(const T &c) {
        Vect v(*this);
        return (c!=0) ? ((size<szSingle) ? v.stmulc(c) : ThreadedCalc(v).evaluate(c, opMUL)) : v.zero();
    }
    Vect operator/(const T &c) {
        Vect v(*this);
        return (c!=0) ? ((size<szSingle) ? v.stdivc(c) : ThreadedCalc(v).evaluate(c, opDIV)) : v.zero(); //  assume v/0=0
    }
    Vect operator^(const T &c) {
        Vect v(*this);
        return (c!=0) ? ((size<szSingle) ? v.stpowc(c) : ThreadedCalc(v).evaluate(c, opPOW)) : v.one();
    }
    
    // Vect (+-*/)= const<T>
    
    Vect &operator+=(const T &c) {
        if (c!=0) { if (size<szSingle)  staddc(c); else ThreadedCalc(*this).evaluate(c, opADD); }
        return *this;
    }
    Vect &operator-=(const T &c) {
        if (c!=0) { if (size<szSingle)  stsubc(c); else ThreadedCalc(*this).evaluate(c, opSUB); }
        return *this;
    }
    Vect &operator*=(const T &c) {
        if (c!=0) { if (size<szSingle)  stmulc(c); else ThreadedCalc(*this).evaluate(c, opMUL); }
        else zero();
        return *this;
    }
    Vect &operator/=(const T &c) {
        if (c!=0) { if (size<szSingle)  stdivc(c); else ThreadedCalc(*this).evaluate(c, opDIV); }
        else zero(); // assume  v/=0 -> 0
        return *this;
    }
    
    Vect operator|(const Vect &other) { // join *this and other
        Vect v(size + other.size); // Vect v = *this << other
        memcpy(v.data,      data,       size*sizeof(T)); // v=*this
        memcpy(v.data+size, other.data, other.size*sizeof(T)); // append other
        return v;
    }
    Vect operator&(Vect &other) { // intersection
        Vect v(size + other.size);
        this->sort(); other.sort();
        auto n = std::set_intersection(begin(), end(), other.begin(), other.end(), v.begin()) - v.begin();
        v.resize(n);
        return v;
    }
    
    Vect operator||(const Vect &other) { //
        return *this;
    }
    Vect operator&&(const Vect &other) {
        return *this;
    }
    
    
    // Vect op(+-*/) Vect
    
    Vect operator+(const Vect &other) {
        assert(size==other.size);
        Vect v(*this);
        return (size<szSingle) ? v.stadd(other) : ThreadedCalc(v).evaluate(other, opADD);
    }
    Vect operator-(const Vect &other) {
        assert(size==other.size);
        Vect v(*this);
        return (size<szSingle) ? v.stsub(other) : ThreadedCalc(v).evaluate(other, opSUB);
    }
    Vect operator*(const Vect &other) {
        assert(size==other.size);
        Vect v(*this);
        return (size<szSingle) ? v.stmul(other) : ThreadedCalc(v).evaluate(other, opMUL);
    }
    Vect operator/(const Vect &other) {
        assert(size==other.size);
        Vect v(*this);
        return (size<szSingle) ? v.stdiv(other) : ThreadedCalc(v).evaluate(other, opDIV);
    }
    Vect operator^(const Vect &other) {
        assert(size==other.size);
        Vect v(*this);
        return (size<szSingle) ? v.stpow(other) : ThreadedCalc(v).evaluate(other, opPOW);
    }
    
    // Vect (+-*/)= Vect
    
    Vect &operator+=(const Vect &other) {
        assert(size==other.size);
        if (size<szSingle) stadd(other); else ThreadedCalc(*this).evaluateMutable(other, opADD);
        return *this;
    }
    Vect &operator-=(const Vect &other) {
        assert(size==other.size);
        if (size<szSingle) stsub(other); else ThreadedCalc(*this).evaluateMutable(other, opSUB);
        return *this;
    }
    Vect &operator*=(const Vect &other) {
        assert(size==other.size);
        if (size<szSingle) stmul(other); else ThreadedCalc(*this).evaluateMutable(other, opMUL);
        return *this;
    }
    Vect &operator/=(const Vect &other) {
        assert(size==other.size);
        if (size<szSingle) stdiv(other); else ThreadedCalc(*this).evaluateMutable(other, opDIV);
        return *this;
    }
    
    // ++, -- -> works on _index
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
    
    // single threaded methods
    
    T seqSum_dont_use() { // don't use as it generates precission errors -> apply kahan algo.
        T s=0;
        for (auto d:*this) s+=d;
        return s; // instead use  return kahanSum(data, 0, size);
    }
    std::pair<T,T>stminmax() {
        T mx=*data, mn=*data;
        for (size_t i=0; i<size; i++) {
            if(data[i]>mx) mx=data[i];
            if(data[i]<mn) mn=data[i];
        }
        return std::pair<T,T>(mn,mx);
    }
    std::tuple<T,T,T>seqMinmaxdiff() {
        auto mm=stminmax();
        return std::tuple<T,T,T>(mm.first,mm.second,mm.second-mm.first);
    }
    T stmin() {
        T mn=*data;
        for (size_t i=0; i<size; i++)
            if(data[i]<mn) mn=data[i];
        return mn;
    }
    T stmax() {
        T mx=*data;
        for (size_t i=0; i<size; i++)
            if(data[i]>mx) mx=data[i];
        return mx;
    }
    size_t stlocate(T c) {
        for (size_t i=0; i<size; i++)
            if (c==data[i]) return i;
        return -1;
    }
    bool stEQ(const Vect&other) {
        bool eq = (size==other.size);
        size_t itemNE=0;
        if(eq)
            for (size_t i = 0; i < size; i++)
                if (data[i] != other.data[i]) { eq=false; itemNE=i; break; }
        return eq;
    }
    static Vect stSeq(T from, T to, T inc) {
        assert((from<to && inc>0) && "seq paramaters error");
        Vect v(sizeOfRange(from, to, inc));
        for (size_t i=0; i<v.length(); i++) v.data[i]=(T)i * inc;
        return v;
    }
    Vect& stscale(T from, T to, T min, T diff) {
        auto pdf=to-from;
        for (auto &d:*this) d=from + ((d-min)/diff) * pdf;
        return *this;
    }
    Vect& stnorm(T min, T diff) {
        for (auto &d:*this) d=(d-min)/diff;
        return *this;
    }
    T stsum() { // Single Thread Kahan summation algorithm, sequential sum generates precission errors
        return kahanSum(data, 0, size);
    }
    Vect stsequence(T from, T to) {
        assert(from<to && "sequence error: from > to");
        if(size) {
            T inc=(T)(to-from)/(T)size;
            for (auto &d:*this) d = inc * (_index++);
        }
        return *this;
    }
    
    Vect& staddc(T c) {   for (auto &d:*this) d+=c;  return*this;}
    Vect& stsubc(T c) {   for (auto &d:*this) d-=c;  return*this;}
    Vect& stmulc(T c) {   for (auto &d:*this) d*=c;  return*this;}
    Vect& stdivc(T c) {   for (auto &d:*this) d/=c;  return*this;}
    Vect& stpowc(T c) {   for (auto &d:*this) d=pow(d,c);  return*this;}
    
    Vect& stadd(const Vect&v) { for (auto &d:*this) d+=v[_index++]; return *this;}
    Vect& stsub(const Vect&v) { for (auto &d:*this) d-=v[_index++]; return *this;}
    Vect& stmul(const Vect&v) { for (auto &d:*this) d*=v[_index++]; return *this;}
    Vect& stdiv(const Vect&v) { for (auto &d:*this) d/=v[_index++]; return *this;}
    Vect& stpow(const Vect&v) { for (auto &d:*this) d=pow(d,v[_index++]); return *this;}
    
    void strnd() {
        for (auto &d:*this) d=(T)rand()/(T)RAND_MAX;
    }
    Vect stfunc(Lambda lambda) { // non mutable -> return a copy applying func
        Vect v(*this);
        for (auto &d:v) d=lambda(d);
        return v;
    }
    Vect& stapply(Lambda lambda) { // apply labda func -> mutable, usage: v.apply(sin)
        for (auto &d:*this) d=lambda(d);
        return *this;
    }
    Vect& stapply(std::function<T()> const& lambda) { // usage: v.apply(funcNoArgs)
        for (auto &d:*this) d=lambda();
        return *this;
    }
    Vect& stshuffle() {
        if (size > 1) {
            for (int i=0; i<size; i++)
                std::swap<T>(data[i], data[(rand() % (i + 1))]);
        }
        return *this;
    }
    Vect stfilter(std::function<bool(T)> const& lambda) { // filter by bool lambda conditional expr.
        Vect v;
        for (auto const d:*this) if (lambda(d)) v<<d;
        return v;
    }
    VectIndex stfilterIndex(std::function<bool(T)> const& lambda) { // indexes of selected items
        VectIndex res;
        for (size_t ix=0; ix<size; ix++)
            if (lambda(data[ix])) res<<ix;
        return res;
    }
};


void testVect();

#endif /* Vect_hpp */
