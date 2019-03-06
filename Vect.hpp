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

template <class T>
class Vect {
    typedef Vect<size_t> VectIndex;
    
private:
    size_t realSize=0;
    const int nchunk=128; //  memory increments on each resize
    
    T *data=nullptr;
    size_t size=0;
    size_t _index=0; // used in iterator
    
public:
    // constructors
    Vect() { }
    
    Vect(size_t size) {
        alloc(size);
    }
    Vect(const Vect &other) { // Vect v0(v1), v2=v1;
        alloc(other.size);
        std::copy(&other.data[0], &other.data[0] + size, &data[0]);
    }
    Vect(size_t size, T v[]){
        alloc(size);
        for (auto &d:*this) d=v[_index++];
    }
    
    ~Vect() { dealloc(); }
    
    static Vect rnd(int size) { // norm random
        Vect v(size);
        return v.apply([]() -> T{ return (T)rand()/(T)RAND_MAX; });
    }
    static Vect seq(size_t size) { // 0,1...n-1
        Vect v;
        for (size_t i=0; i<size; i++) v<<i;
        return v;
    }
    static Vect seq(size_t size, T inc) { // 0,inc, 2*inc, ...(n-1)*inc
        Vect v;
        for (size_t i=0; i<size; i++) v<<(T)i*inc;
        return v;
    }
    static Vect seq(double from, double to, double inc) {
        assert((from<to && inc!=0) && "seq paramaters error");
        Vect v;
        for (double i=from; i<to; i+=inc) v<<i;
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
        apply([]() -> T {return (T)rand();});
    }
    void rand01() {
        apply([]() -> T { return (T)rand()/(T)RAND_MAX; });
    }
    T sum() {
        T s=0;
        for (auto d:*this) s+=d;
        return s;
    }
    std::pair<T,T>minmax() {
        T mx=*data, mn=*data;
        for (size_t i=0; i<size; i++) {
            if(data[i]>mx) mx=data[i];
            if(data[i]<mn) mn=data[i];
        }
        return std::pair<T,T>(mn,mx);
    }
    std::tuple<T,T,T>minmaxdiff() {
        auto mm=minmax();
        return std::tuple<T,T,T>(mm.first,mm.second,mm.second-mm.first);
    }
    T min() {
        T mn=*data;
        for (size_t i=0; i<size; i++)
            if(data[i]<mn) mn=data[i];
        return mn;
    }
    T max() {
        T mx=*data;
        for (size_t i=0; i<size; i++)
            if(data[i]>mx) mx=data[i];
        return mx;
    }
    T avg() {
        return sum()/(T)size;
    }
    T diff() {
        auto mm=minmax();
        return mm.first - mm.second;
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
        assert(!(from>to) && "sequence error: from > to");
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
    
    // apply labda func, i.e: v2.apply(sinf).sort().apply([](float x){ return x*x; });
    Vect apply(std::function<T(T)> const& lambda) { // usage: v.apply(sin)
        for (auto &d:*this) d=lambda(d);
        return *this;
    }
    Vect apply(std::function<T()> const& lambda) { // usage: v.apply(funcNoArgs)
        for (auto &d:*this) d=lambda();
        return *this;
    }
    Vect filter(std::function<bool(T)> const& lambda) { // filter by bool lambda conditional expr.
        Vect res;
        for (auto const d:*this) if (lambda(d)) res<<d;
        return res;
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
        if(eq)
            for (size_t i = 0; i < size; i++)
                if (data[i] != other.data[i]) { eq=false; break; }
        return eq;
    }
    bool operator==(const T c) { // vect == scalar, true if all ==
        bool eq=true;
        for (auto d:*this) if(d!=c) {eq=false; break;}
        return eq;
    }
    bool operator!=(const Vect &other) {
        bool neq = (size!=other.size);
        if(neq) {}
        else {
            for (size_t i = 0; i < size; i++)
                if (data[i] == other.data[i]) { neq=false; break; }
        }
        return neq;
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
        Vect v(size);
        for (size_t i = 0; i < size; i++)
            v.data[i] = data[i] + c;
        return v;
    }
    Vect operator-(const T &c) {
        Vect v(size);
        for (size_t i = 0; i < size; i++)
            v.data[i] = data[i] - c;
        return v;
    }
    Vect operator*(const T &c) {
        Vect v(size);
        for (size_t i = 0; i < size; i++)
            v.data[i] = data[i] * c;
        return v;
    }
    Vect operator/(const T &c) {
        Vect v(size);
        for (size_t i = 0; i < size; i++)
            v.data[i] = data[i] / c;
        return v;
    }
    
    // Vect (+-*/)= const<T>
    
    Vect &operator+=(const T &c) {
        for (auto &d:*this) d += c;
        return *this;
    }
    Vect &operator-=(const T &c) {
        for (auto &d:*this) d -= c;
        return *this;
    }
    Vect &operator*=(const T &c) {
        for (auto &d:*this) d *= c;
        return *this;
    }
    Vect &operator/=(const T &c) {
        if(c!=0)
            for (auto &d:*this) d /= c;
        return *this;
    }
    
    
    // Vect op(+-*/) Vect
    
    Vect operator+(const Vect &other) {
        Vect v(size);
        for (size_t i = 0; i < std::min(size, other.size); i++)
            v.data[i] = data[i] + other.data[i];
        return v;
    }
    Vect operator-(const Vect &other) {
        Vect v(size);
        for (size_t i = 0; i < std::min(size, other.size); i++)
            v.data[i] = data[i] - other.data[i];
        return v;
    }
    Vect operator*(const Vect &other) {
        Vect v(size);
        for (size_t i = 0; i < std::min(size, other.size); i++)
            v.data[i] = data[i] * other.data[i];
        return v;
    }
    Vect operator/(const Vect &other) {
        Vect v(size);
        for (size_t i = 0; i < std::min(size, other.size); i++)
            v.data[i] = data[i] / other.data[i];
        return v;
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
