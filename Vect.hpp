//
//  Vect.hpp
//  general purpose numerical vect class
//  debug array: p *(float(*)[10])data
// usage: see testVect() func

#ifndef Vect_hpp
#define Vect_hpp

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <string>
#include <functional>

using std::copy;

template <class T>
class Vect {
    
private:
    size_t realSize=0;
    const int nchunk=128; //  memory increments on each resize
    
    T *data=nullptr;
    size_t size=0;
    size_t _index=0;
    
public:
    // constructors
    Vect() { }
    
    Vect(size_t size) {
        alloc(size);
    }
    Vect(const Vect &other) { // Vect v0(v1), v2=v1;
        alloc(other.size);
        copy(&other.data[0], &other.data[0] + size, &data[0]);
    }
    Vect(size_t size, T v[]){
        alloc(size);
        for (auto &d:*this) d=v[_index++];
    }
    
    ~Vect() { dealloc(); }
    
    static Vect rnd(int size) { // norm random
        Vect v(size);
        return v.apply([](float x){ return (T)rand()/(T)RAND_MAX; });
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
        if ( (data = (T*)calloc(size, sizeof(T))) == nullptr ) size=0;
        realSize = this->size = size;
    }
    void resize(size_t size) {
        if ( (data = (T*)realloc(data, size*sizeof(T))) == nullptr ) size=0;
        realSize = this->size = size;
    }
    
public:
    size_t length() {
        return size;
    }
    
    void random() {
        apply([](T x){return (T)rand();});
    }
    void rand01() {
        apply([](float x){ return (T)rand()/(T)RAND_MAX; });
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
        
        for (auto &d:*this)
            d = from + ((d-mn) / df) * pdf;
        return *this;
    }
    Vect fill(T value) {
        for (auto &d:*this) d = value;
        return *this;
    }
    Vect sequence(T from, T to) {
        T inc=(T)(to-from)/(T)size;
        for (auto &d:*this) d = inc * (_index++);
        return *this;
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
    Vect apply(std::function<T(T)> const& lambda) { // usage: v.apply(sinf)
        for (auto &d:*this) d=lambda(d);
        return *this;
    }
    
    T*begin() { _index=0; return data; }
    T*end() { return data+size; }
    
    // append, void return as Vect return -> very bad performnce
    // so v0 << v1 << v2 is v0<<v1; v0<<v2;
    void operator << (const T c) {
        append(c);
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
    
    Vect &operator=(const Vect &other) { // asignment
        // check for self-assignment
        if (&other == this) return *this;
        
        if (size != other.size) resize(other.size);
        std::copy(&other.data[0], &other.data[0] + size, &data[0]);
        return *this;
    }
    
    // index mutator
    inline T& operator[](size_t index) { return data[index]; }
    // index accessor
    inline const T& operator[](size_t index) const { return data[index]; };
    
    // boolean,
    bool operator==(const Vect &other) {
        bool eq = (size==other.size);
        if(eq)
            for (size_t i = 0; i < size; i++)
                if (data[i] != other.data[i]) { eq=false; break; }
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
    
    // ++, -- Vect, prefix only,
    // postfix: Vect &operator++(int) not implemented->bad performance
    Vect &operator++() {
        for (auto &d:*this) d++;
        return *this;
    }
    Vect &operator--() {
        for (auto &d:*this) d--;
        return *this;
    }
};


void testVect();

#endif /* Vect_hpp */
