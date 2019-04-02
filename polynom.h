#ifndef POLYNOM_H
#define POLYNOM_H

#include <stdio.h>
#include <iostream>
#include <memory>

#include "monom.h"


// This class works like a linked list:
//
//        Mon---(sig)--->Mon---(sig)--->...---(sig)--->Mon
//         ^                                            ^
//         |                                            |
//         |                                            |
//      first                                       last



class Poly{
public:
    // Basics parameters:
    // first: sahred pointer to first Monom (greater degree)
    // last:  shared pointer to last Monom (less degree)
    std::shared_ptr<Mon> first = std::make_shared<Mon>();
    std::shared_ptr<Mon> last = std::make_shared<Mon>();
    // basic constructors
    // marked explicit for avoid implicit conversions
    Poly(); 
    explicit Poly(Int coef, int deg);
    explicit Poly(int coef, int deg);
    explicit Poly(Int n);
    // reset shared pointers
    ~Poly();
    //Copy constructor
    Poly(const Poly & q);
    //Assignment operator
    Poly operator=(const Poly & q);
    // actually deprecated
    bool isNull();
    int getDeg() const;
    // insert monom on polynom structure, taking care of the position to get de structure ordered by degree
    void ponMonomio(const Mon &m);
    // return leading unit
    inline Int lc() const{
	return Int(first->getCoef());
    }
    // return the monic polynom dividing by its leading coeficient
    Poly normal() const;
    // Used as basic constructor
    Poly operator +(const Poly & p) const;
    Poly operator -(const Poly & p) const;
    Poly operator *(const Mon & m) const;
    Poly operator *(const Poly & p) const;
    Poly operator /(const Int & m) const;
    Poly operator /(const Mon & m) const;
    Poly operator /(const Poly & q) const;
    Poly operator %(const Poly & q) const;
    bool operator ==(const Poly & q) const;
    bool operator !=(const Poly & q) const;
};

// Procedural functions for basic operations
// keeps the order arg1 (op) arg2 (=) result
void sum(const Poly & f, const Poly & g, Poly & result);
void sum(const Poly & f, const Mon & m, Poly & result);
inline void sum(const Mon & m, const Poly & g, Poly & result){
    sum(g, m, result);
}
void sum(const Poly & f, const Int & a, Poly & result);
inline void sum(const Int & a, const Poly & g, Poly & result){
    sum(g, a, result);
}

void sub(const Poly & f, const Poly & g, Poly & result);
void sub(const Poly & f, const Mon & m, Poly & result);
//void sub(const Mon & m, const Poly & f, Poly & result);
inline void sub(const Poly & f, const Int & a, Poly & result){
    sub(f, Mon(a, 0), result);
}
//void sub(const Int & a, const Poly & f, Poly & result);

// in the future the FFT will be implemented to take only nlogn cost
void mult(const Poly & f, const Poly & g, Poly & result);
void mult(const Poly & f, const Mon & m, Poly & result);
inline void mult(const Mon & m, const Poly & f, Poly & result){
    mult(f, m, result);
}
void mult(const Poly & f, const Int & a, Poly & result);
void mult(const Int & a, const Poly & f, Poly & result);

void div(const Poly & f, const Poly & g, Poly & result);
void div(const Poly & f, const Mon & m, Poly & result);
void div(const Poly & f, const Int & a, Poly & result);

void normal(const Poly & f, Poly & result);


//void mult(const Poly & f, const Poly & g, Poly & result);

//int random_number(int low, int high);

void mod_power(const Poly & p, const int n, Poly & result);
Poly mod_power(const Poly & p, const int n);

// Extended Euclides Algorithm
void extEuclidesPoly(const Poly& a, const Poly& b, Poly & auxA, Poly & auxB, Poly & mcd);

// Inverse using Extended Euclides Algorithm
Poly invMultPoly(const Poly & a, const Poly & Fmod);

// GCD(classical version) computation using Extended Euclides Algorithm
void gcd_poly(const Poly & f, const Poly & g, Poly & result);

// Random polynomial
Poly rand_poly(const int deg, const long long int mod);

std::vector<Int> poliToVector_n(const Poly & p, const int n);

// Differential polynomial
Poly diff(const Poly & f);

// auxiliar functions for squarefree decomposition (on Galois Field(finite fields))
std::vector<Poly> fill(std::vector<Poly> vec, const int p);
std::vector<Poly> prod(const std::vector<Poly> u, const std::vector<Poly> v);

// standard output printing
std::ostream& operator << (std::ostream& os, const Poly& p);



#endif
