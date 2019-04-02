#ifndef MONOM_H
#define MONOM_H

#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <memory>

#include "int.h"

// This class assume the operations behind the same modulus

class Mon{
public:
    // Basics parameters:
    // coef: Int coeficient (integer with modulus)
    // deg: integer degree
    // shared pointer to next monom (like linked list) managed by polynom
    Int coef;
    int deg;
    
public:

    // shared pointer to next monom (basics structure for polynomials)
    std::shared_ptr<Mon> sig = std::make_shared<Mon>();

    // Basic constructors
    explicit Mon(Int c,int g,Mon *s);    
    explicit Mon(Int c, int g);
    Mon();
    // only care about the pointer, the others are removed by the default delete operator
    ~Mon();
    // copy constructor
    Mon(const Mon & m);
    // default assignment operator
    Mon & operator=(const Mon & m);
    // inlined for better performance 
    inline int getDeg() const {return deg;};
    inline Int getCoef() const {return coef;};

    // set functions
    inline void set_coef(const Int & a){
	coef.set_mod(a.getMod());
	coef.set_value(a.getValue());
    }

    inline void set_deg(const int n){
	deg = n;
    }
  
    inline bool operator==(const Mon & m) const{
	return (deg == m.deg && coef == m.coef);
    }
  
    inline bool operator!=(const Mon & m) const{
	return (deg != m.deg || coef != m.coef);
    }
    
    Mon operator /(const Int & m) const;
    Mon operator /(const Mon & m) const;
    Mon operator *(const Mon & m) const;
    Mon operator +(const Mon & m) const;
    Mon operator -(const Mon & m) const;
};

// Procedural functions
void sum(const Mon & m, const Mon & n, Mon & result);
void sub(const Mon & m, const Mon & n, Mon & result);

void mult(const Mon & m, const Mon & n, Mon & result);

void mult(const Mon & m, const Int & b, Mon & result);
void mult(const Int & a, const Mon & n, Mon & result);

void div(const Mon & m, const Mon & n, Mon & result);

void div(const Mon & m, const Int & b, Mon & result);


//void power(const Mon & m, const long long signed int n, Mon & result);


// for printing at standard output
std::ostream& operator << (std::ostream& os, const Mon& mon);


#endif
