#include <stdio.h>
#include <iostream>
#include <utility>
#include <assert.h>

#include "monom.h"


Mon::Mon(Int c,int g,Mon *s): coef(c), deg(g), sig(s){
    if(deg < 0){
	coef = Int(0, c.getMod());
	deg = 0;
    }

    if(coef.getValue() == 0){
	deg = 0;
    }
}

Mon::Mon(Int c, int g): coef(c), deg(g), sig(nullptr){
    if(deg < 0){
	coef = Int(0, c.getMod());
	deg = 0;
    }
    
    if(coef.getValue() == 0){
	deg = 0;
    }
}

Mon::Mon(): coef(Int(0,0)), deg(0), sig(nullptr){}

// only care about the pointer, the others are removed by the default delete operator
Mon::~Mon(){
    coef.~Int();
    sig.reset();
}

// copy constructor
Mon::Mon(const Mon & m){
    deg = m.deg;
    coef = m.coef;
    sig = m.sig;
}

// default assignment operator
Mon & Mon::operator=(const Mon & m)= default;

//bool Mon::operator==(const Mon & m) const
//bool Mon::operator!=(const Mon & m) const

Mon Mon::operator/(const Int &m) const{
    assert (m.getValue() != 0);
    Mon result = Mon(coef / m, deg);
    return result;
}

Mon Mon::operator /(const Mon & m) const{
    assert (m.getCoef().getValue() != 0);
    Mon result = Mon(coef / m.coef, deg - m.deg);
    return result;
}

Mon Mon::operator * (const Mon & m) const{
    Mon result = Mon(coef * m.coef, deg + m.deg);
    return result;
}

Mon Mon::operator + (const Mon & m) const{
    Mon result = Mon(coef + m.coef, m.deg);
    return result;
}

Mon Mon::operator - (const Mon & m) const{
    Mon result = Mon(coef - m.coef, m.deg);
    return result;
}   

// Procedural functions
void sum(const Mon & m, const Mon & n, Mon & result){
    assert(m.getDeg() == n.getDeg());
    sum(m.getCoef(), n.getCoef(), result.coef);
    result.set_deg(m.getDeg());
}

void sub(const Mon & m, const Mon & n, Mon & result){
    assert(m.getDeg() == n.getDeg());
    sub(m.getCoef(), n.getCoef(), result.coef);
    result.set_deg(m.getDeg());
}

void mult(const Mon & m, const Mon & n, Mon & result){
    result.set_deg(m.getDeg() + n.getDeg());
    mult(m.getCoef(), n.getCoef(), result.coef);
}

void mult(const Mon & m, const Int & b, Mon & result){
    result.set_deg(m.getDeg());
    mult(m.getCoef(), b, result.coef);
}

void mult(const Int & a, const Mon & n, Mon & result){
    result.set_deg(n.getDeg());
    mult(n.getCoef(), a, result.coef);
}

void div(const Mon & m, const Mon & n, Mon & result){
    assert(m.getDeg() >= n.getDeg());
    result.set_deg(m.getDeg() - n.getDeg());
    div(m.getCoef(), n.getCoef(), result.coef);
}

void div(const Mon & m, const Int & b, Mon & result){
    result.set_deg(m.getDeg());
    div(m.getCoef(), b, result.coef);
}

// for printing at standard output
std::ostream& operator << (std::ostream& os, const Mon & mon){
    os << (mon.getCoef()) <<"x^"<< (mon.getDeg());
    return os;
}
