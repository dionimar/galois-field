#ifndef MOD_ARITH_H
#define MOD_ARITH_H


#include <stdio.h>
#include <iostream>
#include <exception>
#include <cstdlib>
#include <vector>

#include "polynom.h"
//using namespace std;

std::vector<int> prime_divisors(int n);

template<class Type>
Type mcd(Type a, Type b);


template<class Type>
void extEuclides(Type a, Type b, Type & auxA, Type & auxB, Type & mcd);

void eea(Poly f, Poly g, Poly s, Poly t);

/* template<class Type> */
/* Type invMult(Type a, Type b, Type & auxA, Type & auxB); */

template<class Type>
Type garberRemainder(Type u[], Type m[], int l);

int garnerRemainder(int u[], int m[], int l);

Poly& mod_power(Poly & p, int n);

Poly& gcd_poly(Poly & f, Poly & g);

/* Poly mod_composition(const Poly & g, const Poly & h, const Poly & mod); */

#endif
