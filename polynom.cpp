#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <memory>
#include <vector>
#include <random>
#include <utility>
#include <assert.h>

#include "polynom.h"

// basic constructors
Poly::Poly(){
    std::shared_ptr<Mon> first = std::make_shared<Mon>();
    std::shared_ptr<Mon> last = std::make_shared<Mon>();
}

Poly::Poly(Int coef,int deg){
    *first = Mon(coef,deg,nullptr);
    last = first;
}

Poly::Poly(int coef, int deg){
    auto new_coef = new Int(coef, 0);
    *first = Mon(*new_coef,deg,nullptr);
    last=first;
}

// marked explicit for avoid implicit conversions
Poly::Poly(Int n){
    *first = Mon(n,0,nullptr);
    last=first;
}

// reset shared pointers
Poly::~Poly(){
    last.reset();
    first.reset();
}

//Copy constructor
Poly::Poly(const Poly &q){
    std::shared_ptr<Mon> act = q.first;
    while(act != nullptr){
	Mon aux = *act;
        ponMonomio(aux);
        act=act->sig;
    }
}

//Asignment operator
Poly Poly::operator = (const Poly & q){
    last.reset();
    first.reset();

    last = std::move(q.last);
    first = std::move(q.first);
    return *this;
}

// actually deprecated
bool Poly::isNull(){
    return (first == nullptr);
}

int Poly::getDeg()const{
    if(first != nullptr){
        return first->getDeg();
    }
    else return 0;
}

std::vector<Int> poliToVector_n(const Poly & p, const int n){
    // better performance using vector.reserve(n)   ?????
    std::shared_ptr<Mon> aux = p.first;
    int prim = aux->getCoef().getMod();
    std::vector<Int> result = std::vector<Int>(n, Int(0, prim));
    while(aux != nullptr){
	result[aux->getDeg()] = aux->getCoef();
	aux = aux->sig;
    }
    return result;
}

// insert monom on polynom structure, taking care of the position to get de structure ordered by degree
void Poly::ponMonomio(const Mon &m){
    if(first != nullptr){
	if(first->getCoef().getValue() == 0){
	    first.reset();
	    first = std::make_shared<Mon>(m.getCoef(), m.getDeg(), nullptr);
	    last = first;
	}
	
	else{
	    if(m.getCoef().getValue() != 0){
		last->sig = std::make_shared<Mon>(m.getCoef(), m.getDeg(), nullptr);
		last = last->sig;
	    }
	}
    }
    
    else{
	first = std::make_shared<Mon>(m.getCoef(), m.getDeg(), nullptr);
	last = first;
    }
}

// return the leading unit
// Int Poly::lc() const{
//     return Int(first->getCoef());
// }

// return the monic polynom dividing by its leading coeficient
Poly Poly::normal() const{
    if(first->getCoef().getValue() != 0){
	Mon lc_term = Mon(Int(this->first->getCoef().getValue(), this->first->getCoef().getMod()),0);
	return ((*this)/lc_term);
    }
    else{
	return *this;
    }
}

// Used as basic constructor, keeping the structure ordered by degree using ponMonomio
Poly Poly::operator+(const Poly & p) const{
    Poly result;
    std::shared_ptr<Mon> i = first;
    std::shared_ptr<Mon> j = p.first;
    while (i != nullptr  && j != nullptr){
	
        if ((i->getDeg()) > (j->getDeg())){
	    result.ponMonomio(*i);
	    i=i->sig;
        }
	
        else if ((i->getDeg()) < (j->getDeg())){
	    result.ponMonomio(*j);
            j=j->sig;
        }
	
        else{
	    Mon aux = Mon(*i + *j);
	    result.ponMonomio(aux);
            i=i->sig;
            j=j->sig;
        }
    }

    while(i != nullptr){
	Mon aux = *i;
        result.ponMonomio(aux);
        i=i->sig;
    }
    
    while(j != nullptr){
	Mon aux = *j;
        result.ponMonomio(aux);
        j=j->sig;
    }
    return result;
}

Poly Poly::operator-(const Poly & p) const{
    Poly P = p;

    if(this->first == nullptr){
    	Poly cero = Poly(Int(0,P.first->getCoef().getMod()),0);
    	return (cero - P);
    }
    
    else if(P.first == nullptr){
    	return *this;
    }

    else{
	Poly result = Poly(Int(0,p.first->getCoef().getMod()));
	std::shared_ptr<Mon> i = first;
	std::shared_ptr<Mon> j = p.first;
	while (i != nullptr  && j != nullptr){
	
	    if ((i->getDeg()) > (j->getDeg())){
		result.ponMonomio(*i);
		i=i->sig;
	    }
	
	    else if ((i->getDeg()) < (j->getDeg())){
		Mon aux = Mon(Int(0-(j->getCoef().getValue()),j->getCoef().getMod()), j->getDeg());
		result.ponMonomio(aux);
		j=j->sig;
	    }
	
	    else{
     		Mon aux = Mon(*i - *j);
		result.ponMonomio(aux);
		i=i->sig;
		j=j->sig;
	    }
	}

	while(i != nullptr){
	    Mon aux = Mon(*i);
	    result.ponMonomio(aux);
	    i=i->sig;
	}
    
	while(j != nullptr){
	    Mon aux = Mon(Int(0-(j->getCoef().getValue()),j->getCoef().getMod()), j->getDeg());
	    result.ponMonomio(aux);
	    j=j->sig;
	}
	return result;
    }   
}

Poly Poly::operator*(const Mon & m) const{
    Mon cero = Mon(Int(0,m.getCoef().getMod()),0,nullptr);
    Mon uno = Mon(Int(1,m.getCoef().getMod()),0,nullptr);
    Mon M = m;
    Poly result;

    if(M == cero){
	result = Poly(Int(0,m.getCoef().getMod()),0);
	return result;
    }

    if(M == uno){
	return *this;
    }

    else{
	std::shared_ptr<Mon> act = first;
	while(act!=nullptr){
	    Mon aux = Mon((*act)*M);
	    result.ponMonomio(aux);
	    act=act->sig;
	}
	return result;	
    }
}

Poly Poly::operator*(const Poly & p) const{
    Poly cero = Poly(Int(0,this->first->getCoef().getMod()),0);
    Poly uno = Poly(Int(1,this->first->getCoef().getMod()),0);
    Poly result;
    Poly P = p;

    if(P == cero){
	result = cero;
	return result;
    }

    if(P == uno){
	result = *this;
	return result;
    }

    else{
	Poly auxThis = Poly(*this);
	std::shared_ptr<Mon> act = p.first;	
	while(act != nullptr){  
	    Poly aux = Poly(auxThis * (*act));
	    result = result + aux;
	    act = act->sig;
	}
	return result;
    }    
}

// void mult(const Poly & f, const Poly & g, Poly & result){
//     int sum_deg = f.getDeg() + g.getDeg();
//     if(sum_deg <= 2){
// 	Poly F = f;
// 	Poly G = g;
// 	result = F * G;
//     }

//     else{
	
//     }
// }

Poly Poly::operator/(const Int &m) const{
    Poly result;
    std::shared_ptr<Mon> act = first;
    while(act != nullptr){
	Mon aux = Mon((*act)/(m));
        result.ponMonomio(aux);
        act = act->sig;
    }
    return result;
}

Poly Poly::operator /(const Mon & m) const{
    assert (m.getCoef().getValue() != 0);
    Poly result;
    std::shared_ptr<Mon> act = first;
    while(act != nullptr){
	Mon aux = Mon((*act)/(m));
        result.ponMonomio(aux);
        act = act->sig;
    }
    return result;
}

Poly Poly::operator /(const Poly & q) const{
    assert (q.first->getCoef().getValue() != 0);
    Poly Q = q;
    int prim = q.first->getCoef().getMod();
    int n = this->getDeg();
    int m = q.getDeg();
    Poly remainder, quo;
    remainder = *this;
    
    Int u = Int(invMult(q.first->getCoef().getValue(), prim), prim);

    for(int i = n-m; i >= 0; i--){
	if(remainder.getDeg() == m + i){
	    Int aux = remainder.lc() * u;
	    quo = quo + Poly(aux, i);
	    remainder = remainder  - (Poly(aux,i)*Q);
	}
    }
    return quo;
}

Poly Poly::operator %(const Poly & q) const{
    Poly Q = q;
    int prim = q.first->getCoef().getMod();
    int n = this->getDeg();
    int m = q.getDeg();
    Poly remainder, quo;
    remainder = *this;
    
    Int u = Int(invMult(q.first->getCoef().getValue(), prim), prim);

    for(int i = n-m; i >= 0; i--){
	if(remainder.getDeg() == m + i){
	    Int aux = remainder.lc() * u;
	    quo = quo + Poly(aux, i);
	    remainder = remainder  - Poly(aux,i)*Q;
	}
    }
    return remainder;    
}

bool Poly::operator==(const Poly & q) const{   
    bool T = true;
    std::shared_ptr<Mon> me = this->first;
    std::shared_ptr<Mon> p = q.first;
    if(me == nullptr){
    	if(p != nullptr){
    	    T = false;
    	}
    }
    
    else{
	if(p == nullptr){
	    T = false;
	}

	else{
	    while((me != nullptr && p != nullptr) && T){
		T = (*me == *p);
		me = me->sig;
		p = p->sig;
	    }

	    if(((me == nullptr && p == nullptr)) == false){
		T = false;
	    }
	}    
    }    
    return T;
}

bool Poly::operator!=(const Poly & q) const{
    return not(*this==q);
}


// Prodecural functions
void sum(const Poly & f, const Poly & g, Poly & result){
//    result = Poly();
    std::shared_ptr<Mon> i = f.first;
    std::shared_ptr<Mon> j = g.first;
    
    while(i != nullptr  && j != nullptr){
	
        if((i->getDeg()) > (j->getDeg())){
	    result.ponMonomio(*i);
	    i=i->sig;
        }
	
        else if ((i->getDeg()) < (j->getDeg())){
	    result.ponMonomio(*j);
            j=j->sig;
        }
	
        else{
	    Mon aux = Mon(*i + *j);
	    result.ponMonomio(aux);
            i=i->sig;
            j=j->sig;
        }
    }

    while(i != nullptr){
	Mon aux = *i;
        result.ponMonomio(aux);
        i=i->sig;
    }
    
    while(j != nullptr){
	Mon aux = *j;
        result.ponMonomio(aux);
        j=j->sig;
    }
}

void sum(const Poly & f, const Mon & m, Poly & result){
    result = Poly();
    std::shared_ptr<Mon> i = f.first;
    if(i->getDeg() < m.getDeg()){
	
	result.ponMonomio(m);
	
	while(i != nullptr){
	    result.ponMonomio(*i);
	    i = i->sig;
	}
    }

    else{
	while(i != nullptr && i->getDeg() > m.getDeg()){
	    result.ponMonomio(*i);
	    i = i->sig;
	}
	
	if(i->getDeg() == m.getDeg()){
	    result.ponMonomio(*i + m);
	    i = i->sig;
	    while(i != nullptr){
		result.ponMonomio(*i);
		i = i->sig;
	    }
	}
	else{
	    result.ponMonomio(m);
	    while(i != nullptr){
		result.ponMonomio(*i);
		i = i->sig;
	    }
	}
    }
}

void sum(const Poly & f, const Int & a, Poly & result){
//    result = Poly();
    Mon aux = Mon(a, 0);
    std::shared_ptr<Mon> j = f.first;
    
    while(j != nullptr){
	result.ponMonomio(*j);
	j = j->sig;
    }
    
    if(result.last->getDeg() != 0){
	result.ponMonomio(aux);
    }

    else{
	sum(*f.last, aux, *result.last);
    }
}

void sub(const Poly & f, const Poly & g, Poly & result){
//    result = Poly();
    std::shared_ptr<Mon> i = f.first;
    std::shared_ptr<Mon> j = g.first;
    
    while(i != nullptr  && j != nullptr){
	
        if((i->getDeg()) > (j->getDeg())){
	    result.ponMonomio(*i);
	    i=i->sig;
        }
	
        else if ((i->getDeg()) < (j->getDeg())){
	    result.ponMonomio(Mon(Int(0, f.first->getCoef().getMod()), j->getDeg()) - *j);
            j=j->sig;
        }
	
        else{
	    result.ponMonomio(*i - *j);
            i=i->sig;
            j=j->sig;
        }
    }

    while(i != nullptr){
        result.ponMonomio(*i);
        i=i->sig;
    }
    
    while(j != nullptr){
        result.ponMonomio(Mon(Int(0, f.first->getCoef().getMod()), j->getDeg()) - *j);
        j=j->sig;
    }    
}

void sub(const Poly & f, const Mon & m, Poly & result){
    result = Poly();
    std::shared_ptr<Mon> i = f.first;
    if(i->getDeg() < m.getDeg()){
	
	result.ponMonomio(Mon(Int(0, m.getCoef().getMod()), m.getDeg()) - m);
	
	while(i != nullptr){
	    result.ponMonomio(*i);
	    i = i->sig;
	}
    }

    else{
	while(i != nullptr && i->getDeg() > m.getDeg()){
	    result.ponMonomio(*i);
	    i = i->sig;
	}
	
	if(i->getDeg() == m.getDeg()){
	    result.ponMonomio(*i - m);
	    i = i->sig;
	    while(i != nullptr){
		result.ponMonomio(*i);
		i = i->sig;
	    }
	}
	else{
	    result.ponMonomio(m);
	    while(i != nullptr){
		result.ponMonomio(*i);
		i = i->sig;
	    }
	}
    }
}

void mult(const Poly & f, const Mon & m, Poly & result){
    result = Poly();
    std::shared_ptr<Mon> i = f.first;
    while(i != nullptr){
	result.ponMonomio(*i * m);
	i = i->sig;
    }
}

// void mult(const Poly & f, const Poly & g, Poly & result){
//     result = Poly();
//     long long signed int mod = (f.first->getCoef().getMod() > 0) ? (f.first->getCoef().getMod()) : (g.first->getCoef().getMod());
    
//     if(f.first->getCoef().getValue() == 0 && g.first->getCoef().getValue() == 0){
// 	result->first = Mon(Int(0, mod), 0);
//     }

//     else if((f.first->getCoef().getValue() == 1 && f.first->getDeg() == 0) || (f.first->getCoef().getValue() == 1 && f.first->getDeg() == 0)){
// 	result = (f.first->getDeg() == 0) ? g : f;
//     }

//     else{
// //	Poly auxThis = Poly(*this);
// 	std::shared_ptr<Mon> i = f.first;
// //	std::shared_ptr<Mon> j = g.first;
// 	while(i != nullptr){  
// 	    mult(*i, g, result);
// 	    i = i->sig;
// 	}
// 	return result;
//     }    
// }

// Repeated squaring algorithm for calculating the power
void mod_power(const Poly & p, const int n, Poly & result){
//    std::cout<<"power n = "<<n<<std::endl;
    result = p;

    if(n == 0){
	result = Poly(Int(1,p.first->getCoef().getMod()),0);
    }

    else if(n == 1){
	result = p;
    }

    else{
	std::vector<int> bin_n_rep = intToBin(n);
	int k = bin_n_rep.size();
//    Poly aux1;
    
	for(int i = 1; i < k; i++){
//	std::cout<<i<<"/"<<k<<std::endl;
	    if(bin_n_rep[i] == 1){
//	    std::cout<<"    mult deg :"<<aux0.getDeg()<<" P deg :"<<P.getDeg()<<std::endl;
		result = result * result * p;
//	    std::cout<<aux1<<std::endl;
	    }
	    else{
//	    std::cout<<"    mult deg :"<<aux0.getDeg()<<std::endl;
		result = result * result;
//	    std::cout<<aux1<<std::endl;
	    }
//	aux0 = result;
	}
    }
//    result = aux1;
//    std::cout<<"return mod power"<<result<<std::endl;
//    return result;
}

Poly mod_power(const Poly & p, const int n){
    Poly result;
    mod_power(p, n, result);
    return result;
}

// Extended Euclides Algorithm
void extEuclidesPoly(const Poly & a, const Poly & b, Poly & auxA, Poly & auxB, Poly & mcd){

    int mod = a.first->getCoef().getMod();
    Int ro2_aux;
    Poly q, rem;
    
    std::vector<Poly> Ro;
    Ro.push_back(Poly(a.lc(),0));
    Ro.push_back(Poly(b.lc(),0));
    std::vector<Poly> R;
    R.push_back(a.normal());
    R.push_back(b.normal());
    std::vector<Poly> S;
    S.push_back(Poly(Int(invMult(a.lc().getValue(), mod), mod),0));
    S.push_back(Poly(Int(0,mod),0));
    std::vector<Poly> T;
    T.push_back(Poly(Int(0,mod),0));
    T.push_back(Poly(Int(invMult(b.lc().getValue(), mod), mod),0));
    
    int i = 1;
    
    while(R[i].getDeg() > 0){
	q = R[i-1] / R[i];
	rem = R[i-1] % R[i];	
	
	if(rem.lc().getValue() != 0){
	    ro2_aux = rem.lc();
	}
	else{
	    ro2_aux = Int(1, mod);
	}

	Ro.push_back(Poly(ro2_aux,0));
	R.push_back(rem.normal());
	S.push_back((S[i-1] - (q*S[i]))/Ro[i+1]);
	T.push_back((T[i-1] - (q*T[i]))/Ro[i+1]);
		
	i++;
    }

    if(R[i].first->getCoef().getValue() == 0){
	mcd = R[i-1];
	auxA = S[i-1];
	auxB = T[i-1];
    }
    else{
	mcd = R[i];
	auxA = S[i];
	auxB = T[i];
    }
}

// Inverse using Extended Euclides Algorithm
Poly invMultPoly(const Poly & a, const Poly & Fmod){
    Poly auxA, auxB, gcd;
    extEuclidesPoly(a, Fmod, auxA, auxB, gcd);
    if(gcd.first->getCoef().getValue() == 1 && gcd.getDeg() == 0){
	return auxA;
    }
    else{
	std::cout<<"Polynomial inverse ERROR"<<std::endl;
    }
}

// GCD(classical version) computation using Extended Euclides Algorithm
void gcd_poly(const Poly & f, const Poly & g, Poly & result){
    int mod = f.first->getCoef().getMod();
    Int ro2_aux;
    Poly q, rem;
    std::vector<Poly> R;
    R.push_back(f.normal());
    R.push_back(g.normal());
    int i = 1;
    
    while(R[i].getDeg() > 0){
	rem = R[i-1] % R[i];
	R.push_back(rem.normal());	
	i++;
    }

    if(R[i].first->getCoef().getValue() == 0){
	result = R[i-1];
    }
    
    else{
	result = R[i];
    }
}

// generates a random polynomial of degree less than deg and coefficients with modulus mod
Poly rand_poly(const int deg, const long long int mod){
    std::uniform_int_distribution<> dis(0, mod-1);
    Poly result = Poly(Int(1,mod),deg);   
    int i = deg-1;
    
    while(i >= 0){
    	if(deg <= 1){
    	    result = result + Poly(Int(random_number(1, mod-1), mod), deg);
    	}

    	else{
    	    result = result + Poly(Int(random_number(0, mod-1), mod), random_number(0, deg-1));
    	}	
    	i--;
    }
    return result;
}



// Differential polynomial
Poly diff(const Poly & f){
    Poly F = f;
    Poly diff;
    std::shared_ptr<Mon> aux = F.first;
    while(aux != nullptr){
	if(aux->getDeg() == 0){
	    aux = aux->sig;
	}
	else{
	    Int add = (aux->getCoef())*(aux->getDeg());
	    if(add.getValue() != 0){
		diff = diff + Poly(add, aux->getDeg() -1);
	    }
	    aux = aux->sig;
	}
    }
    return diff;
}

// auxiliar functions for squarefree decomposition (on Galois Field(finite fields))
std::vector<Poly> fill(std::vector<Poly> vec, const int p){
    // int p = mod in galField, for simplicity we use this fact
    //INPUT (h1,h2,...,hs)
    //OUTPUT (1,...,1,h1,1,...,1,h2,1,...,1,hs)
    //          ^
    //         p-1 tiems
    
    Poly uno = Poly(Int(1,p),0);
    std::vector<Poly> result = std::vector<Poly>(p*vec.size(),uno);
    int i = 0;
    while(i<=vec.size()){
	if(i !=0){
	    result[p*i-1] = vec[i-1];
	}
	i++;
    }
    return result;
}

std::vector<Poly> prod(const std::vector<Poly> u, const std::vector<Poly> v){
    //requires v.size() >= u.size()
    //OUTPUT     (u[1]*v[1],...,u[i]*v[i],v[i+1],...)
    std::vector<Poly> result = std::vector<Poly>(v.size());
    for(int i = 0; i<v.size(); i++){
	if(i<u.size()){
	    result[i] = u[i]*v[i];
	}
	else{
	    result[i] = v[i];
	}
	
    }
    return result;
}

// standard output printing
std::ostream& operator << (std::ostream& os, const Poly& p){
    std::shared_ptr<Mon> act = p.first;
    
    if (act!=nullptr){
        os<<act->getCoef()<<"x"<<"^"<<act->getDeg();
        act=act->sig;
	
        while(act!=nullptr){
	    if(act->getCoef().getValue()>=0){
		os<<" +";
	    }	    
	    os<<" "<<act->getCoef().getValue()<<"x"<<"^"<<act->getDeg();
            act=act->sig;
        }
    }
    
    else{
	os<<0<<"(nullptr)";
    }
    
    return os;
}
