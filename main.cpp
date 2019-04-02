#include <stdio.h>
#include <iostream>
#include <exception>
#include <cstdlib>
#include <memory>
#include <random>
#include <time.h>

#include "int.h"
#include "monom.h"
#include "polynom.h"
#include "gal_field.h"
#include "matrix.h"


using namespace std;

int main(){
    std::cout<<endl<<endl;

// ---------------------------------------------TESTING---------------------------------------
    long long signed int p = 3;//9422921;
    bool T = true;
    Poly a, b, c, d;
    int MAX_OPS = 100;
    int deg_poly = 10;

    for(int i = 0; i < MAX_OPS; i++){
	clock_t start = clock();
	a = rand_poly(deg_poly, p);

	
	//------------Ext Euclides Testing-------------------
	// std::cout<<"--------------------------"<<std::endl;
	// std::cout<<std::endl;
	// std::cout<<"Testing Extended Euclides Algorithm:"<<std::endl;
	// std::cout<<"GCD("<<a<<", "<<b<<")"<<std::endl;
	// Poly test, auxA, auxB, gcd;
    	// b = rand_poly(deg_poly - 1,p);
   	// extEuclidesPoly(a, b, auxA, auxB, gcd);
   	// test = a*auxA + b*auxB;
	// std::cout<<std::endl;
	// std::cout<<std::endl;
    	// std::cout<<"            a*auxA + b*auxB = "<<test<<std::endl;
    	// std::cout<<"                        gcd = "<<gcd<<std::endl;
	// std::cout<<std::endl;
	// clock_t end = clock();
	// double time = static_cast<double> (end-start) / CLOCKS_PER_SEC * 1000.0;
	// std::cout<<"            time exec.(ms): "<<time<<std::endl;
        //--------------------------------------------------

	//------------Polynomial Inverse Test---------------
	std::cout<<"-----------------------------"<<std::endl;
	std::cout<<std::endl;
	Poly b = rand_poly(deg_poly - 1, p);
	Poly test;
	galField gal = galField(p, 1);
	while(!gal.isIrred(b)){
	    b = rand_poly(deg_poly - 1, p);
	}
	Poly inv = invMultPoly(a, b);
	std::cout<<"Testing Polynomial Inverse of "<<a<<" mod "<<b<<std::endl;
	std::cout<<std::endl;
	std::cout<<"            Polynomial Inverse "<<inv<<std::endl;
	std::cout<<std::endl;
	test = (a * inv) % b;
	std::cout<<"            Mult. TEST: "<<test<<std::endl;
	clock_t end = clock();
	double time = static_cast<double> (end-start) / CLOCKS_PER_SEC * 1000.0;
	std::cout<<"            time exec.(ms): "<<time<<std::endl;
	std::cout<<std::endl;
	//--------------------------------------------------

	//------------Irreducibility Test-------------------
	// Care of big primes and degrees!!
	// std::cout<<"--------------------------"<<std::endl;
	// std::cout<<std::endl;
	// bool is_irred;
	// galField gal = galField(p, 1);
	// is_irred = gal.isIrred(a);
	// std::cout<<"Testing Irreducibility of : "<<a<<std::endl;
	// std::cout<<"            Irred = "<<is_irred<<std::endl;
	// std::cout<<std::endl;
	// clock_t end = clock();
	// double time = static_cast<double> (end-start) / CLOCKS_PER_SEC * 1000.0;
	// std::cout<<"            time exec.(ms): "<<time<<std::endl;
	//--------------------------------------------------

	//------------General Factorization Test------------
	// std::cout<<"--------------------------"<<std::endl;
	// std::cout<<std::endl;
	// galField gal = galField(p, 1);
	// Poly test1 = Poly(Int(1, p), 0);
	// std::cout<<"Factoring Polynomial "<<a<<std::endl;
	// std::cout<<std::endl;
	// auto fact = gal.factorization(a);
	// std::cout<<"            ";
	// for(auto&& irr: fact){
	//     std::cout<<"("<<irr.factor<<")^"<<irr.multiplicity<<" , ";
	//     test1 = test1 * mod_power(irr.factor, irr.multiplicity);
	// }
	// std::cout<<std::endl;
	// std::cout<<std::endl;
	// std::cout<<"            Test PASS: "<<(test1==a)<<std::endl;
	// std::cout<<std::endl;
	// clock_t end = clock();
	// double time = static_cast<double> (end-start) / CLOCKS_PER_SEC * 1000.0;
	// std::cout<<"            time exec.(ms): "<<time<<std::endl;
	// std::cout<<std::endl;
	//--------------------------------------------------

	//------------Berlekamp Factorization Test----------
	// std::cout<<"--------------------------"<<std::endl;
	// std::cout<<std::endl;
	// galField gal = galField(p, 1);	
	// auto fact = gal.berlekamp(a);
	// std::cout<<"Factoring Polynomial "<<a<<std::endl;
	// std::cout<<std::endl;
	// Poly test1 = Poly(Int(1, p), 0);
	// std::cout<<"            ";
	// for(auto&& irr: fact){
	//     std::cout<<"("<<irr<<")"<<", ";
	//     test1 = test1 * irr;
	// }
	// std::cout<<std::endl;
	// std::cout<<std::endl;
	// std::cout<<"            Test PASS: "<<(test1==a)<<std::endl;
	// std::cout<<std::endl;
	// clock_t end = clock();
	// double time = static_cast<double> (end-start) / CLOCKS_PER_SEC * 1000.0;
	// std::cout<<"            time exec.(ms): "<<time<<std::endl;
	// std::cout<<std::endl;
	//--------------------------------------------------
    }
}
