#include <stdio.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <deque>
#include <time.h>
#include <random>

#include "gal_field.h"
#include "matrix.h"

galField::galField(): dim(0), prim(0){
    mod_poli = Poly();
}

galField::galField(long long int p, int d): dim(d), prim(p){
    mod_poli = Poly();
}

galField::galField(long long int p, int d, Poly pol): dim(d), prim(p), mod_poli(pol){}

// // global constants
// const Poly galField::uno = Poly(Int(1, galField::prim), 0);
// const Poly galField::cero = Poly(Int(0, prim), 0);


// Irreducibility test
bool galField::isIrred(const Poly & f){
    Poly F = f;
    Poly uno = Poly(Int(1,prim),0);
    int n = F.getDeg();
    Poly x = Poly(Int(1,prim), 1);
    int pn = pow(prim, n);
    Poly xq;
    mod_power(x, pn, xq);
    xq = xq % F;

    if(xq != x){
	return false;
    }
    
    std::vector<int> primes = prime_divisors(n);
    
    for(auto&& i : primes){
	int aux = n / i;
	aux = pow(prim, aux);
	Poly xqn_aux;
	mod_power(x, aux, xqn_aux);
	xqn_aux = xqn_aux % F;
        xqn_aux = xqn_aux - x;	
	Poly mcd;
	gcd_poly(xqn_aux, F, mcd);

	if(mcd != uno){
	    return false;
	}
    }
    return true;
}

std::vector<Poly> galField::squarefree_factorization(Poly & f){
    Poly cero = Poly(Int(0,prim),0);
    Poly uno = Poly(Int(1,prim),0);
    int exp = pow(prim, dim-1);
    Poly f_diff = diff(f);

    if(f_diff == cero){
	Poly g;
	mod_power(f, exp, g);
	return fill(galField::squarefree_factorization(g),prim);
    }
    
    else{
	std::vector<Poly> result;
	Poly g0;
	gcd_poly(f, f_diff, g0);
	Poly w0 = f / g0;
	Poly w1 = w0;
	Poly g1;

	while(not(w1 == uno)){
	    gcd_poly(g0,w0, w1);
	    g1 = g0 / w1;
	    g0 = g1;
	    result.push_back(w0/w1);
	    w0 = w1;
	}

	if(g1 == uno){
	    return result;
	}
	
	else{
	    Poly g;
	    mod_power(g1, exp, g);
	    return (prod(result, fill(galField::squarefree_factorization(g), prim)));
	}
    }
}

// Distinct degree factorization
std::vector<Poly> galField::dd_factorization(Poly & f){
//    std::cout<<"DISTINCT DEGREE FACTORIZATION..."<<std::endl;
    std::vector<Poly> result;
    Poly F = f;
    Poly x = Poly(Int(1,prim),1);
    Poly uno = Poly(Int(1,prim),0);
    Poly h0 = Poly(Int(1,prim),1);
    Poly f0 = f;
    Poly fi = f;
    int q = pow(prim, dim);
    
    while(fi != uno){
	Poly hii;
	mod_power(h0, q, hii);
	Poly hi = hii % F;
	Poly aux = (hi)-(x);
	Poly gi;
	gcd_poly(aux, f0, gi);
	result.push_back(gi);
	fi = f0 / gi;
	f0 = fi;
	h0 = hi;
    }
    return result;
}

// Equal degree splitting; return a proper factor of f of degree d
Poly galField::eq_degree_splitting(Poly & f, int d){
//    std::cout<<"EQUAL DEGREE SPLITTING..."<<std::endl;
//    clock_t start = clock();
//    auto result = new Poly();
    Poly result;
    Poly F = f;
    Poly uno = Poly(Int(1,prim),0);
    Poly cero = Poly(Int(0,prim),0);    
    Poly a = rand_poly(F.getDeg(), prim);
    
    while(a.getDeg() < 1){
	a = rand_poly(F.getDeg(), prim);
    }

    bool T = true;
    while(T){
	Poly g1;
	gcd_poly(a, F, g1);
//	std::cout<<"splitt g1 = "<<g1<<std::endl;
	
	if(g1 != uno){
	    result = g1;
	    T = false;
	}

	else{
	    int exp = (pow(pow(prim,dim), d) -1) / 2;
	    Poly b;
	    mod_power(a, exp, b);
	    b = b % F;
	    Poly g2;
	    Poly aux = b - uno;
	    gcd_poly(aux, F, g2);
	    
	    if(g2 != cero && g2 != uno && g2 != F){
		result = g2;
		T = false;
	    }

	    else{
		a = rand_poly(F.getDeg(), prim);
		while(a.getDeg() < 1){
		    a = rand_poly(F.getDeg(), prim);
		}
	    }
	}
    }
     // std::cout<<"------------------TIME EQ DEG SPLITTING------------------"<<std::endl;
     // clock_t end = clock();
     // double time = static_cast<double> (end-start) / CLOCKS_PER_SEC * 1000.0;
     // std::cout<<"time exec.(ms): "<<time<<std::endl;
//    std::cout<<*result<<std::endl;
    return result;
}

// Equal degree factorization (algorithm procedure)
void galField::equal_degree_factorization_algorithm(Poly & f, int d, std::vector<Poly> &vec){
    Poly F = f;
//    auto result = new Poly();
    Poly result;
    
    if(f.getDeg() == d){
	result = F;
	vec.push_back(result);
    }
    
    else{
	Poly g = eq_degree_splitting(F, d);
//	std::cout<<g<<" , "<<d<<std::endl;
	equal_degree_factorization_algorithm(g, d, vec);
	Poly div = F / g;
	equal_degree_factorization_algorithm(div, d, vec);
    }
}

// Equal degree factorization (standard call)
std::vector<Poly> galField::equal_degree_factorization(Poly & f, int d){
//    std::cout<<"EQUAL DEGREE FACTORIZATION..."<<std::endl;
//    clock_t start = clock();
    std::vector<Poly> fact;
    equal_degree_factorization_algorithm(f, d, fact);

    // printing factorization at runtime
    // for(auto&& elem : fact){
    // 	std::cout<<elem<<" , ";
    // }
    
    // std::cout<<std::endl;
    // std::cout<<"------------------TIME EQUAL DEGREE------------------"<<std::endl;
    // clock_t end = clock();
    // double time = static_cast<double> (end-start) / CLOCKS_PER_SEC * 1000.0;
    // std::cout<<"time exec.(ms): "<<time<<std::endl;
    return fact;
}

// Factorization process
std::vector<pair> galField::factorization(Poly & f){
    // clock_t start = clock();
    // std::cout<<"-------------------CALL FACTORIZATION:------------------    f= "<<f<<std::endl;
    std::vector<pair> result;
    int q = pow(prim, dim);
    Poly F = f;
    
    if(F.lc().getValue() != 1){
	Poly leading_coef = Poly(F.lc(),0);
	F = F.normal();
	pair leading;
	leading.factor = leading_coef;
	leading.multiplicity = 1;
	result.push_back(leading);
    }


    Poly x = Poly(Int(1, prim), 1);
    Poly uno = Poly(Int(1, prim),0);
    Poly cero = Poly(Int(0, prim),0);
    Poly h0 = x;
    Poly v0 = F.normal();
    Poly vi = v0;
    Poly hi, aux, g;
    int i = 0;
    
    while(vi != uno){
	i++;
	mod_power(h0, q, hi);
	hi = hi % F;
	aux = hi - x;
	gcd_poly(aux, v0, g);
	
	if(g != uno){
//	    std::cout<<g<<" , i = "<<i<<std::endl;
	    std::vector<Poly> eq_deg_fact = equal_degree_factorization(g, i);
	    vi = v0;

	    for(auto&& fac : eq_deg_fact){
		int m = 0;
		while(vi % fac == cero){
		    vi = vi / fac;
		    m++;
		}

		pair add;
		add.factor = fac;
		add.multiplicity = m;

		if(m != 0){
		    result.push_back(add);
		}
	    }
	}
	h0 = hi;
	v0 = vi;
    }
    // std::cout<<"------------------TIME FACTORIZATION------------------"<<std::endl;
    // clock_t end = clock();
    // double time = static_cast<double> (end-start) / CLOCKS_PER_SEC * 1000.0;
    // std::cout<<"time exec.(ms): "<<time<<std::endl;
    return result;
}

// Berlekamp's Algorithm
void galField::berlekamp_splitting(Poly & f, std::vector<Poly> poly_basis, std::vector<pair> fact){
    Poly result;
    int q = pow(prim, dim);
    Poly uno = Poly(Int(1,prim),0);
//    std::uniform_int_distribution<> ran(0, prim - 1);

//    while(fact.size() < poly_basis.size()){
    
	// choose random combination from Berlekamp basis
	Poly v = Poly(Int(0,prim),0);
	for(int i = 0; i < poly_basis.size(); i++){
	    v = v + (Poly(Int(random_number(0, prim-1), prim), 0) * poly_basis[i]);
	}

	Poly gcd1;
	gcd_poly(v, f, gcd1);
	    
	if(gcd1 != uno && gcd1 != f){

	    Poly new_f = f / gcd1;

	    std::cout<<"new_f = "<<new_f<<std::endl;
	    std::cout<<"gcd1 = "<<gcd1<<std::endl;
	    
	    berlekamp_splitting(new_f, poly_basis, fact);
	    berlekamp_splitting(gcd1, poly_basis, fact);
	    
	    std::cout<<gcd1<<std::endl;
	    bool gcd1_founded_in_factorization = false;
	    for(int i = 0; i < fact.size(); i++){
		if(fact[i].factor == gcd1){
		    fact[i].multiplicity++;
		    gcd1_founded_in_factorization = true;
		}
	    }

	    if(!gcd1_founded_in_factorization){
		pair factor;
		factor.factor = gcd1;
		factor.multiplicity = 1;
		fact.push_back(factor);
	    }
	}

	else{
	    Poly b;
	    mod_power(v, pow(prim, dim)/2, b);
	    b = (b % f) - uno;
		
	    Poly gcd2;
	    gcd_poly(b, f, gcd2);
	    
	    if(gcd2 != uno && gcd2 != f){

		Poly new_f = f / gcd2;

		std::cout<<"new_f = "<<new_f<<std::endl;
		std::cout<<"gcd2 = "<<gcd2<<std::endl;
		
		berlekamp_splitting(new_f, poly_basis, fact);
		berlekamp_splitting(gcd2, poly_basis, fact);
		
		bool gcd2_founded_in_factorization = false;
		for(int i = 0; i < fact.size(); i++){
		    if(fact[i].factor == gcd2){
			fact[i].multiplicity++;
			gcd2_founded_in_factorization = true;
		    }
		}

		if(!gcd2_founded_in_factorization){
		    pair factor;
		    factor.factor = gcd2;
		    factor.multiplicity = 1;
		    fact.push_back(factor);
		}
	    }
	}
//    }
}

Poly berlekamp_factor(const Poly & f, std::vector<Poly> basis){
    Poly F = f;
    
    if(basis.size() == 1){
	return F;
    }
    
    else{
	int mod = f.first->getCoef().getMod();
	Poly uno = Poly(Int(1, mod),0);
	for(int i = 0; i < basis.size(); i++){
	    for(int p = 0; p < mod; p++){
		Poly gcd;
		Poly aux = basis[i] - Poly(Int(p, mod),0);
		gcd_poly(aux, F, gcd);
		if(gcd != uno && gcd != F){
		    std::cout<<gcd<<std::endl;
		    return gcd;
		}
	    }
	}
    }
}


// Berlekamp's factorization algorithm
std::vector<Poly> galField::berlekamp(Poly & f){
//    std::cout<<"-------------------CALL FACTORIZATION:------------------    f= "<<f<<std::endl;
    std::vector<pair> result;
    Poly F = f;
    std::vector<Poly> result2;
    // check if f is constant or ax^1
    if(F.getDeg() <= 1){
	pair factor_with_multip;
	factor_with_multip.factor = F;
	factor_with_multip.multiplicity = 1;
	result.push_back(factor_with_multip);
    }

    // if lenght is 1 f is a monom
    else if(F.first->sig == nullptr){
	pair factor_with_multip;
	factor_with_multip.factor = F;
	factor_with_multip.multiplicity = 1;
	result.push_back(factor_with_multip);
    }

    else{
	if(F.lc().getValue() != 1){
	    pair leading_coefficient;
	    leading_coefficient.factor = Poly(F.lc(),0);
	    leading_coefficient.multiplicity = 1;
	    result.push_back(leading_coefficient);
	    F = F.normal();
	}
	
	int n = F.getDeg();
	int q = pow(prim, dim);
	Poly uno = Poly(Int(1,prim),0);
	Poly x = Poly(Int(1,prim),1);
	Poly xq;
	mod_power(x, q, xq);
	xq = xq % F;
	Matrix Q = Matrix(n, n, prim);
    
	for(int i = 0; i < n; i++){
	    Poly xqi;
	    mod_power(x, q*i, xqi);
	    Poly xqi_test = mod_power(x, q*i) % F;
	    xqi = xqi % F;
	    // std::cout<<"xqi "<<xqi<<std::endl;
	    // std::cout<<"xqi test "<<xqi_test<<std::endl;	    
	    Q[i] = poliToVector_n(xqi, n);
	    // for(auto&& elem : Q[i]){
	    // 	std::cout<<elem<<", ";
	    // }
	    // std::cout<<std::endl;
	}
	
	Matrix M = Q - MatrixIn_mod(n, prim);
	M = M.transpose();
	std::vector<std::vector<Int>> basis = M.kernel_basis();
	
	// from basis to poly
	std::vector<Poly> poly_basis;
	for(int i = 0; i < basis.size(); i++){
	    Poly u;
	    for(int j = 0; j < basis[i].size(); j++){
		u = u + Poly(basis[i][j], j);
	    }
	    poly_basis.push_back(u);
	}

	Poly poly_to_factorize = F;
	
	int dim = poly_basis.size();

	result2.push_back(F);

	std::vector<Poly> factors;
	std::vector<Poly> irred;
	factors.push_back(F);
	
	while(result2.size() < dim){
//	    berlekamp_splitting(F, poly_basis, result);

	    // choose random combination from Berlekamp basis
	    //clock_t random_time_start = clock();
	    Poly v = Poly(Int(0,prim),0);
	    for(int i = 0; i < poly_basis.size(); i++){
	    	v = v + (Poly(Int(random_number(0, prim-1), prim), 0) * poly_basis[i]);
	    }
	    
	    // clock_t random_time_end = clock();
	    // double random_time = static_cast<double> (random_time_end - random_time_start)/CLOCKS_PER_SEC*1000.0;
	    // std::cout<<v<<" time spent in ms:  "<<random_time<<std::endl;
	    
	    std::vector<Poly> factors_new;

	    for(int i = 0; i < result2.size(); i++){
		Poly h = result2[i];
		Poly b = v % h;
		// std::cout<<"start mod power"<<std::endl;
		// clock_t time_start_mod_power = clock();
		Poly b_power;
		mod_power(b, (pow(prim, dim)-1)/2, b_power);
		// clock_t time_end_mod_power = clock();
		// double time_mod_power = static_cast<double> (time_end_mod_power - time_start_mod_power)/CLOCKS_PER_SEC*1000.0;
		// std::cout<<"modular power spent time in ms :  "<<time_mod_power<<std::endl;
		b = (b_power - uno) % h;
//		std::cout<<"power"<<b<<std::endl;
		Poly gcd;
		gcd_poly(b, h, gcd);
//		std::cout<<gcd<<std::endl;
		if(gcd == uno || gcd == h){
//		    std::cout<<"factoring: "<<h<<std::endl;
		    factors_new.push_back(h);
		}
		else{
		    Poly aux = h / gcd;
//		    std::cout<<"factoring: "<<aux<<std::endl;
		    factors_new.push_back(aux);
//		    std::cout<<"factoring: "<<gcd<<std::endl;
		    factors_new.push_back(gcd);
		}
	    }
	    result2 = factors_new;
	}
	
    }

    return result2;
}

// standard output printing
std::ostream & operator << (std::ostream &os, const galField & gal){
    os<<"Galois Field: "<<std::endl;
    os<<"        F"<<gal.prim<<"^"<<gal.dim<<" = Zp/("<<gal.mod_poli<<")"<<std::endl;
    return os;
}
