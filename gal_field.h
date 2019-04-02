#ifndef GAL_FIELD
#define GAL_FIELD

#include <stdio.h>
#include <iostream>
#include <vector>
#include "polynom.h"


// auxiliar structure for managing factors and their multiplicities
typedef struct{
    Poly factor;
    int multiplicity;
} pair;

// Actually works only like finite fields
class galField{
public:

    // Basic parameters:
    // mod_poli: modulus polynomial (Galois Extension Generator) NOT implemented yet
    // prim: prime modulus for integers
    // dim: dimension of the extension (must equal with mod_poli degree)
    Poly mod_poli;
    int dim;
    long long int prim;

    // Basic constructors
    galField();
    galField(long long int p, int d);
    galField(long long int p, int d, Poly pol);

    // Irreducibility test
    bool isIrred(const Poly & f);

    std::vector<Poly> squarefree_factorization(Poly & f);

    // Distinct degree factorization
    std::vector<Poly> dd_factorization(Poly & f);

    // Equal degree splitting
    Poly eq_degree_splitting(Poly & f, int d);

    // Equal degree factorization algorithm
    void equal_degree_factorization_algorithm(Poly & f, int d, std::vector<Poly>& vec);

    // Equal degree factorization (standard call)
    std::vector<Poly> equal_degree_factorization(Poly & f, int d);

    // Factorization process
    std::vector<pair> factorization(Poly & f);

    // Berlekamp's splitting algorithm
    void berlekamp_splitting(Poly & f, std::vector<Poly> poly_basis, std::vector<pair> fact);

    // Berlekamp's factorization algorithm
    std::vector<Poly> berlekamp(Poly & f);
    
};

// standard output printing
std::ostream & operator << (std::ostream &os, const galField & gal);


#endif
