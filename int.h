#ifndef INT_H
#define INT_H

#include <stdio.h>
#include <iostream>

#include <vector>


class Int{
private:
    // Basic parameters:
    // val : integer value of Int
    // mod : modulus of the integer
    long long signed int val;
    long long signed int mod;

    // avoid mod < 0
    bool check();
	
public:
	
    class Invalid{};
	
    Int();

    // explicit constructor for prevent implicit conversion
    explicit Int(long long signed int m);
    Int(long long signed int m, long long signed int p);
    ~Int() = default;
    long long signed int getValue() const{return val;};
    long long signed int getMod() const{return mod;};

    // opposite of the integer % mod
    void setNeg();

    // copy constructor
    Int(const Int & a);

    // assignment operator = default
    Int & operator=(const Int & a);

    // set functions
    void set_value(const long long signed int a);
    void set_mod(const long long signed int p);

    bool operator==(const Int & a) const;
    bool operator!=(const Int & a) const;
    
    // for uses like Int(1,7) == 0
    bool operator==(long long signed int a) const ;
    bool operator!=(long long signed int a) const;
    
    Int operator+(const Int & a) const;
    Int operator-(const Int & a) const;
    Int operator*(const Int & a) const;
    Int operator*(long long signed int a) const;
    // return inverse multiplication
    Int operator/(const Int & a) const;
    // repeated squaring algorithm for power
    Int operator^(const long long signed int n);
};

// Procedural functions for basic operations
void sum(const Int & a, const Int & b, Int & result);
void sub(const Int & a, const Int & b, Int & result);
void mult(const Int & a, const Int & b, Int & result);
void power(const Int & a, const long long int n, Int & result);
void inverse_mod(const Int & a, Int & result);
void div(const Int & a, const Int & b, Int & result);

int random_number(int low, int high);

// Extended Euclides Algorithm for integers
void extEuclides(long long signed int a, long long signed int b, long long signed int auxA, long long signed int auxB, long long signed int mcd);

// modular inverse using Extended Euclides Algorithm
long long signed int invMult(const long long signed int a, const long long signed int mod);

// return the binary representation of an integer
std::vector<int> intToBin(const long long signed int a);

long long signed int power(long long signed int a, long long signed int exp);

// Chinesse Remainder Theorem
long long signed int garnerRemainder(long long signed int u[], long long signed int m[], long long signed int l);

std::vector<int> prime_divisors(const long long signed int n);

// for printing in standard output
std::ostream & operator <<(std::ostream &os, const Int & a);

#endif
