# main: main.cpp
# 	clang++ -g -std=gnu++11 int.cpp matrix.cpp monom.cpp polynom.cpp gal_field.cpp main.cpp -o main.out
main: main.cpp
	g++ -g -O2 -std=gnu++11 int.cpp monom.cpp matrix.cpp polynom.cpp gal_field.cpp main.cpp -o main.out
