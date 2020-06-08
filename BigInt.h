#ifndef BIGINT_H
#define BIGINT_H

#include <cmath>
#include <iostream>
#include <map>

class BigInt{
	private:
		std::map<int, int> data;

		std::map<int, int> prime_factorize(int n);

	public:
		BigInt();
		BigInt(int n);

		BigInt& operator*=(int rhs);
		BigInt& operator*=(const BigInt& rhs);

		BigInt& operator/=(int rhs);
		BigInt& operator/=(const BigInt& rhs);

		void print();
		double to_double();
};

#endif

