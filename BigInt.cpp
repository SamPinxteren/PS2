#include "BigInt.h"

std::map<int, int> BigInt::prime_factorize(int n){
	std::map<int, int> result;

	int increments[] = {4, 2, 4, 2, 4, 6, 2, 6};

	while (n % 2 == 0){
		++result[2];
		n /= 2;
	}

	while (n % 3 == 0){
		++result[3];
		n /= 3;
	}

	while (n % 5 == 0){
		++result[5];
		n /= 5;
	}

	int f = 7;
	int i = 0;
	while (f * f <= n){
		if (n % f == 0){
			++result[f];
			n /= f;
		} else {
			f += increments[i % 8];
		}
	}

	if (n != 1) ++result[n];

	return result;
}

BigInt::BigInt(){
}

BigInt::BigInt(int n){
	data = prime_factorize(n);
}

BigInt& BigInt::operator*=(int rhs){
	(*this) *= BigInt(rhs);
	return *this;
}

BigInt& BigInt::operator*=(const BigInt& rhs){
	for (auto const& x: rhs.data){
		data[x.first] += x.second;
	}
	return *this;
}

void BigInt::print(){
	std::cout << "=== BigInt ===" << std::endl;
	for (auto const& x: data){
		std::cout << x.first << " " << x.second << std::endl;
	}
	std::cout << "==============" << std::endl;
}

BigInt& BigInt::operator/=(int rhs){
	(*this) /= BigInt(rhs);
	return *this;
}

BigInt& BigInt::operator/=(const BigInt& rhs){
	for (auto const& x: rhs.data){
		data[x.first] -= x.second;
	}
	return *this;
}

double BigInt::to_double(){
	double result = 1;
	for (auto const& x: data){
		result *= pow(x.first, x.second);
	}
	return result;
}
