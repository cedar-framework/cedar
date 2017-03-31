#ifndef CEDAR_DECOMP_H
#define CEDAR_DECOMP_H

#include <array>
#include <vector>
#include <algorithm>

#include <cedar/types.h>

namespace cedar
{


inline std::vector<int> sieve(int n)
{
	std::vector<bool> prime(n+1);
	std::vector<int> ret;
	std::fill(prime.begin(), prime.end(), true);

	for (int p = 2; p*p <= n; p++) {
		if (prime[p]) {
			for (int i = 2*p; i <= n; i+= p) {
				prime[i] = false;
			}
		}
	}

	for (int p = 2; p <= n; p++) {
		if (prime[p])
			ret.push_back(p);
	}

	return ret;
}


inline std::vector<int> factor(int np)
{
	std::vector<int> facs;
	std::vector<int> primes = sieve(np);

	int p = 0;
	for (unsigned int i = 0; i < primes.size(); i++) {
		if (np % primes[p] == 0) {
			facs.push_back(primes[p]);
			np = np / primes[p];
		} else {
			p++;
		}
	}

	std::sort(facs.begin(), facs.end());
	std::reverse(facs.begin(), facs.end());

	return facs;
}


template<int ND>
std::array<int, ND> grid_decomp(std::array<len_t,ND> n, int np)
{
	std::array<int,ND> decomp;

	for (int i = 0; i < ND; i++) {
		decomp[i] = 1;
	}

	std::vector<int> facs = factor(np);

	len_t curr;
	int ind;
	for (auto fac : facs) {
		curr = 0;
		ind = 0;
		for (int j = 0; j < ND; j++) {
			if (n[j] > curr) {
				curr = n[j];
				ind = j;
			}
		}
		decomp[ind] *= fac;
		n[ind] = n[ind] / fac;
	}

	return decomp;
}

}
#endif
