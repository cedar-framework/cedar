#ifndef BOXMG_DECOMP_H
#define BOXMG_DECOMP_H

#include <array>
#include <vector>
#include <algorithm>

#include <boxmg/types.h>

namespace boxmg
{

// placeholder for future factorization function
std::vector<int> factor(int np)
{
	std::vector<int> facs;
	std::vector<int> primes({2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503});

	int p = 0;
	for (int i = 0; i < primes.size(); i++) {
		if (np % primes[p] == 0) {
			facs.push_back(primes[p]);
			np = np / primes[p];
		} else {
			p++;
		}
	}

	std::sort(facs.begin(), facs.end());

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

	int curr,ind;
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
