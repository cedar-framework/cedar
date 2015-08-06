#include <sys/types.h>
#include <sys/time.h>
#include <iostream>
#include <stdlib.h>
#include <mm_malloc.h>

#include <boxmg/2d/core>


double dtime()
{
	struct timeval tp;
	struct timezone tzp;

	gettimeofday(&tp,&tzp);
	return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}


void time_raw(int size)
{
	int i,j,k;
	double timer;
	int ntests = 10;
	using namespace boxmg;
	using namespace boxmg::bmg2d::core;

	data_t *a = (data_t*) _mm_malloc(size*size*5*sizeof(data_t), sizeof(data_t));
	data_t *v = (data_t*) _mm_malloc(size*size*sizeof(data_t), sizeof(data_t));

	std::cout << "Timing raw array" << std::endl;
	std::cout << "================" << std::endl;

	timer = dtime();
	for (k = 0; k < ntests; k++) {
		for (i = 0; i < size; i++) {
			for (j = 0; j < size; j++) {
				a[(i*size + j)*5] = rand();
				a[(i*size + j)*5 + 1] = rand();
				a[(i*size + j)*5 + 2] = rand();
				a[(i*size + j)*5 + 3] = rand();
				a[(i*size + j)*5 + 4] = rand();
			}
		}
	}
	std::cout << "Write: " << (dtime() - timer) / ntests << std::endl;

	timer = dtime();
	for (k = 0; k < ntests; k++) {

		for (i = 0; i < size; i++) {
			for (j = 0; j < size; j++) {
				v[i*size+j] = a[(i*size+j)*5] - .25*(a[(i*size+j)*5+1]+
													 a[(i*size+j)*5+2]+
													 a[(i*size+j)*5+3]+
													 a[(i*size+j)*5+4]);
			}
		}
	}

	std::cout << "Read: " << (dtime() - timer) / ntests << std::endl << std::endl;

	_mm_free(a);
	_mm_free(v);
}

void time_grid(int size)
{
	int i,j,k;
	double timer;
	int ntests = 10;
	std::cout << "Timing Grid Stencil data structure" << std::endl;
	std::cout << "==================================" << std::endl;
	using namespace boxmg::bmg2d::core;

	GridStencil gs(size,size);
	GridFunc gf(size,size);

	timer = dtime();
	for (k = 0; k < ntests; k++) {
		for (i = 0; i < size; i++) {
			for (j = 0; j < size; j++) {
				gs(i,j,Dir::C) = rand();
				gs(i,j,Dir::W) = rand();
				gs(i,j,Dir::S) = rand();
				gs(i,j,Dir::SW) = rand();
				gs(i,j,Dir::SE) = rand();
			}
		}
	}
	std::cout << "Write: " << (dtime() - timer) / ntests << std::endl;

	timer = dtime();
	for (k = 0; k < ntests; k++) {

		for (i = 0; i < size; i++) {
			for (j = 0; j < size; j++) {
				gf(i,j) = gs(i,j,Dir::C) - .25*(gs(i,j,Dir::W) + gs(i,j,Dir::S) +
												gs(i,j,Dir::SW) + gs(i,j,Dir::SE));
			}
		}
	}

	std::cout << "Read: " << (dtime() - timer) / ntests << std::endl << std::endl;
}


int main(int argc, char *argv[])
{
	int size = 2000;

	using namespace boxmg::bmg2d::core;

	srand(time(NULL));
	time_grid(size);
	time_raw(size);

	return 0;
}
