#include "cedar/util/log.h"
#include <cedar/types.h>
#include <cedar/util/time_log.h>
#include <cedar/halo_exchanger_base.h>
#include <cedar/2d/relax/mpi/ml_shm.h>

extern "C" {

	void print_error(char *string)
	{
		cedar::log::error << string << std::endl;
	}

	void ftimer_begin(char *string)
	{
		cedar::timer_begin(string);
	}

	void ftimer_end(char *string)
	{
		cedar::timer_end(string);
	}

	void halo_exchange(int k, cedar::real_t *q, cedar::halo_exchanger_base *halof)
	{
		halof->exchange_func(k, q);
	}

	void halo_stencil_exchange(int k, cedar::real_t *so, cedar::halo_exchanger_base *halof)
	{
		halof->exchange_sten(k, so);
	}

	void ml_relax_shm_down(cedar::cdr2::mpi::ml_relax_pup * puper, cedar::real_t * rwork, int nlines)
	{
		puper->pack(rwork, nlines);
	}

	void ml_relax_shm_up(cedar::cdr2::mpi::ml_relax_pup * puper, cedar::real_t * rwork, int nlines)
	{
		puper->unpack(rwork, nlines);
	}
}
