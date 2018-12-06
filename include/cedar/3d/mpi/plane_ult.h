#ifndef CEDAR_3D_MPI_PLANE_ULT_H
#define CEDAR_3D_MPI_PLANE_ULT_H

#ifdef PLANE_AGG
#include <abt.h>
#endif
#include <vector>
#include <memory>

#include <cedar/2d/mpi/solver.h>

namespace cedar { namespace cdr3 { namespace mpi {

template<class sten2>
struct ult_params
{
	cdr2::mpi::solver<sten2> *slv;
};


void plane_ult_run_full(void *args);
void plane_ult_run_comp(void *args);

#ifdef PLANE_AGG
template<class sten2>
class plane_ult
{
public:
	plane_ult() : initialized(false) {}

	void init()
	{
		if (ABT_initialized() == ABT_ERR_UNINITIALIZED) {
			ABT_init(0, NULL);
		}

		if (not initialized) {
			ABT_xstream_self(&xstream);
			ABT_xstream_get_main_pools(xstream, 1, &pool);
		}

		initialized = true;
	}


	~plane_ult()
	{
		if (initialized) {
			for (std::size_t i = 0; i < threads.size(); i++) {
				if (thread_created[i])
					ABT_thread_free(&threads[i]);
			}
		}
	}


	void add_plane(cdr2::mpi::solver<sten2> *slv)
	{
		if (not initialized)
			init();
		params.emplace_back();
		thread_created.emplace_back(false);
		threads.emplace_back();
		auto & param = params.back();
		param.slv = slv;
	}


	void start(std::size_t i)
	{
		void (*compute)(void*);
		if (std::is_same<sten2, cdr2::five_pt>::value)
			compute = plane_ult_run_comp;
		else
			compute = plane_ult_run_full;

		if (not thread_created[i]) {
			ABT_thread_create(pool, compute, (void*) &params[i],
			                  ABT_THREAD_ATTR_NULL, &threads[i]);
		} else {
			ABT_thread_revive(pool, compute, (void*) &params[i],
			                  &threads[i]);
		}
	}


	void join(std::size_t i)
	{
		ABT_thread_join(threads[i]);
	}


	std::vector<ABT_thread> * get_threads() { return &threads; }

protected:
	bool initialized;
	ABT_pool pool;
	ABT_xstream xstream;
	std::vector<bool> thread_created;
	std::vector<ABT_thread> threads;
	std::vector<ult_params<sten2>> params;
	bool threads_created;
};
#else
template<class sten2> using plane_ult = int;
#endif

}}}

#endif
