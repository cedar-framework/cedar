#ifndef CEDAR_3D_MPI_MATVEC_H
#define CEDAR_3D_MPI_MATVEC_H

#include <cedar/3d/mpi/types.h>
#include <cedar/kernels/matvec.h>

extern "C" {
	using namespace cedar;
	void MPI_BMG3_SymStd_UTILS_matvec(int kg, real_t *so, real_t *qf,
	                                  real_t *q, len_t ii, len_t jj, len_t kk,
	                                  int nog, int ifd, int nstncl, void *halof);
}

namespace cedar { namespace cdr3 { namespace mpi {

class matvec_f90 : public kernels::matvec<stypes>
{
	void run(const stencil_op<seven_pt> & so,
	         const grid_func & x,
	         grid_func & y) override
	{ this->run_impl(so, x, y); }
	void run(const stencil_op<xxvii_pt> & so,
	         const grid_func & x,
	         grid_func & y) override
	{ this->run_impl(so, x, y); }

	template<class sten>
	void run_impl(const stencil_op<sten> & so,
	              const grid_func & x,
	              grid_func & y)
	{
		int kg, ifd, nstencil, nog;

		auto & sod = const_cast<stencil_op<sten>&>(so);
		grid_func & xd = const_cast<grid_func&>(x);
		grid_topo & topo = sod.grid();

		nog = topo.nlevel();
		kg = topo.level()+1;
		nstencil = stencil_ndirs<sten>::value;
		if (std::is_same<sten, seven_pt>::value)
			ifd = 1;
		else
			ifd = 0;

		MPI_BMG3_SymStd_UTILS_matvec(kg, sod.data(), y.data(), xd.data(),
		                             so.len(0), so.len(1), so.len(2),
		                             nog, ifd, nstencil, halof);
	}
};

}}}
#endif
