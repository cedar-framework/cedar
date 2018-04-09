#ifndef CEDAR_KERNEL_MPI_REGISTRY_H
#define CEDAR_KERNEL_MPI_REGISTRY_H

#include <cedar/kernels/residual.h>
#include <cedar/2d/mpi/types.h>

extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_residual(int k, int kf, int nog,
	                              real_t *SO, real_t *QF, real_t *Q, real_t *RES,
	                              len_t ii, len_t jj, int ifd, int nstencil,
	                              int irelax, int irelax_sym,
	                              int mpicomm);
}

namespace cedar { namespace cdr2 { namespace mpi {

class residual_f90 : public kernels::residual<stypes>
{
	void run(const stencil_op<five_pt> & so,
	         const grid_func & x,
	         const grid_func & b,
	         grid_func & r) override
	{
		this->run_impl(so, x, b, r);
	}
	void run(const stencil_op<nine_pt> & so,
	         const grid_func & x,
	         const grid_func & b,
	         grid_func & r) override
	{
		this->run_impl(so, x, b, r);
	}

	template<class sten>
	void run_impl(const stencil_op<sten> & so,
	              const grid_func & x,
	              const grid_func & b,
	              grid_func & r)
	{
		int k, kf, nog, ifd, nstencil;
		auto & Ad = const_cast<stencil_op<sten> &>(so);
		auto & xd = const_cast<grid_func&>(x);
		auto & bd = const_cast<grid_func&>(b);
		grid_topo & topo = Ad.grid();

		nstencil = stencil_ndirs<sten>::value;
		if (std::is_same<five_pt, sten>::value)
			ifd = 1;
		else
			ifd = 0;

		int irelax = 0;
		int irelax_sym = 0;

		nog = kf = topo.nlevel();
		k = topo.level()+1;

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);
		MPI_BMG2_SymStd_residual(k, kf, nog,
		                         Ad.data(), bd.data(), xd.data(), r.data(),
		                         r.len(0), r.len(1), ifd, nstencil,
		                         irelax, irelax_sym, fcomm);
		halof->exchange_func(k, r.data());
	}
};


}}}
#endif
