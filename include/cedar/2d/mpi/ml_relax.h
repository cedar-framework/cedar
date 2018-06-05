#ifndef CEDAR_2D_KERNEL_RELAX_MPI_ML_H
#define CEDAR_2D_KERNEL_RELAX_MPI_ML_H

#include <stdbool.h>

#include <cedar/2d/ftn/mpi/BMG_parameters_c.h>
#include <cedar/2d/ftn/mpi/BMG_workspace_c.h>
#include <cedar/kernels/line_relax.h>
#include <cedar/2d/mpi/types.h>


extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_SETUP_trisolve(int nproc, len_t nlines, int min_gsz, int nog, int *nspace,
	                               int *nol);
	void BMG2_SymStd_SETUP_ADD_SOR_PTRS(int nproc, len_t nlines, int min_gsz, int nol, int *poffset,
	                                    int *psor_lev);
	void MPI_BMG2_SymStd_SETUP_comm(MPI_Fint *comm, int nol, int rank, int min_gsz);
	void MPI_BMG2_SymStd_SETUP_lines_x(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl);
	void MPI_BMG2_SymStd_SETUP_lines_y(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl);
	void MPI_BMG2_SymStd_relax_lines_x_ml(int k, real_t *so, real_t *qf, real_t *q,
	                                      real_t *sor, real_t *res,
	                                      len_t II, len_t JJ, len_t igs, len_t jgs,
	                                      int nog, int nstencil, int irelax_sym, int updown,
	                                      len_t *datadist, real_t *rwork, len_t nmsgr,
	                                      MPI_Fint mpicomm, MPI_Fint *xcomm, int nolx,
	                                      int *tdsx_sor_ptrs, len_t nsor, real_t *tdg,
	                                      bool *fact_flags, bool factorize, void *halof);
	void MPI_BMG2_SymStd_relax_lines_y_ml(int k, real_t *so, real_t *qf, real_t *q,
	                                      real_t *sor, real_t *B,
	                                      len_t II, len_t JJ, len_t igs, len_t jgs,
	                                      int nog, int nstencil, int irelax_sym, int updown,
	                                      len_t *datadist, real_t *rwork, len_t nmsgr,
	                                      MPI_Fint mpicomm, MPI_Fint *ycomm, int noly,
	                                      int *tdsy_sor_ptrs, len_t nsor, real_t *tdg,
	                                      bool *fact_flags, bool factorize, void *halof);
}

namespace cedar { namespace cdr2 { namespace mpi {

template<relax_dir dir>
class ml_line_relax : public kernels::line_relax<stypes, dir>
{
public:
	ml_line_relax(bool fact=true) : factorize(fact), initialized(false) {}
	~ml_line_relax()
	{
		if (initialized) {
			for (std::size_t i = 0; i < fact_flags.size(); i++)
				delete[] fact_flags[i];
		}
	}

	void set_factorize(bool fact) { this->factorize = fact; }

	void init_ndim(int nproc, int nproc_other, int coord, int coord_other, len_t nlines, len_t nlines_other,
	               int min_gsz, int nog, MPI_Comm comm)
	{
			int workspace_size, nspace;
			MPI_BMG2_SymStd_SETUP_trisolve(nproc, nlines, min_gsz, nog,
			                               &nspace, &this->nlevels);
			comms.init(2, this->nlevels - 1);
			comms.set(MPI_COMM_NULL);
			MPI_Comm gcomm;
			MPI_Comm_split(comm, coord_other, coord, &gcomm);
			comms(0, this->nlevels - 2) = MPI_Comm_c2f(gcomm);
			MPI_BMG2_SymStd_SETUP_comm(comms.data(), this->nlevels, coord+1, min_gsz);
			sor_ptrs.emplace_back(this->nlevels);
			BMG2_SymStd_SETUP_ADD_SOR_PTRS(nproc, nlines, min_gsz, this->nlevels,
			                               &workspace_size, sor_ptrs.back().data());
			workspace_size += (nlines + 2)*(nlines_other + 2)*4;
			tdg.emplace_back(workspace_size);
			// this bound may not be correct!
			int rwork_size = std::max(8 * nlines_other * nproc + 8 * nlines,
			                          8 * nlines * nproc_other + 8 * nlines_other);
			rwork.emplace_back(rwork_size);
			bool *fflags;
			fflags = new bool[2 * nog];
			initialized = true;
			for (auto i : range<std::size_t>(2 * nog))
				fflags[i] = true;
			fact_flags.push_back(fflags);
	}

	virtual void setup(const stencil_op<five_pt> & so,
	                   relax_stencil & sor) override
	{
		this->setup_impl(so, sor);
	}
	virtual void setup(const stencil_op<nine_pt> & so,
	                   relax_stencil & sor) override
	{
		this->setup_impl(so, sor);
	}
	virtual void run(const stencil_op<five_pt> & so,
	                 grid_func & x,
	                 const grid_func & b,
	                 const relax_stencil & sor,
	                 grid_func & res,
	                 cycle::Dir cdir) override
	{
		this->run_impl(so, x, b, sor, res, cdir);
	}
	virtual void run(const stencil_op<nine_pt> & so,
	                 grid_func & x,
	                 const grid_func & b,
	                 const relax_stencil & sor,
	                 grid_func & res,
	                 cycle::Dir cdir) override
	{
		this->run_impl(so, x, b, sor, res, cdir);
	}

	template<class sten>
	void setup_impl(const stencil_op<sten> & so,
	                relax_stencil & sor)
	{
		int nstencil = stencil_ndirs<sten>::value;
		auto & sod = const_cast<stencil_op<sten>&>(so);

		auto & topo = sod.grid();
		int min_gsz = this->params->ml_relax.min_gsz;

		if (dir == relax_dir::x) {
			init_ndim(topo.nproc(0), topo.nproc(1),
			          topo.coord(0), topo.coord(1),
			          topo.nlocal(0) - 2, topo.nlocal(1) - 2,
			          min_gsz, topo.nlevel(), topo.comm);
			MPI_BMG2_SymStd_SETUP_lines_x(sod.data(), sor.data(), so.len(0), so.len(1), nstencil);
		} else {
			init_ndim(topo.nproc(1), topo.nproc(0),
			          topo.coord(1), topo.coord(0),
			          topo.nlocal(1) - 2, topo.nlocal(0) - 2,
			          min_gsz, topo.nlevel(), topo.comm);
			MPI_BMG2_SymStd_SETUP_lines_y(sod.data(), sor.data(), so.len(0), so.len(1), nstencil);
		}
	}

	template<class sten>
	void run_impl(const stencil_op<sten> & so,
	              grid_func & x,
	              const grid_func & b,
	              const relax_stencil & sor,
	              grid_func & res,
	              cycle::Dir cycle_dir)
	{
		int k, kf;
		int updown, nstencil;

		auto & sod = const_cast<stencil_op<sten>&>(so);
		grid_topo & topo = sod.grid();
		relax_stencil & sord = const_cast<relax_stencil&>(sor);
		auto & bd = const_cast<grid_func&>(b);

		k = topo.level()+1;
		kf = topo.nlevel();

		nstencil = stencil_ndirs<sten>::value;

		if (cycle_dir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		// ibc = BMG_BCs_definite;
		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		int dim;
		if (dir == relax_dir::x)
			dim = 0;
		else
			dim = 1;
		len_t * datadist = this->halof->datadist(dim, k);

		std::function<void(int k, real_t *so, real_t *qf, real_t *q,
		                   real_t *sor, real_t *res,
		                   len_t II, len_t JJ, len_t igs, len_t jgs,
		                   int nog, int nstencil, int irelax_sym, int updown,
		                   len_t *datadist, real_t *rwork, len_t nmsgr,
		                   MPI_Fint mpicomm, MPI_Fint *xcomm, int nolx,
		                   int *tdsx_sor_ptrs, len_t nsor, real_t *tdg,
		                   bool *fact_flags, bool factorize, void *halof)> relax_lines;

		if (dir == relax_dir::x)
			relax_lines = MPI_BMG2_SymStd_relax_lines_x_ml;
		else
			relax_lines = MPI_BMG2_SymStd_relax_lines_y_ml;

		relax_lines(k, sod.data(), bd.data(), x.data(), sord.data(), res.data(),
		            so.len(0), so.len(1), topo.is(0), topo.is(1),
		            kf, nstencil, BMG_RELAX_SYM, updown,
		            datadist, rwork[kf-k].data(), rwork[kf-k].len(0),
		            fcomm, this->comms.data(), this->nlevels,
		            sor_ptrs[kf-k].data(), tdg[kf-k].len(0), tdg[kf-k].data(),
		            fact_flags[kf-k], this->factorize, this->halof);
	}

protected:
	int nlevels;                            /** Number of levels in tridiagonal solve */
	std::vector<array<int, 1>> sor_ptrs;    /** Pointers into tridiagonal workspace */
	std::vector<bool*> fact_flags;          /** flags for factorization */
	std::vector<array<real_t, 1>> tdg;      /** Tridiagonal solver workspace array */
	std::vector<array<real_t, 1>> rwork;    /** Buffer workspace for tridiagonal solver */
	array<MPI_Fint, 2> comms;               /** Communicators for ml lines */
	bool factorize;                         /** Flag for local factorization (contrast to elimination) */
	bool initialized;
};

}}}

#endif
