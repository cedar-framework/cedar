#ifndef CEDAR_2D_KERNEL_RELAX_MPI_ML_H
#define CEDAR_2D_KERNEL_RELAX_MPI_ML_H

#include <stdbool.h>

#include "cedar/2d/ftn/BMG_parameters_c.h"
#include "cedar/2d/ftn/mpi/BMG_workspace_c.h"
#include <cedar/2d/mpi/msg_exchanger.h>
#include <cedar/2d/mpi/tausch_exchanger.h>
#include <cedar/2d/relax/mpi/ml_shm.h>


extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_SETUP_trisolve(int nproc, len_t nlines, int min_gsz, int nog, int *nspace,
	                               int *nol);
	void BMG2_SymStd_SETUP_ADD_SOR_PTRS(int nproc, len_t nlines, int min_gsz, int nol, int *poffset,
	                                    int *psor_lev, bool shm_enabled, int node_size, int shm_size);
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
	                                      bool *fact_flags, bool shm_enabled,
	                                      void *puper, void *halof);
}

namespace cedar { namespace cdr2 { namespace mpi {


enum class relax_dir { x, y };

template<class halo_exchanger>
	len_t * get_datadist(halo_exchanger * halof, int k, relax_dir dir);

template<>
	len_t * get_datadist(mpi::msg_exchanger * halof, int k, relax_dir dir)
{
	len_t * ret;
	kernel::impls::MsgCtx * ctx = (kernel::impls::MsgCtx*) halof->context_ptr();

	if (dir == relax_dir::x)
		ret = &ctx->msg_geom.data()[ctx->pLS(ipL_LS_XDataDist,k-1)-1];
	else
		ret = &ctx->msg_geom.data()[ctx->pLS(ipL_LS_YDataDist,k-1)-1];

	return ret;
}

template<>
	len_t * get_datadist(mpi::tausch_exchanger * halof, int k, relax_dir dir)
{
	len_t * ret;
	mpi::line_pkg & line_info = halof->line_data;

	if (dir == relax_dir::x)
		ret = line_info.datadist[0].data();
	else
		ret = line_info.datadist[1].data();

	return ret;
}


class ml_line_relaxer
{
public:
	ml_line_relaxer(relax_dir dir) : dir(dir), initialized(false), shm(false) {}
	~ml_line_relaxer()
	{
		if (initialized) {
			for (std::size_t i = 0; i < fact_flags.size(); i++)
				delete[] fact_flags[i];
		}
	}

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
			                               &workspace_size, sor_ptrs.back().data(), shm, 0, 0);
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


	void init_comms_shm(MPI_Comm comm, int coord, int coord_other, len_t nlines,
	                    int min_gsz, int nog)
	{
		MPI_Comm linecomm; // Fine level communicator for line relaxation
		MPI_Comm_split(comm, coord_other, coord, &linecomm);
		MPI_Comm_split_type(linecomm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shm_comm);

		int shm_size, shm_rank;
		MPI_Comm_rank(shm_comm, &shm_rank);
		MPI_Comm_size(shm_comm, &shm_size);

		{ // split limcomm into roots of shm_comm (node level communicator)
			int linecomm_rank;

			MPI_Comm_rank(linecomm, &linecomm_rank);
			int color = shm_rank;
			MPI_Comm_split(linecomm, color, linecomm_rank, &node_comm);
		}

		int node_size;
		MPI_Comm_size(node_comm, &node_size);

		int nspace;
		MPI_BMG2_SymStd_SETUP_trisolve(node_size, nlines, min_gsz, nog,
		                               &nspace, &this->nlevels);
		// add space for shm level
		nspace += (2 * shm_size + 2) * nlines * 4;
		// add shared memory level
		if (this->nlevels > 2)
			this->nlevels++;

		comms.init(2, this->nlevels - 1);
		comms.set(MPI_COMM_NULL);
		comms(0, this->nlevels - 2) = MPI_Comm_c2f(linecomm);
		comms(1, this->nlevels - 2) = MPI_Comm_c2f(shm_comm);
		if (this->nlevels > 2)
			comms(0, this->nlevels - 3) = MPI_Comm_c2f(node_comm);

		if (shm_rank == 0) {
			int myproc;
			MPI_Comm_rank(node_comm, &myproc);
			if (this->nlevels > 2)
				MPI_BMG2_SymStd_SETUP_comm(comms.data(), this->nlevels-1, myproc+1, min_gsz);
		}
	}


	void init_ndim_shm(int nproc, int nproc_other, int coord, int coord_other,
	                   len_t nlines, len_t nlines_other,
	                   int min_gsz, int nog, MPI_Comm comm)
	{
		int shm_size;
		if (not initialized) {
			init_comms_shm(comm, coord, coord_other, nlines, min_gsz, nog);
			int shm_len = nlines * 8;
			MPI_Win_allocate_shared(shm_len*sizeof(real_t), 1, MPI_INFO_NULL, shm_comm, &shm_buff, &shm_win);
			puper.init(shm_comm, node_comm, shm_win, shm_buff);
		}

		int node_count;
		MPI_Comm_size(node_comm, &node_count);
		MPI_Comm_size(shm_comm, &shm_size);

		int workspace_size;
		sor_ptrs.emplace_back(this->nlevels);
		BMG2_SymStd_SETUP_ADD_SOR_PTRS(node_count, nlines, min_gsz, this->nlevels,
		                               &workspace_size, sor_ptrs.back().data(),
		                               this->shm, node_count, shm_size);
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


	template<class sten>
		void setup(const kernel_params & params,
		           const mpi::stencil_op<sten> & so,
		           relax_stencil & sor)
	{
		int nstencil = stencil_ndirs<sten>::value;
		auto & sod = const_cast<mpi::stencil_op<sten>&>(so);

		auto & topo = sod.grid();
		int min_gsz = params.ml_relax.min_gsz;
		this->shm = params.ml_relax.shm;

		if (dir == relax_dir::x) {
			if (shm) {
				init_ndim_shm(topo.nproc(0), topo.nproc(1),
				              topo.coord(0), topo.coord(1),
				              topo.nlocal(0) - 2, topo.nlocal(1) - 2,
				              min_gsz, topo.nlevel(), topo.comm);
			} else {
				init_ndim(topo.nproc(0), topo.nproc(1),
				          topo.coord(0), topo.coord(1),
				          topo.nlocal(0) - 2, topo.nlocal(1) - 2,
				          min_gsz, topo.nlevel(), topo.comm);
			}
			MPI_BMG2_SymStd_SETUP_lines_x(sod.data(), sor.data(), so.len(0), so.len(1), nstencil);
		} else {
			init_ndim(topo.nproc(1), topo.nproc(0),
			          topo.coord(1), topo.coord(0),
			          topo.nlocal(1) - 2, topo.nlocal(0) - 2,
			          min_gsz, topo.nlevel(), topo.comm);
			MPI_BMG2_SymStd_SETUP_lines_y(sod.data(), sor.data(), so.len(0), so.len(1), nstencil);
		}
	}


	template<class sten, class halo_exchanger>
		void solve(const kernel_params & params,
		           halo_exchanger * halof,
		           const mpi::stencil_op<sten> & so,
		           mpi::grid_func & x,
		           const mpi::grid_func & b,
		           const relax_stencil & sor,
		           mpi::grid_func & res,
		           cycle::Dir cycle_dir)
	{
		using namespace cedar::cdr2;
		int k, kf;
		int updown, nstencil;

		auto & sod = const_cast<mpi::stencil_op<sten>&>(so);
		grid_topo & topo = sod.grid();
		relax_stencil & sord = const_cast<relax_stencil&>(sor);
		mpi::grid_func & bd = const_cast<mpi::grid_func&>(b);

		k = topo.level()+1;
		kf = topo.nlevel();

		nstencil = stencil_ndirs<sten>::value;

		if (cycle_dir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		// ibc = BMG_BCs_definite;
		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);
		MPI_Fint fwin = MPI_Win_c2f(shm_win);

		len_t * datadist = get_datadist<halo_exchanger>(halof, k, dir);

		std::function<void(int k, real_t *so, real_t *qf, real_t *q,
		                   real_t *sor, real_t *res,
		                   len_t II, len_t JJ, len_t igs, len_t jgs,
		                   int nog, int nstencil, int irelax_sym, int updown,
		                   len_t *datadist, real_t *rwork, len_t nmsgr,
		                   MPI_Fint mpicomm, MPI_Fint *xcomm, int nolx,
		                   int *tdsx_sor_ptrs, len_t nsor, real_t *tdg,
		                   bool *fact_flags, bool shm_enabled,
		                   void *puper, void *halof)> relax_lines;

		if (dir == relax_dir::x)
			relax_lines = MPI_BMG2_SymStd_relax_lines_x_ml;

		relax_lines(k, sod.data(), bd.data(), x.data(), sord.data(), res.data(),
		            so.len(0), so.len(1), topo.is(0), topo.is(1),
		            kf, nstencil, BMG_RELAX_SYM, updown,
		            datadist, rwork[kf-k].data(), rwork[kf-k].len(0),
		            fcomm, this->comms.data(), this->nlevels,
		            sor_ptrs[kf-k].data(), tdg[kf-k].len(0), tdg[kf-k].data(),
		            fact_flags[kf-k], shm, &puper, halof);
	}


protected:
	relax_dir dir;                          /** Direction of simultaneously relaxed unknowns */
	int nlevels;                            /** Number of levels in tridiagonal solve */
	std::vector<array<int, 1>> sor_ptrs;    /** Pointers into tridiagonal workspace */
	std::vector<bool*> fact_flags;          /** flags for factorization */
	std::vector<array<real_t, 1>> tdg;      /** Tridiagonal solver workspace array */
	std::vector<array<real_t, 1>> rwork;    /** Buffer workspace for tridiagonal solver */
	array<MPI_Fint, 2> comms;               /** Communicators for ml lines */
	MPI_Comm shm_comm;
	MPI_Comm node_comm;
	MPI_Win shm_win;
	double *shm_buff;
	bool shm;
	ml_relax_pup puper;
	bool initialized;
};

}}}

#endif
