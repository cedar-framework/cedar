#include <cedar/util/timer.h>
#include <cedar/3d/mpi/msg_exchanger.h>


using namespace cedar;
using namespace cedar::cdr3;
using namespace cedar::cdr3::mpi;
using namespace cedar::cdr3::impls;

#include <cedar/types.h>
#include <cedar/2d/ftn/mpi/BMG_workspace_c.h>
#include <cedar/2d/ftn/mpi/BMG_parameters_c.h>


extern "C" {
	using namespace cedar;
	void BMG_get_bc(int, int*);
	void BMG3_SymStd_SETUP_MSG(int *pMSG, int *pMSGSO, len_t *imsg_geom,
	                           len_t nmsgi, int *pSI_MSG, int IBC, len_t *IGRD,
	                           int nog, int nogm, int nproc, int myproc,
	                           len_t *dimx, len_t *dimy, len_t *dimz,
	                           len_t *dimxfine, len_t *dimyfine, len_t *dimzfine,
	                           int *proc_grid, int nproci, int nprocj, int nprock,
	                           int mpicomm);
	void BMG3_SymStd_SETUP_fine_stencil(int kf, real_t *so,
	                                    len_t nlx, len_t nly, len_t nlz,
	                                    int nstencil,
	                                    len_t *iwork, len_t nmsgi, int *pMSGSO,
	                                    real_t *buffer, len_t nmsgr,
	                                    int mpicomm);
	void BMG3_SymStd_UTILS_update_ghosts(int k, real_t *x,
	                                     len_t nx, len_t ny, len_t nz,
	                                     len_t *iwork, len_t nmsgi, int *pMSG,
	                                     real_t *buffer, len_t nmsgr, int nog);
	void BMG3_SymStd_UTILS_update_stencil_ghosts(int k, real_t *so,
	                                             len_t nx, len_t ny, len_t nz, len_t *iwork,
	                                             len_t nmsgi, int *pMSGSO,
	                                             real_t *buffer, len_t nmsgr, int nog);
}

MsgCtx::MsgCtx(grid_topo & topo) :
	pMSG(NBMG_pMSG, topo.nlevel()),
	pLS(NBMG_pLS, topo.nlevel()),
	pMSGSO(NBMG_pMSG, topo.nlevel()),
	proc_grid(topo.nproc(0), topo.nproc(1), topo.nproc(2)),
	proc_coord(topo.nproc(0)*topo.nproc(1)*topo.nproc(2)*3),
	dimxfine(topo.nproc(0)),
	dimyfine(topo.nproc(1)),
	dimzfine(topo.nproc(2)),
	dimx(topo.nproc(0), topo.nlevel()),
	dimy(topo.nproc(1), topo.nlevel()),
	dimz(topo.nproc(2), topo.nlevel())
{
	len_t NLx = topo.nlocal(0) - 2;
	len_t NLy = topo.nlocal(1) - 2;
	len_t NLz = topo.nlocal(2) - 2;
	// !
	// !  Storage for NLx, NLy, NLz arrays for all levels
	// !
	len_t NMSGi = (topo.nproc(0)+topo.nproc(1)+topo.nproc(2))*topo.nlevel();
	// !
	// !  Add storage for MSG workspace
	// !
	NMSGi = NMSGi + 2*(topo.nproc() + 3
	                   + 2*( 6*(NLx+5)*(NLy+5)
	                         + 6*(NLx+5)*(NLz+5)
	                         + 6*(NLy+5)*(NLz+5) )
	                   + topo.nlevel()*( 80 + 20* topo.nproc()));

	NMSGi = NMSGi + 2*(topo.nproc() + 3
	                   + ( 6*(NLx+5)*(NLy+5)
	                       + 6*(NLx+5)*(NLz+5)
	                       + 6*(NLy+5)*(NLz+5) )
	                   + ( 80 + 20* topo.nproc()));
	msg_geom.resize(NMSGi);
	//  !
	//  !  Need to fix the this bound!!
	//  !
	//  !  - MSG part is probably an incorrect extension of 2D
	//  !  - Coarse-grid solve
	//  !

	//  NMSGr = max(
	// &     2*max((NLx+6)*(NLy+6),(NLx+6)*(NLz+6),(NLy+6)*(NLz+6)),
	// &     10*(14*NGx_c*NGy_c*NGz_c+7) )

	p_NLx_kg = 1;
	p_NLy_kg = p_NLx_kg + topo.nproc(0)*topo.nlevel();
	p_NLz_kg = p_NLy_kg + topo.nproc(1)*topo.nlevel();
	pSI_MSG = p_NLz_kg + topo.nproc(2)*topo.nlevel();
	grid_topo cg_topo(topo.igrd_ptr(), 0, topo.nlevel());
	len_t NGx_c = cg_topo.nglobal(0) - 2;
	len_t NGy_c = cg_topo.nglobal(1) - 2;
	len_t NGz_c = cg_topo.nglobal(2) - 2;
	len_t NMSGr = std::max((NLx+6)*(NLy+6), (NLx+6)*(NLz+6));
	NMSGr = std::max(NMSGr, (NLy+6)*(NLz+6));
	NMSGr = std::max(2*NMSGr, 10*(14*NGx_c*NGy_c*NGz_c+7));
	msg_buffer.resize(NMSGr);

	auto pc = proc_coord.begin();
	for (auto k : range(topo.nproc(2))) {
		for (auto j : range(topo.nproc(1))) {
			for (auto i : range(topo.nproc(0))) {
				proc_grid(i,j,k) = k*topo.nproc(1)*topo.nproc(0) + j*topo.nproc(0) + i + 1;
				*pc = i+1;
				pc++;
				*pc = j+1;
				pc++;
				*pc = k+1;
				pc++;
			}
		}
	}

	if (topo.dimyfine.size() > 0) {
		dimyfine = topo.dimyfine;
	} else {
		for (auto j : range(topo.nproc(1))) {
			dimyfine[j] = topo.nlocal(1) - 2;
		}
	}

	if (topo.dimxfine.size() > 0) {
		dimxfine = topo.dimxfine;
	} else {
		for (auto i : range(topo.nproc(0))) {
			dimxfine[i] = topo.nlocal(0) - 2;
		}
	}

	if (topo.dimzfine.size() > 0) {
		dimzfine = topo.dimzfine;
	} else {
		for (auto k : range(topo.nproc(2))) {
			dimzfine[k] = topo.nlocal(2) - 2;
		}
	}
}


namespace cedar { namespace cdr3 { namespace mpi {


void msg_exchanger::setup(std::vector<topo_ptr> topos)
{
	ctx = std::make_unique<MsgCtx>(*topos[0]);
	dims.init(topos[0]->nlevel(), 3);
	coord.init(3);

	auto & topo = *topos[0];
	for (auto i : range<int>(3)) {
		coord(i) = topo.coord(i);
	}

		int rank;
		int ibc;

		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		MPI_Comm_rank(topo.comm, &rank);
		rank++; // Fortran likes to be difficult...

		BMG_get_bc(params->per_mask(), &ibc);

		BMG3_SymStd_SETUP_MSG(ctx->pMSG.data(), ctx->pMSGSO.data(),
		                      ctx->msg_geom.data(), ctx->msg_geom.size(),
		                      &ctx->pSI_MSG, ibc, topo.IGRD(), topo.nlevel(),
		                      topo.nlevel(), topo.nproc(), rank,
		                      ctx->dimx.data(), ctx->dimy.data(), ctx->dimz.data(),
		                      ctx->dimxfine.data(), ctx->dimyfine.data(), ctx->dimzfine.data(),
		                      ctx->proc_grid.data(), topo.nproc(0), topo.nproc(1), topo.nproc(2),
		                      fcomm);
		this->init();
}


void msg_exchanger::init()
{
	for (auto i : range<len_t>(dims.len(0))) {
		dims(i, 0) = ctx->dimx(coord(0), i) + 2;
		dims(i, 1) = ctx->dimy(coord(1), i) + 2;
		dims(i, 2) = ctx->dimz(coord(2), i) + 2;
	}
}



void msg_exchanger::exchange_func(int k, real_t *gf)
{
	timer_begin("halo");
	BMG3_SymStd_UTILS_update_ghosts(k, gf,
	                                dims(k-1,0), dims(k-1,1), dims(k-1,2),
	                                ctx->msg_geom.data(),
	                                ctx->msg_geom.size(), ctx->pMSG.data(), ctx->msg_buffer.data(),
	                                ctx->msg_buffer.size(), dims.len(0));
	timer_end("halo");
}


void msg_exchanger::exchange_sten(int k, real_t * so)
{
	timer_begin("halo-stencil");
	BMG3_SymStd_UTILS_update_stencil_ghosts(k, so, dims(k-1, 0), dims(k-1, 1), dims(k-1, 2),
	                                        ctx->msg_geom.data(), ctx->msg_geom.size(),
	                                        ctx->pMSGSO.data(), ctx->msg_buffer.data(),
	                                        ctx->msg_buffer.size(), dims.len(0));
	timer_end("halo-stencil");
}


void msg_exchanger::run(mpi::grid_func & f)
{
		grid_topo &topo = f.grid();

		BMG3_SymStd_UTILS_update_ghosts(topo.level()+1, f.data(),
		                                f.len(0), f.len(1), f.len(2),
		                                ctx->msg_geom.data(), ctx->msg_geom.size(),
		                                ctx->pMSG.data(), ctx->msg_buffer.data(),
		                                ctx->msg_buffer.size(), topo.nlevel());
}


}}}
