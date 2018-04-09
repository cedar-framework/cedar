#include <cedar/util/timer.h>
#include <cedar/types.h>
#include <cedar/2d/ftn/mpi/BMG_workspace_c.h>
#include <cedar/2d/ftn/mpi/BMG_parameters_c.h>
#include <cedar/2d/mpi/msg_exchanger.h>

extern "C" {
	using namespace cedar;
	void BMG2_SymStd_SETUP_MSG(int *pMSG, int *pMSGSO, len_t *iMSG_Geom,
	                           len_t NMSGi, int *pSI_MSG, int IBC, len_t *IGRD,
	                           int nog, int nogm, int nproc, int myproc,
	                           len_t* dimx, len_t *dimy, len_t *dimxfine, len_t *dimyfine,
	                           int *procgrid, int nproci, int nprocj, int mpicomm);
	void BMG2_SymStd_SETUP_LS(len_t *iWorkMSG, len_t NMSGi, int *pMSG, int *pLS, int *pSI_MSG,
	                          int *procgrid, int nproci, int nprocj, int nog);
	void BMG2_SymStd_SETUP_fine_stencil(int kf, real_t *SO, len_t IIF, len_t JJF, int NStncl,
	                                    len_t *iWork, len_t NMSGi, int *pMSGSO,
	                                    real_t *msg_buffer, len_t NMSGr, int mpicomm);
	void BMG2_SymStd_UTILS_update_ghosts(int K, real_t *x, len_t Nx, len_t Ny, len_t *iWork,
	                                     len_t NMSGi, int *pMSG,
	                                     real_t *buffer, len_t NMSGr, int nog, int mpicomm);
	void BMG2_SymStd_UTILS_update_stencil_ghosts(int K, real_t *x, len_t Nx, len_t Ny, len_t *iWork,
	                                             len_t NMSGi, int *pMSGSO,
	                                             real_t *buffer, len_t NMSGr, int nog, int mpicomm);
	void BMG_get_bc(int, int*);
}


using namespace cedar;
using namespace cedar::cdr2;
using namespace cedar::cdr2::impls;
using namespace cedar::cdr2::mpi;

MsgCtx::MsgCtx(grid_topo & topo) :
	pMSG(NBMG_pMSG, topo.nlevel()),
	pLS(NBMG_pLS, topo.nlevel()),
	pMSGSO(NBMG_pMSG, topo.nlevel()),
	proc_grid(topo.nproc(0), topo.nproc(1)),
	proc_coord(topo.nproc(0)*topo.nproc(1)*2),
	dimxfine(topo.nproc(0)),
	dimyfine(topo.nproc(1)),
	dimx(topo.nproc(0), topo.nlevel()),
	dimy(topo.nproc(1), topo.nlevel()),
	comm(topo.comm)
{
	MPI_Comm_split(topo.comm, topo.coord(1), topo.coord(0),
	               &xlinecomm);
	MPI_Comm_split(topo.comm, topo.coord(0), topo.coord(1),
	               &ylinecomm);
	len_t NLx = topo.nlocal(0) - 2;
	len_t NLy = topo.nlocal(1) - 2;
	// !
	// !  Storage for NLx, NLy, NLz arrays for all levels
	// !
	len_t NMSGi = (topo.nproc(0)+topo.nproc(1))*topo.nlevel();
	// !
	// !  Add storage for MSG workspace (Not a sharp bound)
	// !
	NMSGi = NMSGi + topo.nproc()+2 + topo.nlevel()*(16*(NLx+NLy+6)+18*topo.nproc()+26);
	// !
	// ! Add storage for MSGSO workspace
	// !
	NMSGi = NMSGi + topo.nproc()+2 + topo.nlevel()*(24*(NLx+NLy+4)+18*topo.nproc()+26);
	// !
	// ! Add storage for Line Solves
	// !
	NMSGi = NMSGi + 4*topo.nlevel()*topo.nproc();
	msg_geom.resize(NMSGi);

	// !
	// !  Need to fix the this bound!!
	// !
	// !  - MSG real buffer space 
	// !  - Workspace for coarse-grid solve communication.
	// !
	// NMSGr = std::max(2*std::max(topo.nproc(0)+6, topo.nproc(1)+6),
	//                  10*(5*NGx_c*NGy_c+5),
	//               &              4*NLy*NProcI + 5*NLx, 4*NLx*NProcJ + 5*NLy )

	p_NLx_kg = 1;
	p_NLy_kg = p_NLx_kg + topo.nproc(0)*topo.nlevel();
	pSI_MSG = p_NLy_kg + topo.nproc(1)*topo.nlevel();
	grid_topo cg_topo(topo.igrd_ptr(), 0, topo.nlevel());
	len_t NGx_c = cg_topo.nglobal(0) - 2;
	len_t NGy_c = cg_topo.nglobal(1) - 2;
	len_t NMSGr = std::max(4*NLy*topo.nproc(0) + 5*NLx,
	                       4*NLx*topo.nproc(1) + 5*NLy);
	NMSGr = std::max(NMSGr, 10*(5*NGx_c*NGy_c+5));
	NMSGr = std::max(NMSGr, 2*std::max(NLx+6, NLy+6));
	msg_buffer.resize(NMSGr);

	auto pc = proc_coord.begin();
	for (auto j : range(topo.nproc(1))) {
		dimyfine[j] = topo.nlocal(1) - 2;
		for (auto i : range(topo.nproc(0))) {
			proc_grid(i,j) = j*topo.nproc(0) + i + 1;
			*pc = i+1;
			pc++;
			*pc = j+1;
			pc++;
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
}

namespace cedar { namespace cdr2 { namespace mpi {

void msg_exchanger::setup(std::vector<topo_ptr> topos)
{
	ctx = std::make_unique<MsgCtx>(*topos[0]);
	dims.init(topos[0]->nlevel(), 2);
	coord.init(2);

	auto & topo = *topos[0];
	coord(0) = topo.coord(0);
	coord(1) = topo.coord(1);

	int rank;
	int ibc;
	MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

	MPI_Comm_rank(topo.comm, &rank);
	rank++; // Fortran likes to be difficult...

	BMG_get_bc(params->per_mask(), &ibc);

	BMG2_SymStd_SETUP_MSG(ctx->pMSG.data(), ctx->pMSGSO.data(),
	                      ctx->msg_geom.data(), ctx->msg_geom.size(),
	                      &ctx->pSI_MSG, ibc, topo.IGRD(), topo.nlevel(),
	                      topo.nlevel(), topo.nproc(), rank,
	                      ctx->dimx.data(), ctx->dimy.data(),
	                      ctx->dimxfine.data(), ctx->dimyfine.data(),
	                      ctx->proc_grid.data(), topo.nproc(0), topo.nproc(1),
	                      fcomm);

	BMG2_SymStd_SETUP_LS(ctx->msg_geom.data(), ctx->msg_geom.size(),
	                     ctx->pMSG.data(), ctx->pLS.data(), &ctx->pSI_MSG,
	                     ctx->proc_grid.data(), topo.nproc(0), topo.nproc(1),
	                     topo.nlevel());


	for (auto i : range<len_t>(dims.len(0))) {
		dims(i, 0) = ctx->dimx(coord(0), i) + 2;
		dims(i, 1) = ctx->dimy(coord(1), i) + 2;
	}
}


void msg_exchanger::run(mpi::stencil_op<five_pt> & sop)
{
	grid_topo &topo = sop.grid();
	int nstencil;

	nstencil = 3;

	MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

	timer_begin("halo-stencil");
	BMG2_SymStd_SETUP_fine_stencil(topo.level()+1, sop.data(),
	                               sop.len(0), sop.len(1), nstencil,
	                               ctx->msg_geom.data(), ctx->msg_geom.size(),
	                               ctx->pMSGSO.data(), ctx->msg_buffer.data(),
	                               ctx->msg_buffer.size(), fcomm);
	timer_end("halo-stencil");
}


void msg_exchanger::run(mpi::stencil_op<nine_pt> & sop)
{
	grid_topo &topo = sop.grid();
	int nstencil;

	nstencil = 5;

	MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

	timer_begin("halo-stencil");
	BMG2_SymStd_SETUP_fine_stencil(topo.level()+1, sop.data(),
	                               sop.len(0), sop.len(1), nstencil,
	                               ctx->msg_geom.data(), ctx->msg_geom.size(),
	                               ctx->pMSGSO.data(), ctx->msg_buffer.data(),
	                               ctx->msg_buffer.size(), fcomm);
	timer_end("halo-stencil");
}


void msg_exchanger::run(mpi::grid_func & f)
{
	grid_topo &topo = f.grid();

	MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

	timer_begin("halo");
	BMG2_SymStd_UTILS_update_ghosts(topo.level()+1, f.data(), f.len(0), f.len(1), ctx->msg_geom.data(),
	                                ctx->msg_geom.size(), ctx->pMSG.data(), ctx->msg_buffer.data(),
	                                ctx->msg_buffer.size(), topo.nlevel(), fcomm);
	timer_end("halo");
}


void msg_exchanger::exchange_func(int k, real_t *gf)
{
		MPI_Fint fcomm = MPI_Comm_c2f(this->ctx->comm);

		timer_begin("halo");
		BMG2_SymStd_UTILS_update_ghosts(k, gf, dims(k-1, 0), dims(k-1, 1),
		                                ctx->msg_geom.data(),
		                                ctx->msg_geom.size(), ctx->pMSG.data(), ctx->msg_buffer.data(),
		                                ctx->msg_buffer.size(), dims.len(0), fcomm);
		timer_end("halo");
}


void msg_exchanger::exchange_sten(int k, real_t * so)
{
	MPI_Fint fcomm = MPI_Comm_c2f(this->ctx->comm);

	timer_begin("halo-stencil");
	BMG2_SymStd_UTILS_update_stencil_ghosts(k, so, dims(k-1, 0), dims(k-1, 1),
	                                        ctx->msg_geom.data(), ctx->msg_geom.size(),
	                                        ctx->pMSGSO.data(), ctx->msg_buffer.data(),
	                                        ctx->msg_buffer.size(), dims.len(0), fcomm);
	timer_end("halo-stencil");
}

}}}
