#ifndef CEDAR_3D_MSG_EXCHANGER
#define CEDAR_3D_MSG_EXCHANGER

#include <cedar/2d/ftn/mpi/BMG_workspace_c.h>
#include <cedar/mpi/grid_topo.h>
#include <cedar/3d/mpi/types.h>
#include <cedar/services/halo_exchange.h>

extern "C" {
	using namespace cedar;
	void BMG3_SymStd_SETUP_fine_stencil(int kf, real_t *so,
	                                    len_t nlx, len_t nly, len_t nlz,
	                                    int nstencil,
	                                    len_t *iwork, len_t nmsgi, int *pMSGSO,
	                                    real_t *buffer, len_t nmsgr,
	                                    int mpicomm);
}

namespace cedar { namespace cdr3 {
namespace impls
{
	namespace mpi = cedar::cdr3::mpi;
	struct MsgCtx
	{
		MsgCtx(grid_topo & topo);
		aarray<int,int,2> pMSG; // these should be int, len_t
		aarray<int,int,2> pLS;
		aarray<int,int,2> pMSGSO;
		std::vector<len_t> msg_geom;
		aarray<int,int,3> proc_grid;
		std::vector<int> proc_coord;
		std::vector<len_t> dimxfine;
		std::vector<len_t> dimyfine;
		std::vector<len_t> dimzfine;
		aarray<int,len_t,2> dimx;
		aarray<int,len_t,2> dimy;
		aarray<int,len_t,2> dimz;
		std::vector<real_t> msg_buffer;
		int pSI_MSG;
		int p_NLx_kg, p_NLy_kg, p_NLz_kg;
		/* std::vector<len_t> iworkmsg; */
		/* int *iworkmsg[nmsgi]; */
		/* int nmsgi; */
		/* int pmsg[nbmg_pmsg,nog]; */
		/* int msgbuffer[nmsgr]; */
		/* int nmsgr; */
		len_t nlocal(int kg, int dim, int ijkrank) {
			len_t local_arr_ptr = pMSG(ipL_MSG_LocalArraySize,kg) - 1;  // 1 vs 0 based indexing
			len_t idx = local_arr_ptr;
			ijkrank--; // 1 vs 0 based indexing
			idx += ijkrank*3 + dim;

			return msg_geom[idx];
		}

		len_t cg_nlocal(int dim, int ijkrank) {
			return nlocal(0, dim, ijkrank);
		}
	};
}

namespace mpi {

	class msg_exchanger : public services::halo_exchange<stypes>
	{
		using MsgCtx = impls::MsgCtx;
	public:
		void init();
		MsgCtx & context() { return *ctx; }
		void * context_ptr() { return ctx.get();}

		void setup(std::vector<topo_ptr> topos) override;
		void run(stencil_op<seven_pt> & so) override { this->run_impl(so); }
		void run(stencil_op<xxvii_pt> & so) override { this->run_impl(so); }
		void run(grid_func & gf) override;
		void exchange_func(int k, real_t *gf) override;
		void exchange_sten(int k, real_t * so) override;
		aarray<int, len_t, 2> & leveldims(int k) override{
			if (k == 0)
				return ctx->dimx;
			else if (k == 1)
				return ctx->dimy;
			else
				return ctx->dimz;
		}
		len_t * datadist(int k, int grid) override {
			// placeholder
			return ctx->msg_geom.data();
		}
		std::vector<real_t> & linebuf() override {
			return ctx->msg_buffer;
		}

		template<class sten>
		void run_impl(stencil_op<sten> & sop)
		{
			grid_topo &topo = sop.grid();
			int nstencil = stencil_ndirs<sten>::value;

			MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

			BMG3_SymStd_SETUP_fine_stencil(topo.level()+1, sop.data(),
			                               sop.len(0), sop.len(1), sop.len(2),
			                               nstencil,
			                               ctx->msg_geom.data(), ctx->msg_geom.size(),
			                               ctx->pMSGSO.data(), ctx->msg_buffer.data(),
			                               ctx->msg_buffer.size(), fcomm);
		}

	private:
		std::unique_ptr<MsgCtx> ctx;
		/**
		   Local array extents by grid number needed for MSG ghost
		   exchange calls.  The first dimension is the grid number,
		   the second is the dimension.
		**/
		array<len_t, 2> dims;
		array<len_t, 1> coord;
	};
}}}

#endif
