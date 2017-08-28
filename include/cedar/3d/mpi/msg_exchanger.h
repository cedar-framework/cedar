#ifndef CEDAR_3D_MSG_EXCHANGER
#define CEDAR_3D_MSG_EXCHANGER

#include "cedar/2d/ftn/mpi/BMG_workspace_c.h"
#include <cedar/array.h>
#include <cedar/mpi/grid_topo.h>
#include <cedar/halo_exchanger.h>
#include <cedar/3d/mpi/grid_func.h>

namespace cedar { namespace cdr3 { namespace kernel {
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


	class msg_exchanger : public halo_exchanger
	{
	public:
		msg_exchanger(grid_topo & topo);
		void init();
		MsgCtx & context() { return ctx; }
		void * context_ptr() { return &ctx;}
		virtual void exchange_func(int k, real_t *gf);
		virtual void exchange_sten(int k, real_t * so);

	private:
		MsgCtx ctx;
		/**
		   Local array extents by grid number needed for MSG ghost
		   exchange calls.  The first dimension is the grid number,
		   the second is the dimension.
		**/
		array<len_t, 2> dims;
		array<len_t, 1> coord;
	};


}}}}
#endif
