#ifndef CEDAR_3D_MPI_PLANE_UTIL_H
#define CEDAR_3D_MPI_PLANE_UTIL_H

#include <tuple>
#include <type_traits>

#include <cedar/types.h>
#include <cedar/3d/mpi/types.h>
#include <cedar/2d/mpi/solver.h>
#include <cedar/3d/mpi/kernel_manager.h>
#include <cedar/3d/mpi/plane_ult.h>
#include <cedar/3d/mpi/plane_team.h>
#include <cedar/3d/mpi/plane_util.h>


namespace cedar { namespace cdr3 { namespace mpi {

template<relax_dir rdir>
struct ortho_dir : std::integral_constant<unsigned short, 0> {};
template<>
struct ortho_dir<relax_dir::xy> : std::integral_constant<unsigned short, 4> {};
template<>
struct ortho_dir<relax_dir::xz> : std::integral_constant<unsigned short, 2> {};
template<>
struct ortho_dir<relax_dir::yz> : std::integral_constant<unsigned short, 1> {};


std::tuple<int,MPI_Comm> log_begin(bool log_planes, int ipl, const std::string & suff, MPI_Comm comm);
void log_end(bool log_planes, std::tuple<int,MPI_Comm> saved);

cdr2::mpi::kman_ptr master_kman(config & conf, int nplanes, bool aggregate, plane_team & team);
cdr2::mpi::kman_ptr worker_kman(config & conf, int nplanes, bool aggregate, plane_team & team, int worker_id);

template<class fsten>
void setup_agg_solve(cdr2::mpi::solver<fsten> & slv)
{
	service_manager<cdr2::mpi::stypes> & sman = slv.get_kernels()->services();
	sman.set<services::halo_exchange<cdr2::mpi::stypes>>("plane");
	sman.set<services::message_passing>("plane");
	slv.setup_halo();
}


template<class fsten>
using slv2_ptr = std::unique_ptr<cdr2::mpi::solver<fsten>>;

template<relax_dir rdir, class sten3>
void copy_rhs(const stencil_op<sten3> & so,
              const grid_func & x,
              const grid_func & b,
              cdr2::mpi::grid_func & b2,
              int ipl);

template<relax_dir rdir>
void copy32(const grid_func & x, cdr2::mpi::grid_func & x2, int ipl);

template<relax_dir rdir>
void copy23(const cdr2::mpi::grid_func & x2, grid_func &x, int ipl);


template<relax_dir rdir, class sten3, class sten2>
void relax_planes(const stencil_op<sten3> & so, grid_func & x,
                  const grid_func & b, cycle::Dir cdir,
                  services::halo_exchange<stypes> & halo_service,
                  std::vector<slv2_ptr<sten2>> & planes);


template<relax_dir rdir, class sten3, class sten2>
void relax_planes_agg(const stencil_op<sten3> & so, grid_func & x,
                      const grid_func & b, cycle::Dir cdir,
                      services::halo_exchange<stypes> & halo_service,
                      std::vector<slv2_ptr<sten2>> & planes,
                      std::array<plane_ult<sten2>, 2> & threads);


template<relax_dir rdir>
struct plane_util
{
	static void copy_coeff(const stencil_op<seven_pt> & so3,
	                cdr2::mpi::stencil_op<cdr2::five_pt> & so2, int ipl);
	static void copy_coeff(const stencil_op<xxvii_pt> & so3,
	                       cdr2::mpi::stencil_op<cdr2::nine_pt> & so2, int ipl);
};


}}}

#endif
