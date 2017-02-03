#include <boxmg/2d/kernel/mpi/factory.h>
#include <boxmg/2d/mpi/gallery.h>

namespace boxmg { namespace bmg2d { namespace mpi { namespace gallery {

using namespace boxmg;

stencil_op poisson(topo_ptr grid)
{
	config::reader conf("");
	auto kreg = kernel::mpi::factory::from_config(conf);

	auto so = stencil_op(grid);
	so.set_registry(kreg);
	auto & sten = so.stencil();
	sten.five_pt() = true;
	auto & topo = so.grid();

	sten.set(0);

	real_t nlx = topo.nlocal(0) - 2;
	real_t nly = topo.nlocal(1) - 2;
	real_t ngx = topo.nglobal(0) - 2;
	real_t ngy = topo.nglobal(1) - 2;

	real_t igs = topo.is(0);
	real_t jgs = topo.is(1);

	real_t hx = 1.0 / (topo.nglobal(0) - 1);
	real_t hy = 1.0 / (topo.nglobal(1) - 1);

	real_t nlx_g = nlx + 2;
	real_t nly_g = nly + 2;

	real_t i1 = nlx + 1;
	real_t j1 = nly + 1;
	real_t i2 = nlx;
	real_t j2 = nly;

	real_t igf = igs + i2 - 1;
	real_t jgf = jgs + j2 - 1;

	real_t xh = hy / hx;
	real_t yh = hx / hy;

	real_t ibeg = 1;
	real_t jbeg = 1;

	if (igs == 1)
		ibeg++;
	if (jgs == 1)
		jbeg++;

	real_t iend = nlx_g;
	real_t jend = nly_g;

	if (igf == ngx)
		iend--;
	if (jgf == ngy)
		jend--;

	for (auto j : range<len_t>(jbeg, jend)) {
		for (auto i : range<len_t>(1, iend)) {
			sten(i,j,dir::S) = 1.0 * yh;
		}
	}

	for (auto j : range<len_t>(1, jend)) {
		for (auto i : range<len_t>(ibeg, iend)) {
			sten(i,j,dir::W) = xh;
		}
	}

	for (auto j : range<len_t>(1, j1)) {
		for (auto i : range<len_t>(1, i1)) {
			sten(i,j,dir::C) = 2*xh + 2*yh;
		}
	}

	return so;
}



stencil_op diag_diffusion(topo_ptr grid, real_t dx, real_t dy)
{
	config::reader conf("");
	auto kreg = kernel::mpi::factory::from_config(conf);

	auto so = stencil_op(grid);
	so.set_registry(kreg);
	auto & sten = so.stencil();
	sten.five_pt() = true;
	auto & topo = so.grid();

	sten.set(0);

	real_t nlx = topo.nlocal(0) - 2;
	real_t nly = topo.nlocal(1) - 2;
	real_t ngx = topo.nglobal(0) - 2;
	real_t ngy = topo.nglobal(1) - 2;

	real_t igs = topo.is(0);
	real_t jgs = topo.is(1);

	real_t hx = 1.0 / (topo.nglobal(0) - 1);
	real_t hy = 1.0 / (topo.nglobal(1) - 1);

	real_t nlx_g = nlx + 2;
	real_t nly_g = nly + 2;

	real_t i1 = nlx + 1;
	real_t j1 = nly + 1;
	real_t i2 = nlx;
	real_t j2 = nly;

	real_t igf = igs + i2 - 1;
	real_t jgf = jgs + j2 - 1;

	real_t xh = hy / hx;
	real_t yh = hx / hy;

	real_t ibeg = 1;
	real_t jbeg = 1;

	if (igs == 1)
		ibeg++;
	if (jgs == 1)
		jbeg++;

	real_t iend = nlx_g;
	real_t jend = nly_g;

	if (igf == ngx)
		iend--;
	if (jgf == ngy)
		jend--;

	for (auto j : range<len_t>(jbeg, jend)) {
		for (auto i : range<len_t>(1, iend)) {
			sten(i,j,dir::S) = dy * yh;
		}
	}

	for (auto j : range<len_t>(1, jend)) {
		for (auto i : range<len_t>(ibeg, iend)) {
			sten(i,j,dir::W) = dx * xh;
		}
	}

	for (auto j : range<len_t>(1, j1)) {
		for (auto i : range<len_t>(1, i1)) {
			sten(i,j,dir::C) = 2*dx*xh + 2*dy*yh;
		}
	}

	return so;
}

}}}}
