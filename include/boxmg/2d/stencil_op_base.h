#ifndef BOXMG_STENCIL_OP_BASE
#define BOXMG_STENCIL_OP_BASE

#include <boxmg/types.h>
#include <boxmg/discrete_op.h>
#include <boxmg/kernel_name.h>
#include <boxmg/2d/grid_stencil.h>
#include <typeinfo>

namespace boxmg { namespace bmg2d {

template <class grid_func, class registry>
class stencil_op_base : public discrete_op<grid_func>
{

public:
	stencil_op_base() {}
    stencil_op_base(len_t nx, len_t ny, bool intergrid=false) : gs(nx,ny,1,intergrid) {}
	real_t *data() { return gs.data(); }
	grid_stencil & stencil() { return gs; }
	const grid_stencil & stencil() const { return gs; }
	virtual void set_registry(std::shared_ptr<registry> kreg) { kernels = kreg; }
	virtual std::shared_ptr<registry> get_registry() const { return kernels; }
	template <class stencil_op>
	void apply(const grid_func &x, grid_func &y) const
	{
		if (kernels != nullptr) {
			kernels->run(kernel_name::matvec,
			             static_cast<const stencil_op&>(*this),
			             x,y);
		} else {
			log::error << "Kernel registry for stencil op not set!" << std::endl;
		}
	}


	template <class stencil_op>
	void residual(const grid_func &x, const grid_func &b, grid_func &r) const
	{
		if (kernels != nullptr) {
			kernels->run(kernel_name::residual,
			             static_cast<const stencil_op&>(*this),
			             x,b,r);
		} else {
			log::error << "Kernel registry for stencil op not set!" << std::endl;
		}
	}


	template <class stencil_op>
	grid_func residual(const grid_func &x, const grid_func &b) const
	{
		auto r = grid_func(x.shape(0),x.shape(1));
		this->residual<stencil_op>(x,b,r);

		return r;
	}

protected:
	grid_stencil gs;
	std::shared_ptr<registry> kernels;

};

}}

#endif
