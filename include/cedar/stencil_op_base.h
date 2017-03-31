#ifndef CEDAR_STENCIL_OP_BASE
#define CEDAR_STENCIL_OP_BASE

#include <cedar/types.h>
#include <cedar/discrete_op.h>
#include <cedar/kernel_name.h>
#include <typeinfo>

namespace cedar {

template <class grid_func, class registry, class grid_stencil>
class stencil_op_base : public discrete_op<grid_func>
{

public:
	stencil_op_base() {}
	virtual real_t *data() = 0;
	virtual grid_stencil & stencil() = 0;
	virtual const grid_stencil & stencil() const = 0;
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
		auto r = grid_func::like(x);
		this->residual<stencil_op>(x,b,r);

		return r;
	}

protected:
	std::shared_ptr<registry> kernels;

};

}

#endif
