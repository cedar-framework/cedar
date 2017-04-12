#include <cedar/kernel_params.h>

using namespace cedar;

namespace cedar
{
std::shared_ptr<kernel_params> build_kernel_params(config::reader & conf)
{

	auto params = std::make_shared<kernel_params>();

	for (int i = 0; i < 3; ++i) params->periodic[i] = false;
	params->relax_symmetric = true;
	params->definite = true;

	return params;
}
}
