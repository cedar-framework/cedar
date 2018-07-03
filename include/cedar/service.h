#ifndef CEDAR_SERVICE_H
#define CEDAR_SERVICE_H

#include <memory>
#include <cedar/kernel_params.h>

namespace cedar {

	class service
	{
	public:
		void add_params(std::shared_ptr<kernel_params> params)
		{
			this->params = params;
		}

	protected:
		std::shared_ptr<kernel_params> params;
	};
}

#endif
