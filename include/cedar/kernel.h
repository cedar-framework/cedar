#ifndef CEDAR_KERNEL_H
#define CEDAR_KERNEL_H

#include <memory>
#include <functional>

#include <cedar/config/reader.h>
#include <cedar/kernel_params.h>

namespace cedar {

class kernel_base
{
public:
	virtual ~kernel_base() {};
};


template <typename... Args>
class kernel : public kernel_base
{
private:
	using kern_t = std::function<void(const kernel_params &, Args...)>;
	kern_t kern;
	std::shared_ptr<kernel_params> params;

public:
    kernel(kern_t&& kern, std::shared_ptr<kernel_params> params): kern(std::move(kern)), params(params) {}
	void operator()(Args... args) { if (kern) kern(*params, args...); }
};

}

#endif
