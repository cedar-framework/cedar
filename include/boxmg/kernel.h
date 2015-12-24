#ifndef BOXMG_KERNEL_H
#define BOXMG_KERNEL_H

#include <functional>

namespace boxmg {

class kernel_base
{
public:
	virtual ~kernel_base() {};
};


template <typename... Args>
class kernel : public kernel_base
{
private:
	using kern_t = std::function<void(Args...)>;
	kern_t kern;

public:
    kernel(kern_t&& kern): kern(std::move(kern)) {}
	void operator()(Args... args) { if (kern) kern(args...); }
};

}

#endif
