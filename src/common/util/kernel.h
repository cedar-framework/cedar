#ifndef BOXMG_KERNEL_H
#define BOXMG_KERNEL_H

#include <functional>

namespace boxmg {

class KernelBase
{
public:
	virtual ~KernelBase() {};
};


template <typename... Args>
class Kernel : public KernelBase
{
private:
	using kern_t = std::function<void(Args...)>;
	kern_t kern;

public:
    Kernel(kern_t&& kern): kern(std::move(kern)) {}
	void operator()(Args... args) { if (kern) kern(args...); }
};

}

#endif
