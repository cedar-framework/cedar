#include "boxmg/util/kernel_registry.h"


using namespace boxmg;

void KernelRegistry::set(const std::string & kname, const std::string & kid)
{
	active[kname] = avail[kname].at(kid);
}
