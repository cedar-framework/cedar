#ifndef CEDAR_INTERFACE_TOPO_H
#define CEDAR_INTERFACE_TOPO_H

#include <memory>

#include <cedar/mpi/grid_topo.h>
#include <cedar/interface/object.h>


std::shared_ptr<cedar::grid_topo> cedar_topo_getobj(cedar_topo handle);


#endif
