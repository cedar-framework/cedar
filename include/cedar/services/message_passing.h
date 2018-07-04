#ifndef CEDAR_MESSAGEPASSING_H
#define CEDAR_MESSAGEPASSING_H

#include <mpi.h>

#include <cedar/service.h>

namespace cedar { namespace services {

class message_passing : public service
{
public:
	const static std::string name() { return "message passing"; }
	virtual int comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm) = 0;
};

}}

#endif
