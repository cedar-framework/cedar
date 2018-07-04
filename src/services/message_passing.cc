#include <cedar/services/message_passing.h>

extern "C" {
	void cedar_comm_split(cedar::services::message_passing *mp,
	                      MPI_Fint fcomm, int color, int key, MPI_Fint *fnewcomm)
	{
		MPI_Comm comm = MPI_Comm_f2c(fcomm);
		MPI_Comm newcomm;
		mp->comm_split(comm, color, key, &newcomm);
		*fnewcomm = MPI_Comm_c2f(newcomm);
	}
}

