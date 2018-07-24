#include <cedar/services/message_passing.h>

extern "C" {
	void cedar_comm_split(cedar::services::message_passing *mp,
	                      MPI_Fint fcomm, int color, int key, MPI_Fint *fnewcomm,
	                      int *ierr)
	{
		MPI_Comm comm = MPI_Comm_f2c(fcomm);
		MPI_Comm newcomm;
		*ierr = mp->comm_split(comm, color, key, &newcomm);
		*fnewcomm = MPI_Comm_c2f(newcomm);
	}


	void cedar_gather_impl(cedar::services::message_passing *mp,
	                       const void *sendbuf, int sendcount, MPI_Fint fsendtype,
	                       void *recvbuf, int recvcount, MPI_Fint frecvtype,
	                       int root, MPI_Fint fcomm, int *ierr)
	{
		MPI_Comm comm = MPI_Comm_f2c(fcomm);
		MPI_Datatype sendtype = MPI_Type_f2c(fsendtype);
		MPI_Datatype recvtype = MPI_Type_f2c(frecvtype);

		*ierr = mp->gather(sendbuf, sendcount, sendtype,
		                   recvbuf, recvcount, recvtype,
		                   root, comm);
	}


	void cedar_scatter_impl(cedar::services::message_passing *mp,
	                        const void *sendbuf, int sendcount, MPI_Fint fsendtype,
	                        void *recvbuf, int recvcount, MPI_Fint frecvtype,
	                        int root, MPI_Fint fcomm, int *ierr)
	{
		MPI_Comm comm = MPI_Comm_f2c(fcomm);
		MPI_Datatype sendtype = MPI_Type_f2c(fsendtype);
		MPI_Datatype recvtype = MPI_Type_f2c(frecvtype);

		*ierr = mp->scatter(sendbuf, sendcount, sendtype,
		                    recvbuf, recvcount, recvtype,
		                    root, comm);
	}
}

