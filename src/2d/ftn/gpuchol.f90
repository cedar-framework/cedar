module gpuchol

  use iso_c_binding, only : C_CHAR, C_NULL_CHAR, C_INT, C_DOUBLE
  interface
     subroutine dpotrs_gpu(uplo, N, NRHS, A, lda, B, ldb, info)
       use iso_c_binding, only : C_CHAR, C_NULL_CHAR, C_INT, C_DOUBLE
       character(kind=c_char) :: uplo
       integer(c_int), value :: N, NRHS, lda, ldb
       integer(c_int) :: info
       real(c_double) :: A(*), B(*)
     end subroutine dpotrs_gpu
  end interface

  interface
     subroutine dpotrf_gpu(uplo, N, A, lda, info)
       use iso_c_binding, only : C_CHAR, C_NULL_CHAR, C_INT, C_DOUBLE
       character(kind=c_char) :: uplo
       integer(c_int), value :: N, lda
       integer(c_int) :: info
       real(c_double) :: A(*)
     end subroutine dpotrf_gpu
  end interface

  interface
     subroutine init_lapack_gpu
     end subroutine init_lapack_gpu
  end interface

  interface
     subroutine finalize_lapack_gpu
     end subroutine finalize_lapack_gpu
  end interface

end module gpuchol
