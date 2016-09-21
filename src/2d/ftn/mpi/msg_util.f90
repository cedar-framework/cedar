      subroutine MSG_pause(comm) BIND(C,NAME='MSG_pause')
        implicit none
        include 'mpif.h'
        include 'geom_param_fort_f90.h'
        include 'mpi_param_fort_f90.h'
        include 'MSG_f90.h'

        integer comm

        comm = MSG_COMM

      end subroutine MSG_pause


      subroutine MSG_play(comm) BIND(C,NAME='MSG_play')
        implicit none
        include 'mpif.h'
        include 'geom_param_fort_f90.h'
        include 'mpi_param_fort_f90.h'
        include 'MSG_f90.h'

        integer,value :: comm

        MSG_COMM = comm

      end subroutine MSG_play
