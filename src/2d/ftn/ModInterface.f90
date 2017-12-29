module ModInterface

       use iso_c_binding, only: C_CHAR, C_NULL_CHAR, C_INT, C_DOUBLE, C_PTR, C_BOOL
       interface
         subroutine print_error(string) bind(C, name="print_error")
           use iso_c_binding, only: c_char
           character(kind=c_char) :: string(*)
         end subroutine print_error
       end interface
       !call print_c(C_CHAR_"Hello World"//C_NULL_CHAR)

       interface
          subroutine ftimer_begin(string) bind(C, name="ftimer_begin")
            use iso_c_binding, only : c_char
            character(kind=c_char) :: string(*)
          end subroutine ftimer_begin
       end interface

       interface
          subroutine ftimer_end(string) bind(C, name="ftimer_end")
            use iso_c_binding, only : c_char
            character(kind=c_char) :: string(*)
          end subroutine ftimer_end
       end interface

       integer, parameter :: real_t = C_DOUBLE
       integer, parameter :: len_t = C_INT

       interface
          subroutine halo_exchange(k, Q, halof) bind(C,name="halo_exchange")
            use iso_c_binding, only : C_PTR, C_INT, C_DOUBLE
            integer(C_INT),value :: k
            real(C_DOUBLE) :: Q(*)
            type(C_PTR) :: halof
          end subroutine halo_exchange
       end interface


       interface
          subroutine halo_stencil_exchange(k, SO, halof) bind(C,name="halo_stencil_exchange")
            use iso_c_binding, only : C_PTR, C_INT, C_DOUBLE
            integer(C_INT),value :: k
            real(C_DOUBLE) :: SO(*)
            type(C_PTR) :: halof
          end subroutine halo_stencil_exchange
       end interface


end module ModInterface
