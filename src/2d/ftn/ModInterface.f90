module ModInterface

       use iso_c_binding, only: C_CHAR, C_NULL_CHAR, C_INT, C_DOUBLE
       interface
         subroutine print_error(string) bind(C, name="print_error")
           use iso_c_binding, only: c_char
           character(kind=c_char) :: string(*)
         end subroutine print_error
       end interface
       !call print_c(C_CHAR_"Hello World"//C_NULL_CHAR)

       integer, parameter :: real_t = C_DOUBLE
       integer, parameter :: len_t = C_INT

end module ModInterface
