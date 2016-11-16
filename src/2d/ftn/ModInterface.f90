module ModInterface

       use iso_c_binding, only: C_CHAR, C_NULL_CHAR, C_INT, C_DOUBLE
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

end module ModInterface
