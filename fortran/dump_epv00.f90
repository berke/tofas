module epv00
  implicit none
  contains
    subroutine dump_array(name,x)
      character(len=*), intent(in) :: name
      real(8), intent(in) :: x(:,:,:)
      integer :: m,n,o
      integer :: i,j,k

      m = size(x,1)
      n = size(x,2)
      o = size(x,3)

      ! We want :
      !    Rust    x[i][j][k]
      !
      !    equiv. to
      !
      !    Fortran x(i,j,k)
      !
      ! Therefore if size(x) = [m,n,o]
      ! the type in Rust will be
      !
      !    x : [[[f64;o];n];m]

      write (*,'("pub const ",A," : [[[f64;",I0,"];",I0,"];",I0,"] = [")') name,o,n,m
      do i=1,m
         write (*,'("  [")')
         do j=1,n
            write (*,'("    [")',advance='no')
            do k=1,o
               write (*,'(E23.15)',advance='no') x(i,j,k)
               if (k < o) then
                  write (*,'(",")',advance='no')
               end if
            end do
            write (*,'("], // ",I0,",",I0,",",I0,":",I0)') i,j,1,o
         end do
         write (*,'("  ],")')
      end do
      write (*,'("];")')
    end subroutine dump_array
  end module epv00
program dump_epv00
  call iau_EPV00
end program
