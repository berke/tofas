module test_sofa_m
  implicit none
  double precision, parameter :: DJ00 = 2451545.0, DJC = 36525.0
contains
  subroutine test_pom
    double precision :: xp,yp,sp,rpom(3,3)
    integer, parameter :: nstep = 10
    integer :: istep
    write (*,*) "pub const POM00_DATA : &[(R,R,R,[[R;3];3])] = &[";
    do istep=1,nstep
       call random_number(xp)
       call random_number(yp)
       call random_number(sp)
       call iau_POM00(xp,yp,sp,rpom)
       write (*,*) "(",xp,",",yp,",",sp,",", &
            "[", &
            "[",rpom(1,1),",",rpom(1,2),",",rpom(1,3),"],", &
            "[",rpom(2,1),",",rpom(2,2),",",rpom(2,3),"],", &
            "[",rpom(3,1),",",rpom(3,2),",",rpom(3,3),"],", &
            "]),"
    end do
    write (*,*) "];";
  end subroutine test_pom

  subroutine test_nutation
    double precision :: date1,date2,dpsi,deps
    integer, parameter :: nstep = 10
    integer :: istep
    write (*,*) "pub const NUT00B_DATA : &[(R,R,R,R)] = &[";
    do istep=1,nstep
       date1 = DJ00
       call random_number(date2)
       date2 = date2*DJC
       call iau_NUT00B(date1,date2,dpsi,deps)
       write (*,*) "(",date1,",",date2,",",dpsi,",",deps,"),"
    end do
    write (*,*) "];";
  end subroutine test_nutation

  subroutine write_mat(a)
    double precision, intent(in) :: a(:,:)
    integer :: i,j
    write (*,*) "["
    do i=1,size(a,1)
       write (*,"('[')",advance='no')
       do j=1,size(a,2)
          if (j > 1) write (*,"(',')",advance='no')
          write (*,"(E15.8)",advance='no') a(i,j)
       end do
       write (*,"('],')",advance='no')
    end do
    write (*,*) "]"
  end subroutine write_mat

  subroutine test_bias_and_precession
    double precision :: date1,date2
    double precision :: rb(3,3),rp(3,3),rbp(3,3)
    integer, parameter :: nstep = 10
    integer :: istep
    write (*,*) "pub const BP00_DATA : &[(R,R,[[[R;3];3];3])] = &[";
    do istep=1,nstep
       date1 = DJ00
       call random_number(date2)
       date2 = date2*DJC
       call iau_BP00(date1,date2,rb,rp,rbp)
       write (*,*) "(",date1,",",date2,",["
       call write_mat(rb)
       write (*,*) ","
       call write_mat(rp)
       write (*,*) ","
       call write_mat(rbp)
       write (*,*) "]),"
    end do
    write (*,*) "];";
  end subroutine test_bias_and_precession
end module test_sofa_m

program test_sofa
  use test_sofa_m
  implicit none

  write (*,*) "use crate::common::*;";
  call test_pom
  call test_nutation
  call test_bias_and_precession
end program test_sofa
