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
end module test_sofa_m

program test_sofa
  use test_sofa_m
  implicit none

  write (*,*) "use crate::common::*;";
  call test_pom
  call test_nutation
end program test_sofa
