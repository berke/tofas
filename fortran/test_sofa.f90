program test_sofa
  implicit none
  double precision :: xp,yp,sp,rpom(3,3)
  integer, parameter :: nstep = 10
  integer :: istep

  write (*,*) "use crate::common::*;";
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
end program test_sofa
