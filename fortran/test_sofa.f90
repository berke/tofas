module test_sofa_m
  implicit none
  double precision, parameter :: DJ00 = 2451545.0, DJC = 36525.0, TWO_PI = 6.283185307179586476925287
  integer, parameter :: NSTEP = 1000
contains
  subroutine test_pom
    double precision :: xp,yp,sp,rpom(3,3)
    integer :: istep
    write (*,*) "pub const POM00_DATA : &[(R,R,R,[[R;3];3])] = &["
    do istep=1,NSTEP
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
    write (*,*) "];"
  end subroutine test_pom

  subroutine test_nutation_matrix
    double precision :: epsa,dpsi,deps,rmatn(3,3)
    integer :: istep
    write (*,*) "pub const NUMAT_DATA : &[(R,R,R,[[R;3];3])] = &["
    do istep=1,NSTEP
       call random_number(epsa)
       epsa = epsa*TWO_PI
       call random_number(dpsi)
       dpsi = dpsi*TWO_PI
       call random_number(deps)
       deps = deps*TWO_PI
       call iau_NUMAT(epsa,dpsi,deps,rmatn)
       write (*,*) "(",epsa,",",dpsi,",",deps,","
       call write_mat(rmatn)
       write (*,*) "),"
    end do
    write (*,*) "];"
  end subroutine test_nutation_matrix

  subroutine test_nutation
    double precision :: date1,date2,dpsi,deps
    integer :: istep
    write (*,*) "pub const NUT00B_DATA : &[(R,R,R,R)] = &["
    do istep=1,NSTEP
       date1 = DJ00
       call random_number(date2)
       date2 = date2*DJC
       call iau_NUT00B(date1,date2,dpsi,deps)
       write (*,*) "(",date1,",",date2,",",dpsi,",",deps,"),"
    end do
    write (*,*) "];"
  end subroutine test_nutation

  subroutine write_mat(a)
    double precision, intent(in) :: a(:,:)
    integer :: i,j
    write (*,*) "["
    do i=1,size(a,1)
       write (*,"('[')",advance='no')
       do j=1,size(a,2)
          if (j > 1) write (*,"(',')",advance='no')
          write (*,"(E25.17)",advance='no') a(i,j)
       end do
       write (*,"('],')",advance='no')
    end do
    write (*,*) "]"
  end subroutine write_mat

  subroutine test_celestial_to_terrestrial_from_cio_components
    double precision :: rc2i(3,3),era,rpom(3,3),rc2t(3,3)
    integer :: istep
    write (*,*) "pub const C2TCIO_DATA : &[([[R;3];3],R,[[R;3];3],[[R;3];3])] = &["
    do istep=1,NSTEP
       call random_number(rc2i)
       call random_number(era)
       call random_number(rpom)
       call iau_C2TCIO(rc2i,era,rpom,rc2t)
       write (*,*) "("
       call write_mat(rc2i)
       write (*,*) ",",era,","
       call write_mat(rpom)
       write (*,*) ","
       call write_mat(rc2t)
       write (*,*) "),"
    end do
    write (*,*) "];"
  end subroutine test_celestial_to_terrestrial_from_cio_components

  subroutine test_bias_and_precession
    double precision :: date1,date2
    double precision :: rb(3,3),rp(3,3),rbp(3,3)
    integer :: istep
    write (*,*) "pub const BP00_DATA : &[(R,R,[[[R;3];3];3])] = &["
    do istep=1,NSTEP
       call random_number(date1)
       date1 = date1 + DJ00
       call random_number(date2)
       call iau_BP00(date1,date2,rb,rp,rbp)
       write (*,*) "(",date1,",",date2,",["
       call write_mat(rb)
       write (*,*) ","
       call write_mat(rp)
       write (*,*) ","
       call write_mat(rbp)
       write (*,*) "]),"
    end do
    write (*,*) "];"
  end subroutine test_bias_and_precession

  subroutine test_rotation_angle
    double precision iau_ERA00
    double precision :: ut11,ut12,era
    integer :: istep
    write (*,*) "pub const ERA00_DATA : &[(R,R,R)] = &["
    do istep=1,NSTEP
       ut11 = DJ00
       call random_number(ut12)
       ut12 = ut12*DJC
       era = iau_ERA00(ut11,ut12)
       write (*,*) "(",ut11,",",ut12,",",era,"),"
    end do
    write (*,*) "];"
  end subroutine test_rotation_angle

  subroutine test_precession_rate
    double precision :: tt1,tt2,dpsipr,depspr
    integer :: istep
    write (*,*) "pub const PR00_DATA : &[(R,R,R,R)] = &["
    do istep=1,NSTEP
       call random_number(tt1)
       tt1 = DJ00 + DJC*tt1
       call random_number(tt2)
       call iau_PR00(tt1,tt2,dpsipr,depspr)
       write (*,*) "(",tt1,",",tt2,",",dpsipr,",",depspr,"),"
    end do
    write (*,*) "];"
  end subroutine test_precession_rate

  subroutine test_precession_nutation
    double precision :: tt1,tt2,rbpn(3,3)
    integer :: istep

    write (*,*) "pub const PNM00B_DATA : &[(R,R,[[R;3];3])] = &["
    do istep=1,NSTEP
       call random_number(tt1)
       tt1 = DJ00 + DJC*tt1
       call random_number(tt2)
       call iau_PNM00B(tt1,tt2,rbpn)
       write (*,*) "(",tt1,",",tt2,","
       call write_mat(rbpn)
       write (*,*) "),"
    end do
    write (*,*) "];"
  end subroutine test_precession_nutation

  subroutine test_celestial_to_intermediate
    double precision :: tt1,tt2,rc2i(3,3)
    integer :: istep

    write (*,*) "pub const C2I00B_DATA : &[(R,R,[[R;3];3])] = &["
    do istep=1,NSTEP
       call random_number(tt1)
       tt1 = DJ00 + DJC*tt1
       call random_number(tt2)
       call iau_C2I00B(tt1,tt2,rc2i)
       write (*,*) "(",tt1,",",tt2,","
       call write_mat(rc2i)
       write (*,*) "),"
    end do
    write (*,*) "];"
  end subroutine test_celestial_to_intermediate

  subroutine test_locator
    double precision :: tt1,tt2,x,y,s00
    double precision iau_S00
    integer :: istep
    
    write (*,*) "pub const S00_DATA : &[(R,R,R,R,R)] = &["
    do istep=1,NSTEP
       call random_number(tt1)
       tt1 = DJ00 + DJC*tt1
       call random_number(tt2)
       call random_number(x)
       call random_number(y)
       s00 = iau_S00(tt1,tt2,x,y)
       write (*,*) "(",tt1,",",tt2,",",x,",",y,",",s00,"),"
    end do
    write (*,*) "];"
  end subroutine test_locator

  subroutine test_celestial_to_terrestrial
    double precision :: tt1,tt2,tdb1,tdb2,ut11,ut12,xp,yp,era,rc2t(3,3)
    double precision iau_ERA00
    double precision, parameter :: DTR = 0.0,DT = 0.0
    integer :: err,istep

    xp = 0
    yp = 0
    write (*,*) "pub const C2T00B_DATA : &[(R,R,R,R,R,R,R,R,[[R;3];3])] = &["
    do istep=1,NSTEP
       call random_number(tt1)
       tt1 = DJ00 + DJC*tt1
       call random_number(tt2)
       call iau_TTTDB(tt1,tt2,dtr,tdb1,tdb2,err)
       call iau_TTUT1(tt1,tt2,dt,ut11,ut12,err)
       era = iau_ERA00(ut11,ut12)
       call iau_C2T00B(tt1,tt2,ut11,ut12,xp,yp,rc2t)

       write (*,*) "(",tt1,",",tt2,",",dtr,",",tdb1,",",tdb2,",",dt,",",ut11,",",ut12,","
       call write_mat(rc2t)
       write (*,*) "),"
    end do
    write (*,*) "];"
  end subroutine test_celestial_to_terrestrial
end module test_sofa_m

program test_sofa
  use test_sofa_m
  implicit none

  write (*,*) "use crate::common::*;"
  call test_pom
  call test_nutation
  call test_nutation_matrix
  call test_bias_and_precession
  call test_rotation_angle
  call test_precession_rate
  call test_precession_nutation
  call test_celestial_to_intermediate
  call test_celestial_to_terrestrial_from_cio_components
  call test_locator
  call test_celestial_to_terrestrial
end program test_sofa
