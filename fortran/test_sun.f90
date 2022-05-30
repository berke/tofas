program test_sun
  implicit none
  double precision, parameter :: &
       as2r = 4.848136811095359935889141d-6, &
       D2PI = 6.283185307179586476925287D0
  INTEGER :: IY, IM, ID, IH, MIN, I, J
  DOUBLE PRECISION :: SEC, XP, YP, DUT1, DDP80, DDE80, DX00, DY00, &
       DX06, DY06, DJMJD0, DATE, TIME, UTC, DAT, TAI, TT, TUT, UT1, &
       RP(3,3), DP80, DE80, DPSI, DEPS, EPSA, RN(3,3), RNPB(3,3), EE, GST, &
       RC2TI(3,3), RPOM(3,3), RC2IT(3,3), X, Y, S, RC2I(3,3), ERA, DP00, &
       DE00, RB(3,3), RPB(3,3), V1(3), V2(3), DDP00, DDE00, &
       RC2T(3,3),PVH(3,2),PVB(3,2)
  DOUBLE PRECISION iau_OBL80, iau_EQEQ94, iau_ANP, iau_GMST82, iau_ERA00, iau_SP00, iau_EE00, iau_GMST00

  IY = 2007
  IM = 4
  ID = 5
  IH = 12
  MIN = 0
  SEC = 0D0

  XP = 0.0349282D0 * AS2R
  YP = 0.4833163D0 * AS2R

  DUT1 = -0.072073685D0
  ! Nutation corrections wrt IAU 1976/1980 (mas->radians).
  DDP80 = -55.0655D0 * AS2R/1000D0
  DDE80 = -6.3580D0 * AS2R/1000D0
  ! CIP offsets wrt IAU 2000A (mas->radians).
  DX00 = 0.1725D0 * AS2R/1000D0
  DY00 = -0.2650D0 * AS2R/1000D0
  ! CIP offsets wrt IAU 2006/2000A (mas->radians).
  DX06 = 0.1750D0 * AS2R/1000D0
  DY06 = -0.2259D0 * AS2R/1000D0
  ! TT (MJD).
  CALL iau_CAL2JD ( IY, IM, ID, DJMJD0, DATE, J )
  TIME = ( 60D0*(60D0*DBLE(IH) + DBLE(MIN)) + SEC ) / 86400D0
  UTC = DATE + TIME
  CALL iau_DAT ( IY, IM, ID, TIME, DAT, J )
  TAI = UTC + DAT/86400D0
  TT = TAI + 32.184D0/86400D0
  ! UT1.
  TUT = TIME + DUT1/86400D0
  UT1 = DATE + TUT

  call iau_C2T00B( DJMJD0,TT, DJMJD0,UT1, XP,YP, RC2T)

  print *,'TT=',TT
  print *,'TUT=',TUT
  print *,'UT1=',UT1
  print *,'RC2T='
  do i=1,3
     do j=1,3
        print *,RC2T(i,j)
     end do
  end do


  call iau_EPV00( DJMJD0,TT, PVH,PVB, J )
  print *,'PH=',PVH(:,1)

end program test_sun
