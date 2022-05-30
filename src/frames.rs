//fn celestial_to_true()

// Very precise version
//
// c2i00a
//   pnm00a
//     pn00a
//       nut00a
//         fal03
//         faf03
//         faom03
//         fame03
//         fave03
//         fae03
//         fama03
//         faju03
//         fasa03
//         faur03
//         fapa03
//       pn00
//         pr00
//         obl80
//         numat
//   c2ibpn
//     bpn2xy
//     c2ixy
//       c2ixys

// *  2) The matrix RC2I is the first stage in the transformation from
// *     celestial to terrestrial coordinates:
// *
// *        [TRS]  =  RPOM * R_3(ERA) * RC2I * [CRS]
// *
// *               =  RC2T * [CRS]
// *
// *     where [CRS] is a vector in the Geocentric Celestial Reference
// *     System and [TRS] is a vector in the International Terrestrial
// *     Reference System (see IERS Conventions 2003), ERA is the Earth
// *     Rotation Angle and RPOM is the polar motion matrix.

// c2i00b
//   pnm00b
//     pn00b
//       nut00b
//       pn00
//         pr00
//         obl80
//         bp00
//           bi00
//           pr00
//           cr
//   c2ibpn
//     bpn2xy
//     c2ixy
//       c2ixys
use crate::common::*;
use crate::time::{TT,DJ00,DJC};

/// Form the celestial-to-intermediate matrix for a given date using the
/// IAU 2000B precession-nutation model.
///
/// Source: c2i00b.for
pub fn celestial_to_intermediate(tt:TT)->Mat3 {
    let rbpn = precession_nutation(tt); // pnm00b
    celestial_to_intermediate_with_bpn(tt,&rbpn) //c2ibpn
}

/// Form the matrix of precession-nutation for a given date (including
/// frame bias), equinox-based, IAU 2000B model.
///
/// Source: pnm00b.for
pub fn precession_nutation(tt:TT)->Mat3 {
    let pn : PrecessionNutation = tt.into();
    pn.rbpn
}

#[derive(Copy,Clone,Debug)]
pub struct PrecessionNutation {
    /// Nutation
    pub dpsi:R,

    /// Mean obliquity
    epsa:R,

    /// Frame bias matrix
    rb:Mat3,

    /// Precession matrix
    rp:Mat3,

    /// Bias-precession matrix
    rbp:Mat3,

    /// Nutation matrix
    rn:Mat3,

    /// GCRS-to-true matrix
    rbpn:Mat3
}

//TT((date1,date2))

impl From<TT> for PrecessionNutation {
    /// Source: pn00b.for
    fn from(tt:TT)->Self {
	let (dpsi,deps) = nutation(tt); // nut00b
	Self::from_tt_with_nutation(tt,dpsi,deps) // pn00
    }
}

impl PrecessionNutation {
    /// Source: pn00.for
    pub fn from_tt_with_nutation(tt:TT,dpsi:R,deps:R)->Self {
	let (_dpsipr,depspr) = precession_rate(tt); // pr00
	let epsa = mean_obliquity(tt) + depspr; // obl80
	let (rb,rp,rbp) = bias_and_precession(tt); // bp00
	let rn = nutation_matrix(epsa,dpsi,deps); // numat
	let rbpn = rn.compose(&rbp);
	Self {
	    dpsi,
	    epsa,
	    rb,
	    rp,
	    rbp,
	    rn,
	    rbpn
	}
    }
}

const DPPLAN : R = - 0.135 * DMAS2R;
const DEPLAN : R =   0.388 * DMAS2R;

const NALS : &[[i8;5]] = &[
    [ 0,    0,    0,    0,    1],
    [ 0,    0,    2,   -2,    2],
    [ 0,    0,    2,    0,    2],
    [ 0,    0,    0,    0,    2],
    [ 0,    1,    0,    0,    0],
    [ 0,    1,    2,   -2,    2],
    [ 1,    0,    0,    0,    0],
    [ 0,    0,    2,    0,    1],
    [ 1,    0,    2,    0,    2],
    [ 0,   -1,    2,   -2,    2],
    [ 0,    0,    2,   -2,    1],
    [-1,    0,    2,    0,    2],
    [-1,    0,    0,    2,    0],
    [ 1,    0,    0,    0,    1],
    [-1,    0,    0,    0,    1],
    [-1,    0,    2,    2,    2],
    [ 1,    0,    2,    0,    1],
    [-2,    0,    2,    0,    1],
    [ 0,    0,    0,    2,    0],
    [ 0,    0,    2,    2,    2],
    [ 0,   -2,    2,   -2,    2],
    [-2,    0,    0,    2,    0],
    [ 2,    0,    2,    0,    2],
    [ 1,    0,    2,   -2,    2],
    [-1,    0,    2,    0,    1],
    [ 2,    0,    0,    0,    0],
    [ 0,    0,    2,    0,    0],
    [ 0,    1,    0,    0,    1],
    [-1,    0,    0,    2,    1],
    [ 0,    2,    2,   -2,    2],
    [ 0,    0,   -2,    2,    0],
    [ 1,    0,    0,   -2,    1],
    [ 0,   -1,    0,    0,    1],
    [-1,    0,    2,    2,    1],
    [ 0,    2,    0,    0,    0],
    [ 1,    0,    2,    2,    2],
    [-2,    0,    2,    0,    0],
    [ 0,    1,    2,    0,    2],
    [ 0,    0,    2,    2,    1],
    [ 0,   -1,    2,    0,    2],
    [ 0,    0,    0,    2,    1],
    [ 1,    0,    2,   -2,    1],
    [ 2,    0,    2,   -2,    2],
    [-2,    0,    0,    2,    1],
    [ 2,    0,    2,    0,    1],
    [ 0,   -1,    2,   -2,    1],
    [ 0,    0,    0,   -2,    1],
    [-1,   -1,    0,    2,    0],
    [ 2,    0,    0,   -2,    1],
    [ 1,    0,    0,    2,    0],
    [ 0,    1,    2,   -2,    1],
    [ 1,   -1,    0,    0,    0],
    [-2,    0,    2,    0,    2],
    [ 3,    0,    2,    0,    2],
    [ 0,   -1,    0,    2,    0],
    [ 1,   -1,    2,    0,    2],
    [ 0,    0,    0,    1,    0],
    [-1,   -1,    2,    2,    2],
    [-1,    0,    2,    0,    0],
    [ 0,   -1,    2,    2,    2],
    [-2,    0,    0,    0,    1],
    [ 1,    1,    2,    0,    2],
    [ 2,    0,    0,    0,    1],
    [-1,    1,    0,    1,    0],
    [ 1,    1,    0,    0,    0],
    [ 1,    0,    2,    0,    0],
    [-1,    0,    2,   -2,    1],
    [ 1,    0,    0,    0,    2],
    [-1,    0,    0,    1,    0],
    [ 0,    0,    2,    1,    2],
    [-1,    0,    2,    4,    2],
    [-1,    1,    0,    1,    1],
    [ 0,   -2,    2,   -2,    1],
    [ 1,    0,    2,    2,    1],
    [-2,    0,    2,    2,    2],
    [-1,    0,    0,    0,    2],
    [ 1,    1,    2,   -2,    2],
];

const CLS : &[[R;6]] = &[
    [-172064161.0,-174666.0,33386.0,92052331.0,9086.0,15377.0],
    [-13170906.0,-1675.0,-13696.0,5730336.0,-3015.0,-4587.0],
    [-2276413.0,-234.0,2796.0,978459.0,-485.0,1374.0],
    [2074554.0,207.0,-698.0,-897492.0,470.0,-291.0],
    [1475877.0,-3633.0,11817.0,73871.0,-184.0,-1924.0],
    [-516821.0,1226.0,-524.0,224386.0,-677.0,-174.0],
    [711159.0,73.0,-872.0,-6750.0,0.0,358.0],
    [-387298.0,-367.0,380.0,200728.0,18.0,318.0],
    [-301461.0,-36.0,816.0,129025.0,-63.0,367.0],
    [215829.0,-494.0,111.0,-95929.0,299.0,132.0],
    [128227.0,137.0,181.0,-68982.0,-9.0,39.0],
    [123457.0,11.0,19.0,-53311.0,32.0,-4.0],
    [156994.0,10.0,-168.0,-1235.0,0.0,82.0],
    [63110.0,63.0,27.0,-33228.0,0.0,-9.0],
    [-57976.0,-63.0,-189.0,31429.0,0.0,-75.0],
    [-59641.0,-11.0,149.0,25543.0,-11.0,66.0],
    [-51613.0,-42.0,129.0,26366.0,0.0,78.0],
    [45893.0,50.0,31.0,-24236.0,-10.0,20.0],
    [63384.0,11.0,-150.0,-1220.0,0.0,29.0],
    [-38571.0,-1.0,158.0,16452.0,-11.0,68.0],
    [32481.0,0.0,0.0,-13870.0,0.0,0.0],
    [-47722.0,0.0,-18.0,477.0,0.0,-25.0],
    [-31046.0,-1.0,131.0,13238.0,-11.0,59.0],
    [28593.0,0.0,-1.0,-12338.0,10.0,-3.0],
    [20441.0,21.0,10.0,-10758.0,0.0,-3.0],
    [29243.0,0.0,-74.0,-609.0,0.0,13.0],
    [25887.0,0.0,-66.0,-550.0,0.0,11.0],
    [-14053.0,-25.0,79.0,8551.0,-2.0,-45.0],
    [15164.0,10.0,11.0,-8001.0,0.0,-1.0],
    [-15794.0,72.0,-16.0,6850.0,-42.0,-5.0],
    [21783.0,0.0,13.0,-167.0,0.0,13.0],
    [-12873.0,-10.0,-37.0,6953.0,0.0,-14.0],
    [-12654.0,11.0,63.0,6415.0,0.0,26.0],
    [-10204.0,0.0,25.0,5222.0,0.0,15.0],
    [16707.0,-85.0,-10.0,168.0,-1.0,10.0],
    [-7691.0,0.0,44.0,3268.0,0.0,19.0],
    [-11024.0,0.0,-14.0,104.0,0.0,2.0],
    [7566.0,-21.0,-11.0,-3250.0,0.0,-5.0],
    [-6637.0,-11.0,25.0,3353.0,0.0,14.0],
    [-7141.0,21.0,8.0,3070.0,0.0,4.0],
    [-6302.0,-11.0,2.0,3272.0,0.0,4.0],
    [5800.0,10.0,2.0,-3045.0,0.0,-1.0],
    [6443.0,0.0,-7.0,-2768.0,0.0,-4.0],
    [-5774.0,-11.0,-15.0,3041.0,0.0,-5.0],
    [-5350.0,0.0,21.0,2695.0,0.0,12.0],
    [-4752.0,-11.0,-3.0,2719.0,0.0,-3.0],
    [-4940.0,-11.0,-21.0,2720.0,0.0,-9.0],
    [7350.0,0.0,-8.0,-51.0,0.0,4.0],
    [4065.0,0.0,6.0,-2206.0,0.0,1.0],
    [6579.0,0.0,-24.0,-199.0,0.0,2.0],
    [3579.0,0.0,5.0,-1900.0,0.0,1.0],
    [4725.0,0.0,-6.0,-41.0,0.0,3.0],
    [-3075.0,0.0,-2.0,1313.0,0.0,-1.0],
    [-2904.0,0.0,15.0,1233.0,0.0,7.0],
    [4348.0,0.0,-10.0,-81.0,0.0,2.0],
    [-2878.0,0.0,8.0,1232.0,0.0,4.0],
    [-4230.0,0.0,5.0,-20.0,0.0,-2.0],
    [-2819.0,0.0,7.0,1207.0,0.0,3.0],
    [-4056.0,0.0,5.0,40.0,0.0,-2.0],
    [-2647.0,0.0,11.0,1129.0,0.0,5.0],
    [-2294.0,0.0,-10.0,1266.0,0.0,-4.0],
    [2481.0,0.0,-7.0,-1062.0,0.0,-3.0],
    [2179.0,0.0,-2.0,-1129.0,0.0,-2.0],
    [3276.0,0.0,1.0,-9.0,0.0,0.0],
    [-3389.0,0.0,5.0,35.0,0.0,-2.0],
    [3339.0,0.0,-13.0,-107.0,0.0,1.0],
    [-1987.0,0.0,-6.0,1073.0,0.0,-2.0],
    [-1981.0,0.0,0.0,854.0,0.0,0.0],
    [4026.0,0.0,-353.0,-553.0,0.0,-139.0],
    [1660.0,0.0,-5.0,-710.0,0.0,-2.0],
    [-1521.0,0.0,9.0,647.0,0.0,4.0],
    [1314.0,0.0,0.0,-700.0,0.0,0.0],
    [-1283.0,0.0,0.0,672.0,0.0,0.0],
    [-1331.0,0.0,8.0,663.0,0.0,4.0],
    [1383.0,0.0,-2.0,-594.0,0.0,-2.0],
    [1405.0,0.0,4.0,-610.0,0.0,2.0],
    [1290.0,0.0,0.0,-556.0,0.0,0.0],
];

const NLS : usize = 77;
const U2R : R = DAS2R/1e7;

/// Source: nut00b.for
pub fn nutation(TT((date1,date2)):TT)->(R,R) {
    let t = (( date1 - DJ00 ) + date2) / DJC;
    let el = ((485868.249036 + 1717915923.2178 * t) % TURNAS) * DAS2R;
    let elp = ((1287104.79305 + 129596581.0481 * t) % TURNAS) * DAS2R;
    let f   = ((335779.526232 + 1739527262.8478 * t) % TURNAS) * DAS2R;
    let d   = ((1072260.70369 + 1602961601.2090 * t) % TURNAS) * DAS2R;
    let om  = ((450160.398036 - 6962890.5431 * t) % TURNAS) * DAS2R;

    let mut dp = 0.0;
    let mut de = 0.0;

// *  Summation of luni-solar nutation series (in reverse order).
    for i in (0..NLS).rev() {
	let arg = (NALS[i][0] as R * el  +
		   NALS[i][1] as R * elp +
		   NALS[i][2] as R * f   +
		   NALS[i][3] as R * d   +
		   NALS[i][4] as R * om) % TWO_PI;
	let sarg = sin(arg);
	let carg = cos(arg);

	dp += ( CLS[i][0] + CLS[i][1] * t ) * sarg + CLS[i][2] * carg;
	de += ( CLS[i][3] + CLS[i][4] * t ) * carg + CLS[i][5] * sarg;
    }

    let dpsils = dp * U2R;
    let depsls = de * U2R;

    let dpsipl = DPPLAN;
    let depspl = DEPLAN;

    let dpsi = dpsils + dpsipl;
    let deps = depsls + depspl;

    (dpsi,deps)
}

/// Source: pr00.for
pub fn precession_rate(TT((date1,date2)):TT)->(R,R) {
    const PRECOR : R = -0.29965 * DAS2R;
    const OBLCOR : R = -0.02524 * DAS2R;
    let t = ( ( date1 - DJ00 ) + date2 ) / DJC;
    let dpsipr = PRECOR * t;
    let depspr = OBLCOR * t;
    (dpsipr,depspr)
}

/// Source: obl80.for
pub fn mean_obliquity(TT((date1,date2)):TT)->R {
    let t = ( ( date1 - DJ00 ) + date2 ) / DJC;
    DAS2R * ( 84381.448 +
              ( -46.8150 +
		 ( -0.00059 +
                    0.001813 * t ) * t ) * t )
}

/// Returns (RB,RP,RBP)
/// Source: bp00.for
pub fn bias_and_precession(tt@TT((date1,date2)):TT)->(Mat3,Mat3,Mat3) {
    let t = ( ( date1 - DJ00 ) + date2 ) / DJC;
    let (dpsibi,depsbi,dra0) = frame_bias();

    const EPS0 : R = 84381.448 * DAS2R;

    let psia77 = ( 5038.7784 + ( -1.07259 + ( -0.001147 ) * t ) * t ) * t * DAS2R;
    let oma77 = EPS0 + (( 0.05127 + ( -0.007726 ) * t ) * t ) * t * DAS2R;
    let chia = ( 10.5526 + ( -2.38064 + ( -0.001125 ) * t ) * t ) * t * DAS2R;

    let (dpsipr,depspr) = precession_rate(tt); // pr00

    let psia = psia77 + dpsipr;
    let oma = oma77 + depspr;

    let rb = 
	Mat3::identity()
	.compose(&Mat3::rotation(0,-depsbi))
	.compose(&Mat3::rotation(1,-dpsibi*sin(EPS0)))
	.compose(&Mat3::rotation(2,dra0));

    let rp =
	Mat3::identity()
	.compose(&Mat3::rotation(2,chia))
	.compose(&Mat3::rotation(0,-oma))
	.compose(&Mat3::rotation(2,-psia))
	.compose(&Mat3::rotation(0,EPS0));

    let rbp = rp.compose(&rb);

    (rb,rp,rbp)
}

/// Source: bi00.for
pub fn frame_bias()->(R,R,R) {
    const DPBIAS : R = -0.041775 * DAS2R;
    const DEBIAS : R = -0.0068192 * DAS2R;
    const DRA0 : R = -0.0146 * DAS2R;
    (DPBIAS,DEBIAS,DRA0)
}

/// Source: numat.for
pub fn nutation_matrix(epsa:R,dpsi:R,deps:R)->Mat3 {
    Mat3::identity()
	.compose(&Mat3::rotation(0,-(epsa + deps)))
	.compose(&Mat3::rotation(2,-dpsi))
	.compose(&Mat3::rotation(0,epsa))
}

/// Source: c2ibpn.for
pub fn celestial_to_intermediate_with_bpn(tt:TT,rbpn:&Mat3)->Mat3 {
    let (x,y) = xy_from_bpn(rbpn); // bpn2xy
    celestial_to_intermediate_with_bpn_and_xy(tt,x,y)
}

/// Source: c2ixy.for
pub fn celestial_to_intermediate_with_bpn_and_xy(tt:TT,x:R,y:R)->Mat3 {
    let s00 = cio_locator(tt,x,y); // s00
    celestial_to_intermediate_from_xys(x,y,s00) // c2ixys
}

/// Source: c2ixys.for
pub fn celestial_to_intermediate_from_xys(x:R,y:R,s:R)->Mat3 {
    let r2 = x*x + y*y;
    let e =
	if r2 > 0.0 {
            atan2(y, x)
	} else {
	    0.0
	};
    let d = atan(sqrt(r2/(1.0 - r2)));
    Mat3::identity()
	.compose(&Mat3::rotation(2,-(e + s)))
	.compose(&Mat3::rotation(1,d))
	.compose(&Mat3::rotation(2,e))
}

/// Source: bpn2xy.for
pub fn xy_from_bpn(rbpn:&Mat3)->(R,R) {
    (rbpn[2][0],rbpn[2][1])
}

/// XXX: Approximation
/// The CIO locator s is the difference between the right ascensions
/// of the same point in two systems: the two systems are the GCRS and
/// the CIP,CIO, and the point is the ascending node of the CIP
/// equator.  The quantity s remains below 0.1 arcsecond throughout
/// 1900-2100.
///
/// Source: s00.for
pub fn cio_locator(_tt:TT,_x:R,_y:R)->R {
    0.0
}

// /// Source: c2ibn.for
// pub fn celestial_to_intemediate_with_bpn(tt:TT,rbpn:Mat3)->Mat3 {
    
// }
