use crate::common::*;
use crate::time::{TDB,UT1,DJY,DJ00};
use crate::epv00_data::*;

/// Astronomical unit (m, IAU 2012)
pub const AUM : R = 149597870.7e3;

const AM12 : R =  0.000000211284;
const AM13 : R = -0.000000091603;
const AM21 : R = -0.000000230286;
const AM22 : R =  0.917482137087;
const AM23 : R = -0.397776982902;
const AM32 : R =  0.397776982902;
const AM33 : R =  0.917482137087;

const NE0X : usize = 501;
const NE0Y : usize = 501;
const NE0Z : usize = 137;
// const ME0 : usize = NE0X;
const NE1X : usize = 79;
const NE1Y : usize = 80;
const NE1Z : usize = 12;
// const ME1 : usize = NE1Y;
const NE2X : usize = 5;
const NE2Y : usize = 5;
const NE2Z : usize = 3;
// const ME2 : usize = NE2X;
const NS0X : usize = 212;
const NS0Y : usize = 213;
const NS0Z : usize = 69;
// const MS0 : usize = NS0Y;
const NS1X : usize = 50;
const NS1Y : usize = 50;
const NS1Z : usize = 14;
// const MS1 : usize = NS1X;
const NS2X : usize = 9;
const NS2Y : usize = 9;
const NS2Z : usize = 2;
// const MS2 : usize = NS2X;

const NE0 : [usize;3] = [NE0X, NE0Y, NE0Z];
const NE1 : [usize;3] = [NE1X, NE1Y, NE1Z];
const NE2 : [usize;3] = [NE2X, NE2Y, NE2Z];
const NS0 : [usize;3] = [NS0X, NS0Y, NS0Z];
const NS1 : [usize;3] = [NS1X, NS1Y, NS1Z];
const NS2 : [usize;3] = [NS2X, NS2Y, NS2Z];

#[derive(Debug,Clone)]
pub struct PosVel {
    pub p:Vec3,
    pub v:Vec3
}

#[derive(Debug,Clone)]
pub struct EarthPosVel {
    pub heliocentric:PosVel,
    pub barycentric:PosVel,
    pub warning:Option<EarthPosVelWarning>
}

#[derive(Debug,Clone)]
pub enum EarthPosVelWarning {
    DateOutOfRange
}

impl From<TDB> for EarthPosVel {
    /// Earth position and velocity, heliocentric and barycentric, with
    /// respect to the Barycentric Celestial Reference System.
    ///
    /// The date should be in the 1900-2100 AD range
    ///
    /// Source: IAU SOFA epv00.for
    fn from(TDB((date1,date2)):TDB)->Self {
	let t = ( ( date1 - DJ00 ) + date2 ) / DJY;
	let t2 = t*t;
	let warning =
	    if abs(t) < 100.0 {
		None
	    } else {
		Some(EarthPosVelWarning::DateOutOfRange)
	    };

	let mut ph = [0.0;3];
	let mut vh = [0.0;3];
	let mut pb = [0.0;3];
	let mut vb = [0.0;3];
	
	for k in 0..3 {
	    let mut xyz = 0.0;
	    let mut xyzd = 0.0;

	    // Sun to Earth, T^0 terms.
	    for j in 0..NE0[k] {
		let a = E0[0][j][k];
		let b = E0[1][j][k];
		let c = E0[2][j][k];
		let p = b + c*t;
		xyz += a*cos(p);
		xyzd -= a*c*sin(p);
	    }

	    // Sun to Earth, T^1 terms.
	    for j in 0..NE1[k] {
		let a = E1[0][j][k];
		let b = E1[1][j][k];
		let c = E1[2][j][k];
		let ct = c*t;
		let p = b + ct;
		let cp = cos(p);
		xyz += a*t*cp;
		xyzd += a*(cp - ct*sin(p));
	    }

	    // Sun to Earth, T^2 terms.
	    for j in 0..NE2[k] {
		let a = E2[0][j][k];
		let b = E2[1][j][k];
		let c = E2[2][j][k];
		let ct = c*t;
		let p = b + ct;
		let cp = cos(p);
		xyz += a*t2*cp;
		xyzd += a*t*(2.0*cp - ct*sin(p));
	    }

	    ph[k] = xyz;
	    vh[k] = xyzd / DJY;

	    // SSB to Sun, T^0 terms.
	    for j in 0..NS0[k] {
		let a = S0[0][j][k];
		let b = S0[1][j][k];
		let c = S0[2][j][k];
		let p = b + c*t;
		xyz += a*cos(p);
		xyzd -= a*c*sin(p);
	    }

	    // SSB to Sun, T^1 terms.
	    for j in 0..NS1[k] {
		let a = S1[0][j][k];
		let b = S1[1][j][k];
		let c = S1[2][j][k];
		let ct = c*t;
		let p = b + ct;
		let cp = cos(p);
		xyz += a*t*cp;
		xyzd += a*(cp - ct*sin(p));
	    }

	    // SSB to Sun, T^2 terms.
	    for j in 0..NS2[k] {
		let a = S2[0][j][k];
		let b = S2[1][j][k];
		let c = S2[2][j][k];
		let ct = c*t;
		let p = b + ct;
		let cp = cos(p);
		xyz += a*t2*cp;
		xyzd += a*t*(2.0*cp - ct*sin(p));
	    }

	    pb[k] = xyz;
	    vb[k] = xyzd / DJY;
	}

	// Rotate from ecliptic to ICRF coordinates and return the results.
	let mut x;
	let mut y;
	let mut z;

	x = ph[0];
	y = ph[1];
	z = ph[2];
	let ph = [     x + AM12*y + AM13*z,
		  AM21*x + AM22*y + AM23*z,
		  AM32*y + AM33*z];
	x = vh[0];
	y = vh[1];
	z = vh[2];
	let vh = [     x + AM12*y + AM13*z,
		  AM21*x + AM22*y + AM23*z,
		  AM32*y + AM33*z];
	x = pb[0];
	y = pb[1];
	z = pb[2];
	let pb = [      x + AM12*y + AM13*z,
	           AM21*x + AM22*y + AM23*z,
	                    AM32*y + AM33*z];
	x = vb[0];
	y = vb[1];
	z = vb[2];
	let vb = [       x + AM12*y + AM13*z,
	            AM21*x + AM22*y + AM23*z,
	                     AM32*y + AM33*z];

	EarthPosVel {
	    heliocentric:PosVel{ p:ph,v:vh },
	    barycentric:PosVel{ p:pb,v:vb },
	    warning
	}
    }
}

/// Earth rotation angle (IAU 2000 model)
///
/// Source: era00.for
pub fn rotation_angle(UT1((dj1,dj2)):UT1)->R {
    let (d1,d2) = if dj1 < dj2 { (dj1,dj2) } else { (dj2,dj1) };
    let t = d1 + ( d2 - DJ00 );
    let f = (d1 % 1.0) + (d2 % 1.0);
    
    anp(TWO_PI * ( f + 0.7790572732640 + 0.00273781191135448 * t))
}

/// Form the matrix of polar motion for a given date, IAU 2000.
///    XP,YP      coordinates of the pole (radians)
///    SP         the TIO locator s' (radians)
///
/// Source: pom00.for
pub fn polar_motion_matrix(xp:R,yp:R,sp:R)->Mat3 {
    Mat3::identity()
	.compose(&Mat3::rotation(0,-yp))
	.compose(&Mat3::rotation(1,-xp))
	.compose(&Mat3::rotation(2,sp))
}
