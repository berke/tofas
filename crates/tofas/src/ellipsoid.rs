use std::fmt::{Display,Formatter};
use crate::common::*;
    
pub struct Ellipsoid {
    /// Equatorial radius [m]
    pub a:R,

    /// Flattening, around 0.00335 ~ 1/298 for Earth
    pub f:R
}

/// The WGS84 reference ellipsoid.
/// Source: eform.for
pub const WGS84 : Ellipsoid = Ellipsoid {
    a:6378137.0,
    f:1.0 / 298.257223563
};

/// A point expressed in geodetic coordinates
#[derive(Clone,Copy,Debug,PartialEq)]
pub struct Geodetic {
    /// Longitude, positive Eastwards [rad]
    pub elong:R,

    /// Geodetic latitude [rad]
    pub phi:R,

    /// Height above ellipsoid [m]
    pub height:R
}

impl Display for Geodetic {
    fn fmt(&self,f:&mut Formatter<'_>)->Result<(),std::fmt::Error> {
	// +23.456789,-111.111111,+12345.67
	write!(f,"{:+010.6},{:+011.6},{:+09.2}",
	       anp(self.phi)/DEGREE,
	       anp(self.elong)/DEGREE,
	       self.height)
    }
}

impl Geodetic {
    pub fn zenith(&self)->Vec3 {
	[cos(self.phi)*cos(self.elong),
	 cos(self.phi)*sin(self.elong),
	 sin(self.phi)
	]
    }
}

custom_error!{pub EllipsoidError
	      InvalidRadius        = "invalid radius",
	      InvalidFlattening    = "invalid flattening value",
	      BadEccentricity      = "could not compute EC2",
	      BadCoordinates       = "bad coordinates"
}

pub struct EllipsoidConverter {
    a:R,
    f:R,
    aeps2:R,
    e2:R,
    e4t:R,
    ec2:R,
    ec:R,
    b:R
}

impl EllipsoidConverter {
    pub fn new(ellipsoid:&Ellipsoid)->Result<Self,EllipsoidError> {
	let &Ellipsoid{ a,f } = ellipsoid;

	if f < 0.0 || f >= 1.0 {
	    return Err(EllipsoidError::InvalidFlattening);
	}

	if a < 0.0 {
	    return Err(EllipsoidError::InvalidRadius);
	}

	let aeps2 = a*a*1e-32;
	let e2 = (2.0 - f)*f;
	let e4t = e2*e2*1.5;
	let ec2 = 1.0 - e2;
	if ec2 <= 0.0 {
	    return Err(EllipsoidError::BadEccentricity);
	}
	let ec = sqrt(ec2);
	let b = a*ec;
	Ok(Self{
	    a,
	    f,
	    aeps2,
	    e2,
	    e4t,
	    ec2,
	    ec,
	    b
	})
    }

    /// Transform geocentric coordinates to geodetic.
    ///
    /// Based on: gc2gde.for
    pub fn geocentric_to_geodetic(&self,xyz:&Vec3)->Geodetic {
	let &Self { a,aeps2,e2,e4t,ec2,ec,b,.. } = self;
	let &[x,y,z] = xyz;
	let p2 = sq(x) + sq(y);
	let elong =
	    if p2 > 0.0 {
		atan2(y,x)
	    } else {
		0.0
	    };

	let absz = abs(z);

	let (phi,height) =
	    if p2 > aeps2 {
		let p = sqrt(p2);
		let s0 = absz/a;
		let pn = p/a;
		let zc = ec*s0;

		let c0 = ec*pn;
		let c02 = c0*c0;
		let c03 = c02*c0;
		let s02 = s0*s0;
		let s03 = s02*s0;
		let a02 = c02+s02;
		let a0 = sqrt(a02);
		let a03 = a02*a0;
		let d0 = zc*a03 + e2*s03;
		let f0 = pn*a03 - e2*c03;

		let b0 = e4t*s02*c02*pn*(a0-ec);
		let s1 = d0*f0 - b0*s0;
		let cc = ec*(f0*f0-b0*c0);

		let phi = atan(s1/cc);
		let s12 = s1*s1;
		let cc2 = cc*cc;
		let height = (p*cc+absz*s1-a*sqrt(ec2*s12+cc2))/sqrt(s12+cc2);

		(phi,height)
	    } else {
		(PI/2.0,absz - b)
	    };

	let phi = if z < 0.0 { -phi } else { phi };

	Geodetic{ elong,phi,height }
    }

    /// Transform geodetic coordinates to geocentric.
    ///
    /// Based on: gd2cde.for
    pub fn geodetic_to_geocentric(&self,gd:&Geodetic)->Result<Vec3,EllipsoidError> {
	let &Self { a,f,.. } = self;
	let &Geodetic{ elong,phi,height } = gd;
	let sp = sin(phi);
	let cp = cos(phi);
	let w = 1.0 - f;
	let w = w*w;
	let d = cp*cp + w*sp*sp;
	if d > 0.0 {
	    let ac = a / sqrt(d);
	    let wac = w * ac;
	    let r = ( ac + height ) * cp;
	    let xyz = [
		r * cos(elong),
		r * sin(elong),
		( wac + height ) * sp
	    ];
	    Ok(xyz)
	} else {
	    Err(EllipsoidError::BadCoordinates)
	}
    }
}
