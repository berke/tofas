use std::fmt::Display;

use tofas::{
    common::*,
    time::{TT,TAI,TDB,UT1,UTC},
    frames,
    ellipsoid::{EllipsoidConverter,Geodetic,Geodetic360,WGS84},
    earth::{self,EarthPosVel},
    calendar::{GregorianDate,HMS},
    delta_at::DeltaAt,
};

#[derive(Clone,Debug)]
pub struct SunAngleParameters {
    pub date:GregorianDate,
    pub time:HMS,
    pub position:Geodetic360
}

impl Default for SunAngleParameters {
    fn default()->Self {
	Self {
	    date:GregorianDate {
		year:2007,
		month:4,
		day:5
	    },
	    time:HMS {
		hour:0,
		minute:0,
		second:0.0,
	    },
	    position:Geodetic360 {
		height:0.0,
		lat:0.0,
		lon:-120.0
	    }
	}
    }
}

pub struct SunAngleCalculator {
    jd0:f64,
    jd1_date:f64,
    dat:f64,
    fr:f64,
    dut1:f64,
    dtr:f64,
    xp:f64,
    yp:f64,
    p:[f64;3],
    zen:[f64;3]
}

#[derive(Clone,Debug)]
pub struct SunAngleResultBundle<'a,'b> {
    pub parameters:&'a SunAngleParameters,
    pub result:&'b SunAngleResult
}

#[derive(Clone,Debug)]
pub struct SunAngleResult {
    pub jd0:f64,
    pub jd1:f64,
    pub delta_s:f64,
    pub p:[f64;3],
    pub utc:UTC,
    pub dat:f64,
    pub tai:TAI,
    pub tt:TT,
    pub ut1:UT1,
    pub era:f64,
    pub earth:[f64;3],
    pub sun_e:[f64;3],
    pub sza:f64
}

impl<'a,'b> Display for SunAngleResultBundle<'a,'b> {
    fn fmt(&self,fmt:&mut std::fmt::Formatter<'_>)->Result<(),std::fmt::Error> {
	let &SunAngleParameters {
	    date,
	    time,
	    position
	} = self.parameters;
	let &SunAngleResult {
	    jd0,jd1,delta_s,p,utc,dat,tai,tt,ut1,era,earth,sun_e,sza,..
	} = self.result;
	let sza_d = sza / DEGREE;
	writeln!(fmt,"Date:                {}",date)?;
	writeln!(fmt,"Time:                {} {:+13.6}s",time,delta_s)?;
	writeln!(fmt,"Julian date:         {:.6} = {:.6} + {:.6}",
		 jd0 + jd1,jd0,jd1)?;
	writeln!(fmt,"Position:            {}",position)?;
	writeln!(fmt,"                     X={:16.1} Y={:16.1} Z={:16.1}",
		 p[0],p[1],p[2])?;
	writeln!(fmt,"UTC:                 {:.6}",utc.total())?;
	writeln!(fmt,"Delta AT:            {:.6}",dat)?;
	writeln!(fmt,"TAI:                 {:.6}",tai.total())?;
	writeln!(fmt,"TT:                  {:.6}",tt.total())?;
	writeln!(fmt,"UT1:                 {:.6}",ut1.total())?;
	writeln!(fmt,"ERA:                 {:.6}°",era/DEGREE)?;
	writeln!(fmt,"Earth position (AU): X={:+15.13} Y={:+15.13} Z={:+15.13}",
		 earth[0],earth[1],earth[2])?;
	writeln!(fmt,"Sun position (m):    X={:+16.0} Y={:+16.0} Z={:+16.0}",
		 sun_e[0],sun_e[1],sun_e[2])?;
	writeln!(fmt,"Sun Zenith angle:    {:7.2}° or elevation: {:7.2}°",
		 sza_d,90.0 - sza_d)?;
	Ok(())
    }
}

impl SunAngleCalculator {
    pub fn new(parameters:&SunAngleParameters)->Self {
	let &SunAngleParameters {
	    date,
	    time,
	    position
	} = parameters;

	// See example 5.1 in sofa_pn_f.pdf (p.18)
	let dtr = 0.0; // XXX

	let xp = 0.0349282 * AS2R; // Good for 2007
	let yp = 0.4833163 * AS2R;

	// UT1 - UTC
	let dut1 = -0.072073685;

	let fr = time.to_fraction_of_day();
	let (jd0,jd1_date) = date.to_julian();
	let dat = date.delta_at(fr).expect("Delta AT error").unwrap_or(0.0);

	let wgs84 = EllipsoidConverter::new(&WGS84).expect("Invalid ellipsoid");
	let p_gd : Geodetic = position.into();
	let p = wgs84.geodetic_to_geocentric(&p_gd).expect("Invalid position");
	let zen = p_gd.zenith();

	Self {
	    jd0,
	    jd1_date,
	    fr,
	    dat,
	    dut1,
	    dtr,
	    xp,
	    yp,
	    p,
	    zen
	}
    }

    pub fn compute(&self,delta_s:f64)->SunAngleResult {
	let &Self {
	    jd0,
	    jd1_date,
	    fr,
	    dat,
	    dut1,
	    dtr,
	    xp,
	    yp,
	    p,
	    zen
	} = self;

	let fr = fr + delta_s / 86400.0;

	let jd1 = jd1_date + fr;
	let utc = UTC((jd0,jd1));
	let tai = TAI::from_utc_delta_at(utc,dat);
	
	let tt : TT = tai.into();

	// let ut1 = UT1::from_tt(tt,dut1);
	let ut1 = UT1((jd0,jd1 + dut1/86400.0));
	let era = earth::rotation_angle(ut1);

	let tdb = TDB::from_tt(tt,dtr);

	let epv : EarthPosVel = tdb.into();
	let c2t = frames::celestial_to_terrestrial(tt,ut1,xp,yp);
	let earth = epv.heliocentric.p;
	let sun = earth.neg().scale(earth::AUM);
	let sun_e = c2t.transpose().apply(sun);

	let sza = zen.angle(sun_e);

	SunAngleResult {
	    jd0,
	    jd1,
	    delta_s,
	    p,
	    utc,
	    dat,
	    tt,
	    tai,
	    ut1,
	    era,
	    earth,
	    sun_e,
	    sza
	}
    }
}
