use crate::common::*;
use std::fmt::{Display,Formatter};

custom_error!{pub CalendarError
	      BadYear   = "bad year",
	      BadMonth  = "bad month",
	      BadDay    = "bad day",
	      BadJulian = "bad Julian day",
	      BadFract  = "bad fraction of day"
}

#[derive(Clone,Copy,Debug,PartialEq)]
pub struct GregorianDate {
    pub year:i32,
    pub month:i32,
    pub day:i32
}

#[derive(Clone,Copy,Debug)]
pub struct HMS {
    pub hour:u8,
    pub minute:u8,
    pub second:f64,
}

#[derive(Clone,Copy,Debug)]
pub struct GregorianDateHMS {
    pub date:GregorianDate,
    pub hms:HMS
}

pub const MJD_ZERO : R = 2400000.5;

pub const YEAR_MIN : i32 = -4799;
pub const DJ_MIN : R = -32000.0; // -68569.5;
pub const DJ_MAX : R = 1e9;
const MTAB : [i32;12] = [31,28,31,30,31,30,31,31,30,31,30,31];

impl GregorianDate {
    /// Construct a GregorianData, validating the year, month
    /// and day.
    ///
    /// Source: cal2jd.for
    pub fn new(year:i32,month:i32,day:i32)->Result<Self,CalendarError> {
	if year < YEAR_MIN {
	    return Err(CalendarError::BadYear);
	}

	if month < 1 || month > 12 {
	    return Err(CalendarError::BadMonth);
	}

	let mut ndays = MTAB[month as usize - 1];

	if month == 2 {
	    if year % 4 == 0 {
		ndays = 29
	    }
	    if year % 100 == 0 && year % 400 != 0 {
		ndays = 28;
	    }
	}

	if day < 1 || day > ndays {
	    return Err(CalendarError::BadDay);
	}

	Ok(Self{ year,month,day })
    }

    /// Gregorian Calendar to Julian date.
    ///
    /// Source: cal2jd.for
    pub fn to_julian(&self)->(R,R) {
	let &Self{ year,month,day } = self;
	let my = ( month - 14 ) / 12;
	let iypmy = year + my;
	let djm = ( 1461 * ( iypmy + 4800 ) ) / 4
	    + (  367 * ( month - 2 - 12*my ) ) / 12
	    - (    3 * ( ( iypmy + 4900 ) / 100 ) ) / 4
	    + day - 2432076;
	(MJD_ZERO,djm as R)
    }

    /// Julian Date to Gregorian year, month, day, and fraction of a day.
    ///
    /// Source: jd2cal.for
    pub fn from_julian(dj1:R,dj2:R)->Result<(Self,R),CalendarError> {
	let dj = dj1 + dj2;
	if dj < DJ_MIN || dj > DJ_MAX {
	    return Err(CalendarError::BadJulian);
	}

	let (d1,d2) =
	    if dj1 >= dj2 {
		(dj1,dj2)
	    } else {
		(dj2,dj1)
	    };
	let d2 = d2 - 0.5;

	let f1 = d1 % 1.0;
	let f2 = d2 % 1.0;
	let mut f = (f1 + f2) % 1.0;
	if f < 0.0 {
	    f += 1.0;
	}
	let d = round(d1 - f1) + round(d2 - f2) + round(f1 + f2 - f);
	let jd = d as i32 + 1;

	let l = jd + 68569;
	let n = ( 4*l ) / 146097;
	let l = l - ( 146097*n + 3 ) / 4;
	let i = ( 4000 * (l+1) ) / 1461001;
	let l = l - ( 1461*i ) / 4 + 31;
	let k = ( 80*l ) / 2447;
	let id = l - ( 2447*k ) / 80;
	let l = k / 11;
	let im = k + 2 - 12*l;
	let iy = 100 * ( n-49 ) + i + l;
	let fd = f;

	Ok((Self {
	    year:iy,
	    month:im,
	    day:id },
	    fd))
    }
}

impl Display for GregorianDate {
    fn fmt(&self,f:&mut Formatter<'_>)->Result<(),std::fmt::Error> {
	write!(f,"{:04}-{:02}-{:02}",self.year,self.month,self.day)
    }
}

impl HMS {
    pub fn new(hour:u8,minute:u8,second:f64)->Self {
	Self { hour,minute,second }
    }
    
    pub fn from_fraction_of_day(f:f64)->Result<Self,CalendarError> {
	if f < 0.0 || f > 1.0 {
	    return Err(CalendarError::BadFract)
	}
	let s = 86400.0 * f;
	let t = (s).trunc() as u32;
	let hour = (t / 3600) as u8;
	let minute = ((t % 3600) / 60) as u8;
	let second = (t % 60) as f64 + s.fract();
	Ok(Self {
	    hour,
	    minute,
	    second
	})
    }

    pub fn to_fraction_of_day(&self)->f64 {
	(self.second + 60.0*(self.minute as f64 + 60.0*self.hour as f64))/86400.0
    }
}

impl GregorianDateHMS {
    pub fn from_julian(dj1:R,dj2:R)->Result<GregorianDateHMS,CalendarError> {
	let (date,fod) = GregorianDate::from_julian(dj1,dj2)?;
	let hms = HMS::from_fraction_of_day(fod)?;
	Ok(Self { date,hms })
    }

    pub fn to_julian(&self)->(R,R) {
	let (dj1,dj2) = self.date.to_julian();
	let fod = self.hms.to_fraction_of_day();
	(dj1,dj2 + fod)
    }
}

impl Display for HMS {
    fn fmt(&self,f:&mut Formatter<'_>)->Result<(),std::fmt::Error> {
	write!(f,"{:02}:{:02}:{:0>9.6}",self.hour,self.minute,self.second)
    }
}

impl Display for GregorianDateHMS {
    fn fmt(&self,f:&mut Formatter<'_>)->Result<(),std::fmt::Error> {
	write!(f,"{} {}",self.date,self.hms)
    }
}
