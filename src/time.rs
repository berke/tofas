use crate::common::*;

/// GPS -> TAI -> TT -> TDB -> EPV00
/// GPS -> TAI -> TT -> UT1 -> ERA00

#[derive(Copy,Clone,Debug)]
pub struct UT1(pub (R,R));

#[derive(Copy,Clone,Debug)]
pub struct TT(pub (R,R));

#[derive(Copy,Clone,Debug)]
pub struct TAI(pub (R,R));

#[derive(Copy,Clone,Debug)]
pub struct UTC(pub (R,R));

#[derive(Copy,Clone,Debug)]
pub struct TDB(pub (R,R));

pub const D2S : R = 86400.0;
pub const DJY : R = 365.25;
pub const DJ00 : R = 2451545.0;
pub const DJC : R = 36525.0;

impl TT {
    /// Universal Time, UT1, to Terrestrial Time, TT
    ///
    /// Source: ut1tt.for
    pub fn from_ut1(UT1((ut11,ut12)):UT1,dt:R)->Self {
	let dtd = dt/D2S;
	let (tt1,tt2) =
	    if ut11 > ut12 {
		(ut11,ut12 + dtd)
	    } else {
		(ut11 + dtd,ut12)
	    };
	Self((tt1,tt2))
    }

    pub fn total(self)->R {
	let TT((tt1,tt2)) = self;
	tt1 + tt2
    }
}

impl UTC {
    pub fn total(self)->R {
	let UTC((utc1,utc2)) = self;
	utc1 + utc2
    }
}

impl TAI {
    pub fn from_utc_delta_at(UTC((utc1,utc2)):UTC,dat:f64)->Self {
	let da = dat/86400.0;
	let (tai1,tai2) =
	    if utc1 > utc2 {
		(utc1,utc2 + da)
	    } else {
		(utc1 + da,utc2)
	    };
	Self((tai1,tai2))
    }
    
    pub fn total(self)->R {
	let TAI((tai1,tai2)) = self;
	tai1 + tai2
    }
}

pub const DTAT : R = 32.184/86400.0;

impl From<TAI> for TT {
    fn from(TAI((tai1,tai2)):TAI)->Self {
	let (tt1,tt2) =
	    if tai1 > tai2 {
		(tai1,tai2 + DTAT)
	    } else {
		(tai1 + DTAT,tai2)
	    };
	Self((tt1,tt2))
    }
}

impl TDB {
    /// Time scale transformation: Terrestrial Time, TT, to
    /// Barycentric Dynamical Time, TDB.
    ///
    /// Source: tttdb.for
    pub fn from_tt(TT((tt1,tt2)):TT,dtr:R)->Self {
	let dtrd = dtr/D2S;
	let (tdb1,tdb2) =
	    if tt1 > tt2 {
		(tt1,tt2 + dtrd)
	    } else {
		(tt1 + dtrd,tt2)
	    };
	Self((tdb1,tdb2))
    }
}

impl UT1 {
    /// Time scale transformation:  Terrestrial Time, TT, to Universal Time,
    /// UT1.
    ///
    /// The argument [dt] is classical Delta T.
    ///
    /// Source: ttut1.for
    pub fn from_tt(TT((tt1,tt2)):TT,dt:R)->Self {
	let dtd = dt/D2S;
	let (ut11,ut12) =
	    if tt1 > tt2 {
		(tt1,tt2 - dtd)
	    } else {
		(tt1 - dtd,tt2)
	    };
	Self((ut11,ut12))
    }

    pub fn total(self)->R {
	let UT1((ut11,ut12)) = self;
	ut11 + ut12
    }
}
