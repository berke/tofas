use crate::common::*;

/// GPS -> TAI -> TT -> TDB -> EPV00
/// GPS -> TAI -> TT -> UT1 -> ERA00

#[derive(Copy,Clone,Debug)]
pub struct UT1(pub (R,R));

#[derive(Copy,Clone,Debug)]
pub struct TT(pub (R,R));

#[derive(Copy,Clone,Debug)]
pub struct TAI(pub (R,R));

// #[derive(Copy,Clone,Debug)]
// pub struct UTC((R,R));

#[derive(Copy,Clone,Debug)]
pub struct TDB(pub (R,R));

const D2S : R = 86400.0;

impl TT {
    /// Universal Time, UT1, to Terrestrial Time, TT
    ///
    /// Source: ut1tt.for
    fn from_ut1(UT1((ut11,ut12)):UT1,dt:R)->Self {
	let dtd = dt/D2S;
	let (tt1,tt2) =
	    if ut11 > ut12 {
		(ut11,ut12 + dtd)
	    } else {
		(ut11 + dtd,ut12)
	    };
	Self((tt1,tt2))
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
}
