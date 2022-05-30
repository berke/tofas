#[test]
fn test_ellipsoid() {
    use crate::ellipsoid::*;
    use crate::common::*;

    let ec = EllipsoidConverter::new(&WGS84).unwrap();
    for _ in 0..1000 {
	let elong = (-180.0 + 360.0*fastrand::f64())*DEGREE;
	let phi = (-90.0 + 90.0*fastrand::f64())*DEGREE;
	let height = 10000.0*fastrand::f64();
	let gd = Geodetic{ elong,
			    phi,
			    height };
	let gc = ec.geodetic_to_geocentric(&gd).unwrap();
	let gd_bis = ec.geocentric_to_geodetic(&gc);
	let e = abs(gd.elong - gd_bis.elong) +
	    abs(gd.phi - gd_bis.phi) +
	    abs(gd.height - gd_bis.height);
	if e > sqrt(EPSILON) {
	    eprintln!("gd = {gd:#?}");
	    eprintln!("gc = {gc:#?}");
	    eprintln!("gd_bis = {gd_bis:#?}");
	    eprintln!("e = {e}");
	    panic!("Conversion error");
	}
    }
}

#[test]
fn test_calendar() {
    use crate::calendar::*;
    use crate::common::*;

    for (dj1,year,month,day,f1) in [
	(2436116.31,1957,10,4,0.81)
    ] {
	let gd1 = GregorianDate::new(year,month,day).unwrap();
	let (gd2,f2) = GregorianDate::from_julian(dj1,0.0).unwrap();
	let (dj2a,dj2b) = gd1.to_julian();
	let dj2 = dj2a + dj2b + f1;
	let e = dj1 - dj2;
	if e >= sqrt(EPSILON) {
	    eprintln!("gd1 = {gd1:?}, f1 = {f1}");
	    eprintln!("gd2 = {gd2:?}, f2 = {f2}");
	    eprintln!("dj2 = {dj2} = {dj2a} + {dj2b}, e = {e}");
	    panic!("Known conversion failed");
	}
    }

    for iter in 0..10000000 {
	let dj1a = floor(DJ_MIN + (MJD_ZERO + 36500.0) * fastrand::f64()) + 0.5;
	let f1 = fastrand::f64();
	let dj1 = dj1a + f1;
	let (gd,f2) = GregorianDate::from_julian(dj1a,f1).unwrap();
	let (dj2a,dj2b) = gd.to_julian();
	let dj2 = dj2a + dj2b + f2;
	let e = abs(dj1 - dj2);
	if e >= sqrt(EPSILON) {
	    eprintln!("Iteration {iter}");
	    eprintln!("dj1 = {dj1}");
	    eprintln!("gd = {gd:?}");
	    eprintln!("f1 = {f1:?}, f2 = {f2:?}");
	    eprintln!("dj2 = {dj2} = {dj2a} + {dj2b}, e = {e}");
	    panic!("Conversion failed");
	}
    }
}

#[test]
fn test_earth() {
    use crate::earth::{self,DJ00,EarthPositionAndVelocity};
    use crate::time::{TT,UT1,TDB};
    use crate::common::*;

    let dj0 = DJ00;
    let dj1 = dj0 + 36525.0;
    let dtr = 0.0; // XXX
    let dt = 0.0; // XXX
    for _iter in 0..100 {
	// TDB dates
	let tt1 = dj0;
	let tt2 = (dj1-dj0)*fastrand::f64();
	let tt = TT((tt1,tt2));
	let tdb = TDB::from_tt(tt,dtr);
	let ut1 = UT1::from_tt(tt,dt);
	let epv : EarthPositionAndVelocity = tdb.into();
	assert!(epv.warning.is_none());
	let a = earth::rotation_angle(ut1);
	println!("{:?} {:?} {}",epv.heliocentric,epv.barycentric,a);
    }

    let a0 = earth::rotation_angle(UT1::from_tt(TT((dj0,0.0)),dt));
    let a1 = earth::rotation_angle(UT1::from_tt(TT((dj0,1.0)),dt));
    println!("Rotation rate: {} degree/day",(a1 - a0)/DEGREE);

}
