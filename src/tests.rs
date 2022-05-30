use crate::{
    common::*,
    earth::{self,EarthPosVel},
    time::{TT,UT1,TDB,DJ00},
    ellipsoid::*,
    calendar::*,
    frames,
    test_data::*
};

const S00_ZERO : bool = false;

#[test]
fn test_ellipsoid() {

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
    let xp = 0.0;
    let yp = 0.0;
    let dj0 = DJ00;
    let dtr = 0.0; // XXX
    let dt = 0.0; // XXX

    for iday in 0..1 {
	for ihour in 0..24 {
	    // TDB dates
	    let tt1 = dj0;
	    let tt2 = iday as R + ihour as R / 24.0;
	    let tt = TT((tt1,tt2));
	    let tdb = TDB::from_tt(tt,dtr);
	    let ut1 = UT1::from_tt(tt,dt);
	    let epv : EarthPosVel = tdb.into();
	    assert!(epv.warning.is_none());
	    // let a = earth::rotation_angle(ut1);
	    // println!("{:?} {:?} {}",epv.heliocentric,epv.barycentric,a);
	    let c2t = frames::celestial_to_terrestrial(tt,ut1,xp,yp);
	    // println!("Barycentric coordinates of Earth: {:?}",epv.barycentric.p);
	    let p = epv.barycentric.p;
	    let sun = p.neg(); // scale(-earth::AUM);
	    let sun_e = c2t.apply(sun);
	    println!("{iday:3} {ihour:2} {} {} {}",sun_e[0],sun_e[1],sun_e[2]);
	}
    }

    let a0 = earth::rotation_angle(UT1::from_tt(TT((dj0,0.0)),dt));
    let a1 = earth::rotation_angle(UT1::from_tt(TT((dj0,1.0)),dt));
    println!("Rotation rate: {} degree/day",(a1 - a0)/DEGREE);
}

fn compare_matrices(name:&str,a:&Mat3,b:&Mat3,tol:R) {
    for i in 0..3 {
	for j in 0..3 {
	    let e = abs(a[i][j] - b[i][j]);
	    if e > tol {
		println!("MISMATCH IN {name}");
		println!("{:#?}",a);
		println!("vs");
		println!("{:#?}",b);
		println!("At indices {i},{j} : {} vs {}, error {e:e} > {tol:e}",a[i][j],b[i][j]);
		panic!("Mismatch in {name} indices {i},{j}");
	    }
	}
    }
}

fn compare_numbers(name:&str,a:R,b:R,tol:R) {
    let e = abs(a - b);
    if e > tol {
	panic!("Mismatch in {name} got {a} vs {b} error {e} > {tol}")
    }
}

#[test]
fn test_polar_motion_matrix() {
    let tol = EPSILON;
    for &(xp,yp,sp,rpom1) in POM00_DATA {
	let rpom2 = earth::polar_motion_matrix(xp,yp,sp);
	compare_matrices("POM00",&rpom1,&rpom2,tol);
    }
}

#[test]
fn test_rotation() {
    for _ in 0..100 {
	let theta0 = fastrand::f64()*TWO_PI;
	let theta1 = fastrand::f64()*TWO_PI;
	let theta2 = fastrand::f64()*TWO_PI;
	let rot0 = Mat3::rotation(0,theta0);
	let rot1 = Mat3::rotation(0,theta1);
	let rot2 = Mat3::rotation(0,theta2);
	let a = rot0.compose(&rot1.compose(&rot2));
	let irot0 = Mat3::rotation(0,-theta0);
	let irot1 = Mat3::rotation(0,-theta1);
	let irot2 = Mat3::rotation(0,-theta2);
	let ia = irot2.compose(&irot1.compose(&irot0));
	let id1 = a.compose(&ia);
	let id2 = Mat3::identity();
	compare_matrices("ROT",&id1,&id2,sqrt(EPSILON));
    }
}

#[test]
fn test_nutation_matrix() {
    for &(epsa,dpsi,deps,rmatn1) in NUMAT_DATA {
	let rmatn2 = frames::nutation_matrix(epsa,dpsi,deps);
	compare_matrices("rmatn",&rmatn1,&rmatn2,EPSILON);
    }
}

#[test]
fn test_nutation() {
    for &(date1,date2,dpsi1,deps1) in NUT00B_DATA {
	let (dpsi2,deps2) = frames::nutation(TT((date1,date2)));
	compare_numbers("dpsi",dpsi1,dpsi2,EPSILON);
	compare_numbers("deps",deps1,deps2,EPSILON);
    }
}

fn sel<T:Copy>(c:bool,x1:T,x0:T)->T { if c { x1 } else { x0 } }

#[test]
fn test_bias_and_precession() {
    let tol = sel(S00_ZERO,EPSILON,2.0*EPSILON);
    for &(date1,date2,[rb1,rp1,rbp1]) in BP00_DATA {
	let (rb2,rp2,rbp2) = frames::bias_and_precession(TT((date1,date2)));
	compare_matrices("RB",&rb1,&rb2,tol);
	compare_matrices("RP",&rp1,&rp2,tol);
	compare_matrices("RBP",&rbp1,&rbp2,tol);
    }
}

#[test]
fn test_rotation_angle() {
    let tol = EPSILON;
    for &(ut11,ut12,era1) in ERA00_DATA {
	let era2 = earth::rotation_angle(UT1((ut11,ut12)));
	compare_numbers("era",era1,era2,tol);
    }
}

#[test]
fn test_precession_rate() {
    let tol = EPSILON;
    for &(tt1,tt2,dpsipr1,depspr1) in PR00_DATA {
	let (dpsipr2,depspr2) = frames::precession_rate(TT((tt1,tt2)));
	compare_numbers("dpsipr",dpsipr1,dpsipr2,tol);
	compare_numbers("depspr",depspr1,depspr2,tol);
    }
}

#[test]
fn test_precession_nutation() {
    let tol = 2.0*EPSILON;
    for &(tt1,tt2,rbpn_1) in PNM00B_DATA {
	let tt = TT((tt1,tt2));
	let rbpn_2 = frames::precession_nutation(tt);
	compare_matrices("RBPN",&rbpn_1,&rbpn_2,tol);
    }
}

#[test]
fn test_celestial_to_intermediate() {
    let tol = sel(S00_ZERO,EPSILON,1e-6);
    for &(tt1,tt2,rc2i_1) in C2I00B_DATA {
	let tt = TT((tt1,tt2));
	let rc2i_2 = frames::celestial_to_intermediate(tt);
	compare_matrices("RC2I",&rc2i_1,&rc2i_2,tol);
    }
}

#[test]
fn test_celestial_to_terrestrial_from_cio_components() {
    let tol = EPSILON;
    for &(rc2i,era,rpom,rc2t1) in C2TCIO_DATA {
	let rc2t2 = frames::celestial_to_terrestrial_from_cio_components(&rc2i,era,&rpom);
	compare_matrices("RC2T",&rc2t1,&rc2t2,tol);
    }
}

#[test]
fn test_celestial_to_terrestrial() {
    let tol = EPSILON;
    let tol_mat = sel(S00_ZERO,EPSILON,1e-6);
    let xp = 0.0;
    let yp = 0.0;
    for &(tt1,tt2,dtr,tdb1_1,tdb2_1,dt,ut11_1,ut12_1,rc2t_1) in C2T00B_DATA {
	let tt = TT((tt1,tt2));
	let TDB((tdb1_2,tdb2_2)) = TDB::from_tt(tt,dtr);
	let ut1_2 @ UT1((ut11_2,ut12_2)) = UT1::from_tt(tt,dt);
	compare_numbers("tdb1",tdb1_1,tdb1_2,tol);
	compare_numbers("tdb2",tdb2_1,tdb2_2,tol);
	compare_numbers("ut11",ut11_1,ut11_2,tol);
	compare_numbers("ut12",ut12_1,ut12_2,tol);
	let rc2t_2 = frames::celestial_to_terrestrial(tt,ut1_2,xp,yp);
	compare_matrices("RC2T",&rc2t_1,&rc2t_2,tol_mat);
    }
}
