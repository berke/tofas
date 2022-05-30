use tofas::{
    common::*,
    time::{TT,TAI,TDB,UT1},
    frames,
    ellipsoid::{EllipsoidConverter,Geodetic,WGS84},
    earth::{self,EarthPosVel},
    calendar::GregorianDate
};

fn main() {
    let args : Vec<String> = std::env::args().collect();
    let f = |i:usize| {
	if i >= args.len() {
	    panic!("Need at least {i} arguments");
	}
	&args[i]
    };
    let g = |i:usize| {
	f(i).parse::<i32>().unwrap_or_else(|_| panic!("Invalid argument #{i}"))
    };
    let h = |i:usize| {
	f(i).parse::<f64>().unwrap_or_else(|_| panic!("Invalid argument #{i}"))
    };
    let year = g(1);
    let month = g(2);
    let day = g(3);
    let hour = g(4);
    let minute = g(5);
    let second = g(6);
    let lat = h(7);
    let lon = h(8);
    let height = h(9);

    // See example 5.1 in sofa_pn_f.pdf (p.18)
    let dtr = 0.0; // XXX

    let xp = 0.0349282 * AS2R; // Good for 2007
    let yp = 0.4833163 * AS2R;

    // UT1 - UTC
    let dut1 = -0.072073685;

    let fr = (second as f64 + 60.0*(minute as f64 + 60.0*hour as f64))/86400.0;
    let gd = GregorianDate::new(year,month,day).expect("Invalid date");
    let (jd0,jd1) = gd.to_julian();
    let jd1 = jd1 + fr;

    let wgs84 = EllipsoidConverter::new(&WGS84).expect("Invalid ellipsoid");
    let p_gd = Geodetic{ elong:lon*DEGREE, phi:lat*DEGREE, height };
    let p = wgs84.geodetic_to_geocentric(&p_gd).expect("Invalid position");
    let zen = p_gd.zenith();

    // Now pretend this is TAI
    let tai = TAI((jd0,jd1));
    let tt : TT = tai.into();

    let ut1 = UT1::from_tt(tt,dut1);
    let era = earth::rotation_angle(ut1);

    let tdb = TDB::from_tt(tt,dtr);

    let epv : EarthPosVel = tdb.into();
    let c2t = frames::celestial_to_terrestrial(tt,ut1,xp,yp);
    let earth = epv.heliocentric.p;
    let sun = earth.neg().scale(earth::AUM);
    let sun_e = c2t.transpose().apply(sun);

    let sza = zen.angle(sun_e);
    let sza_d = sza/DEGREE;

    println!("Date:         {year:04}-{month:02}-{day:02}");
    println!("Time:         {hour:02}:{minute:02}:{second:02}");
    println!("Julian date:  {:.6} = {:.6} + {:.6}",jd0 + jd1,jd0,jd1);
    println!("Position:     {lat:9.4} N, {lon:9.4} E, height {height:6.2} m");
    println!("Position:     X={:16.1} Y={:16.1} Z={:16.1}",p[0],p[1],p[2]);
    println!("TT:           {:.6}",tt.total());
    println!("UT1:          {:.6}",ut1.total());
    println!("ERA:          {:.6}°",era/DEGREE);
    println!("Sun position: X={} Y={} Z={}",
	     earth[0],earth[1],earth[2]);
    println!("Sun position: X={:16.1} Y={:16.1} Z={:16.1}",
	     sun_e[0],sun_e[1],sun_e[2]);
    println!("Sun Zenith angle: {:7.2}° or elevation: {:7.2}°",sza_d,90.0 - sza_d);
}
