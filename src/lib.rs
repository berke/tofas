mod common;
pub mod ellipsoid;
    
#[cfg(test)]
mod tests {
    #[test]
    fn test_ellipsoid() {
	use crate::ellipsoid::*;
	use crate::common::*;
	let ec = EllipsoidConverter::new(&WGS84).unwrap();
	let gd = Geodetic{ elong:1.0*DEGREE,
			   phi:43.21*DEGREE,
			   height:1234.5 };
	let gc = ec.geodetic_to_geocentric(&gd).unwrap();
	let gd_bis = ec.geocentric_to_geodetic(&gc);
	let e = abs(gd.elong - gd_bis.elong) +
	        abs(gd.phi - gd_bis.phi) +
	        abs(gd.height - gd_bis.height);
	println!("gd = {gd:#?}");
	println!("gc = {gc:#?}");
	println!("gd_bis = {gd_bis:#?}");
	println!("e = {e}");
	assert!(e < sqrt(EPSILON));
    }
}
