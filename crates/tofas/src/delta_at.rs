use crate::common::*;
use crate::calendar::GregorianDate;

// Calculate delta(AT) = TAI - UTC
//
// Source: dat.for

#[derive(Debug)]
struct DatData {
    year:i32,
    month:i32,
    dat:f64,
    drift:Option<[f64;2]>
}

const DAT_DATA : &[DatData] = &[
    DatData { year:1960,month:1,dat:1.417818,drift:Some([37300.0,0.001296]) },
    DatData { year:1961,month:1,dat:1.422818,drift:Some([37300.0,0.001296]) },
    DatData { year:1961,month:8,dat:1.372818,drift:Some([37300.0,0.001296]) },
    DatData { year:1962,month:1,dat:1.845858,drift:Some([37665.0,0.0011232]) },
    DatData { year:1963,month:11,dat:1.945858,drift:Some([37665.0,0.0011232]) },
    DatData { year:1964,month:1,dat:3.240130,drift:Some([38761.0,0.001296]) },
    DatData { year:1964,month:4,dat:3.340130,drift:Some([38761.0,0.001296]) },
    DatData { year:1964,month:9,dat:3.440130,drift:Some([38761.0,0.001296]) },
    DatData { year:1965,month:1,dat:3.540130,drift:Some([38761.0,0.001296]) },
    DatData { year:1965,month:3,dat:3.640130,drift:Some([38761.0,0.001296]) },
    DatData { year:1965,month:7,dat:3.740130,drift:Some([38761.0,0.001296]) },
    DatData { year:1965,month:9,dat:3.840130,drift:Some([38761.0,0.001296]) },
    DatData { year:1966,month:1,dat:4.313170,drift:Some([39126.0,0.002592]) },
    DatData { year:1968,month:2,dat:4.213170,drift:Some([39126.0,0.002592]) },
    DatData { year:1972,month:1,dat:10.0,drift:None },
    DatData { year:1972,month:7,dat:11.0,drift:None },
    DatData { year:1973,month:1,dat:12.0,drift:None },
    DatData { year:1974,month:1,dat:13.0,drift:None },
    DatData { year:1975,month:1,dat:14.0,drift:None },
    DatData { year:1976,month:1,dat:15.0,drift:None },
    DatData { year:1977,month:1,dat:16.0,drift:None },
    DatData { year:1978,month:1,dat:17.0,drift:None },
    DatData { year:1979,month:1,dat:18.0,drift:None },
    DatData { year:1980,month:1,dat:19.0,drift:None },
    DatData { year:1981,month:7,dat:20.0,drift:None },
    DatData { year:1982,month:7,dat:21.0,drift:None },
    DatData { year:1983,month:7,dat:22.0,drift:None },
    DatData { year:1985,month:7,dat:23.0,drift:None },
    DatData { year:1988,month:1,dat:24.0,drift:None },
    DatData { year:1990,month:1,dat:25.0,drift:None },
    DatData { year:1991,month:1,dat:26.0,drift:None },
    DatData { year:1992,month:7,dat:27.0,drift:None },
    DatData { year:1993,month:7,dat:28.0,drift:None },
    DatData { year:1994,month:7,dat:29.0,drift:None },
    DatData { year:1996,month:1,dat:30.0,drift:None },
    DatData { year:1997,month:7,dat:31.0,drift:None },
    DatData { year:1999,month:1,dat:32.0,drift:None },
    DatData { year:2006,month:1,dat:33.0,drift:None },
    DatData { year:2009,month:1,dat:34.0,drift:None },
    DatData { year:2012,month:7,dat:35.0,drift:None },
    DatData { year:2015,month:7,dat:36.0,drift:None },
    DatData { year:2017,month:1,dat:37.0,drift:None },
];

custom_error! { pub DeltaAtError
		BadFract  = "bad fraction of day" }

pub trait DeltaAt {
    fn delta_at(&self,fd:f64)->Result<Option<f64>,DeltaAtError>;
}

impl DeltaAt for GregorianDate {
    // Returns None for a pre-UTC year
    fn delta_at(&self,fd:f64)->Result<Option<f64>,DeltaAtError> {
	let iy = self.year;
	let im = self.month;
	
	if fd < 0.0 || fd > 1.0 {
	    return Err(DeltaAtError::BadFract);
	}

	if iy < DAT_DATA[0].year {
	    return Ok(None);
	}

	let m = 12*iy + im;

	let idat =
	    DAT_DATA
	    .iter()
	    .rposition(|DatData { year,month,.. }| 12*year + month < m)
	    .unwrap_or(0);

	let DatData { dat,drift,.. } = DAT_DATA[idat];

	let mut da = dat;

	if let Some(d) = drift {
	    let (_djm0,djm) = self.to_julian();
	    da += (djm + fd - d[0])*d[1];
	}

	Ok(Some(da))
    }
}
