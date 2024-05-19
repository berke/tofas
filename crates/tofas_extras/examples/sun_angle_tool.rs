use anyhow::{Result,bail};
use pico_args::Arguments;

use tofas::{
    calendar::{GregorianDate,HMS},
    ellipsoid::Geodetic360,
};
use tofas_extras::sun_angle::{
    SunAngleParameters,
    SunAngleCalculator,
    SunAngleResultBundle
};

fn main()->Result<()> {
    let mut args = Arguments::from_env();

    let year : i32 = args.opt_value_from_str("--year")?.unwrap_or(2007);
    let month : i32 = args.opt_value_from_str("--month")?.unwrap_or(4);
    let day : i32 = args.opt_value_from_str("--day")?.unwrap_or(5);
    let hour : u8 = args.opt_value_from_str("--hour")?.unwrap_or(0);
    let minute : u8 = args.opt_value_from_str("--min")?.unwrap_or(0);
    let second : f64 = args.opt_value_from_str("--sec")?.unwrap_or(0.0);
    let lat : f64 = args.opt_value_from_str("--lat")?.unwrap_or(0.0);
    let lon : f64 = args.opt_value_from_str("--lon")?.unwrap_or(-120.0);
    let height : f64 = args.opt_value_from_str("--height")?.unwrap_or(0.0);

    let delta0 : f64 = args.opt_value_from_str("--delta0")?.unwrap_or(-600.0);
    let delta1 : f64 = args.opt_value_from_str("--delta1")?.unwrap_or(600.0);
    let delta_step : f64 = args.opt_value_from_str("--delta_step")?
	.unwrap_or(60.0);
    let scan = args.contains("--scan");
    if !args.finish().is_empty() {
	bail!("Unhandled extra arguments");
    }

    let date = GregorianDate::new(year,month,day)?;
    let time = HMS { hour,minute,second };
    let position = Geodetic360 { lat,lon,height };

    let parameters = SunAngleParameters {
	date,
	time,
	position
    };

    let calc = SunAngleCalculator::new(&parameters);

    if scan {
	println!("Will scan from {:+13.6}s to {:+13.6}s in steps of {:13.6}",
		 delta0,delta1,delta_step);
    }

    let mut delta = delta0;
    while delta <= delta1 {
	if scan {
	    println!();
	}
	let result = calc.compute(delta);
	let bundle = SunAngleResultBundle {
	    parameters:&parameters,
	    result:&result
	};
	print!("{}",bundle);
	if !scan {
	    break;
	}
	delta += delta_step;
    }

    Ok(())
}
