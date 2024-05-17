use anyhow::{Result,bail};
use pico_args::Arguments;

use tofas::calendar::{GregorianDate,HMS};

fn gregorian(mut args:Arguments)->Result<()> {
    let year : i32 = args.value_from_str("--year")?;
    let month : i32 = args.value_from_str("--month")?;
    let day : i32 = args.value_from_str("--day")?;
    let gd = GregorianDate::new(year,month,day)?;
    println!("Gregorian date: {}",gd);
    let (jd1,jd2) = gd.to_julian();
    println!("Julian date: {:.1}",jd1 + jd2);
    Ok(())
}

fn julian(mut args:Arguments)->Result<()> {
    let jd1 : f64 = args.value_from_str("--jd")?;
    let jd2 : f64 = args.opt_value_from_str("--jd2")?.unwrap_or(0.0);
    if !args.finish().is_empty() {
	bail!("Invalid extra arguments");
    }
    println!("Julian date: {:18.14}",jd1 + jd2);
    let (gd,fod) = GregorianDate::from_julian(jd1,jd2)?;
    println!("Gregorian date: {}",gd);
    println!("Fraction of day: {}",fod);
    let hms = HMS::from_fraction_of_day(fod)?;
    println!("Time: {}",hms);
    Ok(())
}

fn main()->Result<()> {
    let mut args = Arguments::from_env();

    match args.subcommand()?.as_ref().map(|x| x.as_str()) {
	Some("gregorian"|"g") => gregorian(args),
	Some("julian"|"j") => julian(args),
	_ => bail!("Specify gregorian or julian")
    }
}
