use anyhow::{Result,bail};
use pico_args::Arguments;

use tofas::calendar::{GregorianDate,HMS};

fn gregorian(mut args:Arguments)->Result<()> {
    let year : i32 = args.value_from_str("--year")?;
    let month : i32 = args.value_from_str("--month")?;
    let day : i32 = args.value_from_str("--day")?;
    let hour : u8 = args.opt_value_from_str("--hour")?.unwrap_or(12);
    let minute : u8 = args.opt_value_from_str("--minute")?.unwrap_or(0);
    let second : f64 = args.opt_value_from_str("--second")?.unwrap_or(0.0);
    let gd = GregorianDate::new(year,month,day)?;
    println!("Gregorian date: {}",gd);
    let hms = HMS::new(hour,minute,second);
    println!("Time: {}",hms);
    let fd = hms.to_fraction_of_day();
    let (jd1,jd2) = gd.to_julian();
    println!("Julian date: {:.6}",jd1 + jd2 + fd);
    println!("Fraction of day: {:.6}",fd);
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
