mod sun_angle;

use sun_angle::{
    SunAngleParameters,
    SunAngleCalculator,
    SunAngleResultBundle
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
    let parameters = SunAngleParameters {
	year:g(1),
	month:g(2),
	day:g(3),
	hour:g(4),
	minute:g(5),
	second:g(6),
	lat:h(7),
	lon:h(8),
	height:h(9)
    };

    let mut calc = SunAngleCalculator::default();

    let result = calc.compute(&parameters);

    let bundle = SunAngleResultBundle {
	parameters:&parameters,
	result:&result
    };

    print!("{}",bundle);
}
