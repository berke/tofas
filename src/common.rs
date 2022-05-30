pub use num::traits::real::Real;
pub use custom_error::custom_error;

pub type R = f64;

pub fn sqrt<T:Real>(x:T)->T { x.sqrt() }
pub fn sq<T:Real>(x:T)->T { x*x }
pub fn sin<T:Real>(x:T)->T { x.sin() }
pub fn cos<T:Real>(x:T)->T { x.cos() }
pub fn atan<T:Real>(x:T)->T { x.atan() }
pub fn atan2<T:Real>(y:T,x:T)->T { y.atan2(x) }
pub fn abs<T:Real>(x:T)->T { x.abs() }
pub fn floor<T:Real>(x:T)->T { x.floor() }
pub fn round<T:Real>(x:T)->T { x.round() }

pub const PI : R = std::f64::consts::PI;
pub const DEGREE : R = PI/180.0;
pub const EPSILON : R = R::EPSILON;

pub const TWO_PI : R = 6.283185307179586476925287;

pub fn anp(a:R)->R {
    // Normalize angle to range [0,2π[
    let mut w = a % TWO_PI;
    if w < 0.0 {
	w += TWO_PI;
    }
    w
}
