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

pub const PI:R = std::f64::consts::PI;
pub const DEGREE:R = PI/180.0;
pub const EPSILON:R = R::EPSILON;
