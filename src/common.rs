pub use num::traits::real::Real;
pub use custom_error::custom_error;

pub type R = f64;
pub type Vec3 = [R;3];
pub type Mat3 = [Vec3;3];

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
pub const DAS2R : R = 4.848136811095359935899141e-6;
pub const DMAS2R : R = DAS2R / 1e3;
pub const TURNAS : R = 1296000.0;

pub fn anp(a:R)->R {
    // Normalize angle to range [0,2π[
    let mut w = a % TWO_PI;
    if w < 0.0 {
	w += TWO_PI;
    }
    w
}

pub trait Vector {
    fn scale(self,c:R)->Self;
    fn add(self,b:Self)->Self;
    fn sub(self,b:Self)->Self;
    fn neg(self)->Self;
}

impl Vector for Vec3 {
    fn neg(self)->Self {
	[-self[0],
	 -self[1],
	 -self[2]]
    }

    fn scale(self,c:R)->Self {
	[c*self[0],
	 c*self[1],
	 c*self[2]]
    }
    
    fn add(self,b:Self)->Self {
	[self[0] + b[0],
	 self[1] + b[1],
	 self[2] + b[2]]
    }

    fn sub(self,b:Self)->Self {
	[self[0] - b[0],
	 self[1] - b[1],
	 self[2] - b[2]]
    }
}

pub trait Matrix {
    type Vector;
    fn zero()->Self;
    fn identity()->Self;
    fn apply(&self,x:Self::Vector)->Self::Vector;
    fn compose(&self,b:&Self)->Self;
    fn transpose(&self)->Self;
    fn rotation(axis:usize,theta:R)->Self;
}

impl Matrix for Mat3 {
    type Vector = Vec3;

    fn zero()->Self {
	[[0.0,0.0,0.0],
	 [0.0,0.0,0.0],
	 [0.0,0.0,0.0]]
    }
    
    fn identity()->Self {
	[[1.0,0.0,0.0],
	 [0.0,1.0,0.0],
	 [0.0,0.0,1.0]]
    }

    fn apply(&self,x:Vec3)->Vec3 {
	[self[0][0]*x[0] + self[1][0]*x[1] + self[2][0]*x[2],
	 self[0][1]*x[0] + self[1][1]*x[1] + self[2][1]*x[2],
	 self[0][2]*x[0] + self[1][2]*x[1] + self[2][2]*x[2]]
    }

    fn compose(&self,b:&Self)->Self {
	let mut c = Self::zero();
	for i in 0..3 {
	    for j in 0..3 {
		let mut s = 0.0;
		for k in 0..3 {
		    s += self[i][k]*b[k][j];
		}
		c[i][j] = s;
	    }
	}
	c
    }

    fn transpose(&self)->Self {
	let mut c = Self::zero();
	for i in 0..3 {
	    for j in 0..3 {
		c[i][j] = self[j][i];
	    }
	}
	c
    }

    fn rotation(axis:usize,theta:R)->Self {
	let c = cos(theta);
	let s = sin(theta);
	let i0 = (axis + 1) % 3;
	let i1 = (axis + 2) % 3;
	let mut r = Self::zero();;
	r[i0][i0] = c;
	r[i0][i1] = s;
	r[i1][i0] = -s;
	r[i1][i1] = c;
	r[axis][axis] = 1.0;
	r
    }
}
