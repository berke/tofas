use crate::common::*;

/// Mean anomaly of the Moon (IERS Conventions 2003).
/// Source: fal03.for
pub fn l03(t:R)->R {
    (485868.249036 +
     t*(1717915923.2178 +
	t*(31.8792 +
	   t*(0.051635 +
	      t*(-0.00024470)))) % TURNAS) * AS2R
}

/// Mean anomaly of the Sun (IERS Conventions 2003).
/// Source: falp03.for
pub fn lp03(t:R)->R {
    (1287104.793048 +
      t*(129596581.0481 +
	  t*(- 0.5532 +
		t*(0.000136 +
		    t*(-0.00001149)))) % TURNAS) * AS2R
}

/// Mean longitude of the Moon minus mean longitude of the ascending
/// node.
/// Source: faf03.for
pub fn f03(t:R)->R {
    (335779.526232 +
     t*(1739527262.8478 +
	t*(-12.7512 +
	   t*(-0.001037 +
	      t*(0.00000417 )))) % TURNAS) * AS2R
}

/// Mean elongation of the Moon from the Sun.
/// Source: fad03.for
pub fn d03(t:R)->R {
    (1072260.703692 +
     t*(1602961601.2090 +
	t*(-6.3706 +
	   t*(0.006593 +
	      t*(-0.00003169)))) % TURNAS) * AS2R
}

/// Mean longitude of the Moon's ascending node.
/// Source: faom03.for
pub fn om03(t:R)->R {
    (450160.398036 +
     t*(-6962890.5431 +
	t*(7.4722 +
	   t*(0.007702 +
	      t*(-0.00005939 )))) % TURNAS) * AS2R
}

/// Mean longitude of Venus.
/// Source: fave03.for
pub fn ve03(t:R)->R {
    (3.176146697 + 1021.3285546211 * t) % TWO_PI
}

/// Mean longitude of Earth.
/// Source: fae03.for
pub fn e03(t:R)->R {
    (1.753470314 + 628.3075849991 * t) % TWO_PI
}

/// General accumulated precession in longitude.
pub fn pa03(t:R)->R {
    (0.024381750 + 0.00000538691 * t) * t
}
