use crate::constants::*;
use crate::var::*;
use std::f64::consts::PI;

#[derive(Debug, Clone)]
pub struct Beta {
    /// 1-loop order
    //b1_lH: f64,
    //b1_yt: f64,
    //b1_g1: f64,
    //b1_g2: f64,
    //b1_g3: f64,
    //gamma1: f64,
    /// 2-loop order
    //b2_lH: f64,
    //b2_yt: f64,
    //b2_g1: f64,
    //b2_g2: f64,
    //b2_g3: f64,
    //gamma2: f64,
    /// Total
    b_lH: f64,
    b_yt: f64,
    b_g1: f64,
    b_g2: f64,
    b_g3: f64,
    gamma: f64,
}

impl Beta {
    pub fn new(R: &RGE, mt: f64, xi: f64) -> Beta {
        let hg = 2f64.sqrt() / R.yt * mt * R.t.exp();
	let sh = (1. + xi * hg.powi(2)/M_PLANK_R.powi(2)) / (1. + (1. + 6. * xi) * xi * hg.powi(2)/M_PLANK_R.powi(2));

	// 1-loop Beta Function
	let b1lH: f64 = 6. * (1. + 3. * sh.powi(2)) * R.lH.powi(2) + 12. * R.lH * R.yt.powi(2) - 6. * R.yt.powi(4) - 3. * R.lH * (3. * R.g2.powi(2) + R.g1.powi(2)) + 3./8. * (2. * R.g2.powi(4) + (R.g1.powi(2) + R.g2.powi(2)).powi(2));
	let b1yt: f64 = R.yt * ((23./6. + 2./3. * sh) * R.yt.powi(2) - (8. * R.g3.powi(2) + 9./4. * R.g2.powi(2) + 17./12. * R.g1.powi(2)));
	let b1g1: f64 = (81. + sh) / 12. * R.g1.powi(3);
	let b1g2: f64 = (sh - 39.) / 12. * R.g2.powi(3);
	let b1g3: f64 = -7. * R.g3.powi(3);
	let gamma1: f64 = -(9./4. * R.g2.powi(2) + 3./4. * R.g1.powi(2) - 3. * R.yt.powi(2));

        let gamma: f64 = 1. / (16. * PI.powi(2)) * gamma1;

        let g = make_beta(gamma);

        Beta {
            b_lH: g(b1lH, 0f64),
            b_yt: g(b1yt, 0f64),
            b_g1: g(b1g1, 0f64),
            b_g2: g(b1g2, 0f64),
            b_g3: g(b1g3, 0f64),
            gamma: 1f64,
        }
    }
}

fn make_beta(g: f64) -> Box<Fn(f64, f64) -> f64> {
    Box::new(move |x: f64, y: f64| {
        let temp: f64 = 1. / (16. * PI.powi(2)) * x + 1. / (16. * PI.powi(2)).powi(2) * y;
        temp / (1. + g)
    })
}
