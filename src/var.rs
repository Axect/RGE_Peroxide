use crate::constants::*;

pub struct RGE {
    pub t: f64,
    pub lH: f64,
    pub yt: f64,
    pub g1: f64,
    pub g2: f64,
    pub g3: f64,
    pub phi: f64,
    pub G: f64,
}

pub struct Cosmo {
    pub V: f64,
    pub eps: f64,
    pub eta: f64,
    pub A: f64,
}

impl RGE {
    pub fn new(mt: f64) -> RGE {
        let t = 0f64;
        let yt = 0.93690 + 0.00556 * (mt - 173.34) - 0.00003 * (M_H - 125.15) - 0.00042 * (ALPHA_MZ - 0.1184) / 0.0007;
    
        RGE {
            t: 0f64,
            lH: 0.12604 + 0.00206*(M_H - 125.15) - 0.00004*(mt - 173.34),
            yt: (0.93690 + 0.00556*(mt - 173.34) - 0.00003*(M_H - 125.15) - 0.00042*(ALPHA_MZ - 0.1184) / 0.0007),
            g1: 0.35830 + 0.00011*(mt - 173.34) - 0.00020*(M_W - 80.384) / 0.014,
            g2: 0.64779 + 0.00004*(mt - 173.34) + 0.00011*(M_W - 80.384) / 0.014,
            g3: 1.1666 + 0.00314*(ALPHA_MZ - 0.1184) / 0.007 - 0.00046*(mt - 173.34),
            phi: 2f64.sqrt() / yt * mt * t.exp(),
            G: 1f64,
        }
    }
}






