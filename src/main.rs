extern crate rge_peroxide;

use rge_peroxide::*;

fn main() {
    let r: RGE = RGE::new(175.);
    let a: Beta = Beta::new(&r, 175., 10.);
    println!("{:?}", a);
}
