//! Functions to generate point clouds

use std::{f64::consts::PI, ops::Range};

use rand::Rng;

/// Helper: add i.i.d. uniform noise in `noise_range` to every coordinate in-place
fn add_noise(points: &mut [Vec<f64>], noise_range: Option<Range<f64>>) {
    if let Some(r) = noise_range {
        let mut rng = rand::thread_rng();
        for p in points.iter_mut() {
            for coord in p.iter_mut() {
                *coord += rng.gen_range(r.clone());
            }
        }
    }
}

/// Return `m` points evenly spaced on the unit circle (optionally with uniform noise per coord)
pub fn unit_circle(m: usize, noise_range: Option<Range<f64>>) -> Vec<Vec<f64>> {
    let circpoint = |k: usize| {
        let theta = k as f64 * 2.0 * PI / m as f64;
        vec![theta.cos(), theta.sin()]
    };
    let mut pts: Vec<_> = (0..m).map(circpoint).collect();
    add_noise(&mut pts, noise_range);
    pts
}

/// Return `m` points uniformly distributed on the unit 2-sphere surface in R^3
///
/// Sampling: theta ~ U[0, 2π), z ~ U[-1, 1], x = √(1-z²) cos theta, y = √(1-z²) sin theta.
pub fn unit_sphere_surface(m: usize, noise_range: Option<Range<f64>>) -> Vec<Vec<f64>> {
    let mut rng = rand::thread_rng();
    let mut pts = Vec::with_capacity(m);
    for _ in 0..m {
        let theta = rng.gen_range(0.0..(2.0 * PI));
        let z = rng.gen_range(-1.0..1.0);
        let r = (1.0 - z * z).sqrt();
        let x = r * theta.cos();
        let y = r * theta.sin();
        pts.push(vec![x, y, z]);
    }
    add_noise(&mut pts, noise_range);
    pts
}

/// Return `m` points uniformly distributed *inside* a unit 3D ball (not just the surface)
///
/// Sampling: direction from unit sphere; radius = U[0,1]^(1/3) for volume-uniformity.
pub fn unit_ball(m: usize, noise_range: Option<Range<f64>>) -> Vec<Vec<f64>> {
    let mut rng = rand::thread_rng();
    let mut pts = Vec::with_capacity(m);
    for _ in 0..m {
        // direction
        let theta = rng.gen_range(0.0..(2.0 * PI));
        let z = rng.gen_range(-1.0..1.0);
        let r_xy = (1.0 - z * z).sqrt();
        let dir = [r_xy * theta.cos(), r_xy * theta.sin(), z];

        // radius for uniform volume
        let r = rng.gen::<f64>().cbrt();

        pts.push(vec![r * dir[0], r * dir[1], r * dir[2]]);
    }
    add_noise(&mut pts, noise_range);
    pts
}

/// Return `m` points on a torus surface embedded in R^3
///
/// Parametrization: u,v ~ U[0,2π)
/// (x, y, z) = ((R + r cos v) cos u, (R + r cos v) sin u, r sin v)
/// `major` = R (distance from center to tube center), `minor` = r (tube radius)
pub fn torus(m: usize, major: f64, minor: f64, noise_range: Option<Range<f64>>) -> Vec<Vec<f64>> {
    assert!(major > 0.0 && minor > 0.0, "major and minor radii must be positive");
    let mut rng = rand::thread_rng();
    let mut pts = Vec::with_capacity(m);
    for _ in 0..m {
        let u = rng.gen_range(0.0..(2.0 * PI));
        let v = rng.gen_range(0.0..(2.0 * PI));
        let cx = (major + minor * v.cos()) * u.cos();
        let cy = (major + minor * v.cos()) * u.sin();
        let cz = minor * v.sin();
        pts.push(vec![cx, cy, cz]);
    }
    add_noise(&mut pts, noise_range);
    pts
}