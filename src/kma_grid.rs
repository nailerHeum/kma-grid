use std::f64;

const GRID_LENGTH: f64 = 5.0; // km
const EARTH_RADIUS_IN_GRID: f64 = 6371.00877 / GRID_LENGTH; // km

const DEGREE_TO_RADIAN: f64 = 0.017453293; // pi / 180
const RADIAN_TO_DEGREE: f64 = 57.29577951; // 180 / pi

const STANDARD_PARALLEL1: f64 = 30.0 * DEGREE_TO_RADIAN; // radian
const STANDARD_PARALLEL2: f64 = 60.0 * DEGREE_TO_RADIAN; // radian

const REFERENCE_LONGITUDE: f64 = 126.0 * DEGREE_TO_RADIAN; // radian
const REFERENCE_LATITUDE: f64 = 38.0 * DEGREE_TO_RADIAN; // radian

const REFERENCE_X: u8 = 43; // x coordinate of the reference grid
const REFERENCE_Y: u8 = 136; // y coordinate of the reference grid

pub struct KmaGrid {
    x: u8, // x coordinate of the grid, 0 ~ 149
    y: u8, // y coordinate of the grid, 0 ~ 253
}

struct LccConstants {
    n: f64,
    f: f64,
    rho_zero: f64,
}

impl LccConstants {
    fn get_constants() -> Self {
        let n = (STANDARD_PARALLEL1.cos() / STANDARD_PARALLEL2.cos()).ln()
            / (f64::tan(f64::consts::PI * 0.25 + 0.5 * STANDARD_PARALLEL2)
                / f64::tan(f64::consts::PI * 0.25 + 0.5 * STANDARD_PARALLEL1))
            .ln();

        let f = (f64::consts::PI * 0.25 + 0.5 * STANDARD_PARALLEL1)
            .tan()
            .powf(n)
            * STANDARD_PARALLEL1.cos()
            / n;

        let rho_zero = EARTH_RADIUS_IN_GRID * f
            / (f64::consts::PI * 0.25 + 0.5 * REFERENCE_LATITUDE)
                .tan()
                .powf(n);
        LccConstants { n, f, rho_zero }
    }
}

impl KmaGrid {
    // https://en.wikipedia.org/wiki/Lambert_conformal_conic_projection#Transformation
    // Lambert conformal conic projection

    pub fn from_gcs(longitude: f64, latitude: f64) -> KmaGrid {
        let LccConstants { n, f, rho_zero } = LccConstants::get_constants();
        let rho = EARTH_RADIUS_IN_GRID * f
            / (f64::consts::PI * 0.25 + 0.5 * latitude * DEGREE_TO_RADIAN)
                .tan()
                .powf(n);

        let theta = {
            let raw_theta = longitude * DEGREE_TO_RADIAN - REFERENCE_LONGITUDE;
            if raw_theta > f64::consts::PI {
                raw_theta - 2.0 * f64::consts::PI
            } else if raw_theta < -f64::consts::PI {
                raw_theta + 2.0 * f64::consts::PI
            } else {
                raw_theta
            }
        };

        let x = (rho * (theta * n).sin() + (REFERENCE_X as f64)).round() as u8;
        let y = (rho_zero - rho * (theta * n).cos() + (REFERENCE_Y as f64)).round() as u8;

        return KmaGrid { x, y };
    }

    pub fn to_gcs(self: Self) -> (f64, f64) {
        let LccConstants { n, f, rho_zero } = LccConstants::get_constants();
        let xn = (self.x - REFERENCE_X) as f64;
        let yn = rho_zero - (self.y - REFERENCE_Y) as f64;

        let ra = (xn.powi(2) + yn.powi(2)).sqrt();

        let latitude_radian =
            2.0 * (EARTH_RADIUS_IN_GRID * f / ra).powf(1.0 / n).atan() - f64::consts::PI * 0.5;

        let theta = if n.abs() == 0.0 {
            0.0
        } else if yn.abs() == 0.0 {
            f64::consts::PI * 0.5
        } else {
            yn.atan2(xn)
        };

        let longitude_radian = theta / n + REFERENCE_LONGITUDE;

        return (
            longitude_radian * RADIAN_TO_DEGREE,
            latitude_radian * RADIAN_TO_DEGREE,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    mod convert_to_xy {
        use super::*;

        #[test]
        fn assert_almost_eq() {
            let grid = KmaGrid::from_gcs(126.0, 38.0);
            assert_eq!(grid.x, REFERENCE_X);
            assert_eq!(grid.y, REFERENCE_Y);
        }
    }
}
