use std::f64;
use std::f64::consts;

const PI2               : f64 = 2.0 * consts::PI;
const PI_DBY_180        : f64 = consts::PI / 180.0;
const D180_DBY_PI       : f64 = 180.0 / consts::PI;

pub const VNC_DEF_IT_LIMIT  : i32 = 2000;
pub const VNC_DEF_EPSILON   : f64 = 1E-12;

const NLM_A             : f64 = 1.0;
const NLM_B             : f64 = 0.5;
const NLM_R             : f64 = 0.5;
const NLM_Q             : f64 = 0.5;
const NLM_G             : f64 = 2.0;

const NLM_DEF_IT_LIMIT  : i32 = 1200;
const NLM_DEF_PREC_THRLD: f64 = 1E-8;


pub struct EllipsoidDescriptor {
    mjsa_m: f64,
    ifltn: f64,
}

pub struct Ellipsoid {
    mjsa_m: f64,  /* major semi-axis */
    ifltn: f64,   /* inverse flattening */
    
    fltn: f64,    /* flattening */  
    mnsa_m: f64,  /* minor semi-axis */      
    ecnt: f64,    /* eccentricity */
    ecnt_sq: f64, /* eccentricity squared */   
} 

impl Ellipsoid {    
    pub fn new(mjsa_m: f64, ifltn: f64) -> Ellipsoid {
        
        if mjsa_m <= 0.0 {
            panic!("Specified major semi-axis 'mjsa_m' should be greater than zero");
        }
        if ifltn <= 0.0 {
            panic!("Specifed inverse flattening 'ifltn' should be greater than zero")
        }        

        let fltn = 1.0 / ifltn;
        let mnsa_m = mjsa_m * (1.0 - fltn);
        let ecnt = (mjsa_m.powi(2) - mnsa_m.powi(2)) / mjsa_m.powi(2);
    
        Ellipsoid {
            mjsa_m,
            ifltn,
            fltn,
            mnsa_m,
            ecnt,
            ecnt_sq: ecnt.powi(2),
        }
    }
    pub fn from_descriptor(ed: &EllipsoidDescriptor) -> Ellipsoid {
        Ellipsoid::new(ed.mjsa_m, ed.ifltn)        
    }
}

pub const WGS72_ELLIPSOID_DESCRIPTOR: EllipsoidDescriptor = EllipsoidDescriptor{ mjsa_m: 6378135.0, ifltn: 298.26 };
pub const WGS84_ELLIPSOID_DESCRIPTOR: EllipsoidDescriptor = EllipsoidDescriptor{ mjsa_m: 6378137.0, ifltn: 298.257223563 };
pub const GRS80_ELLIPSOID_DESCRIPTOR: EllipsoidDescriptor = EllipsoidDescriptor{ mjsa_m: 6378137.0, ifltn: 298.257222100882711 };
pub const PZ90_ELLIPSOID_DESCRIPTOR: EllipsoidDescriptor = EllipsoidDescriptor{ mjsa_m: 6378136.0, ifltn: 298.257839303 };
pub const IERS_ELLIPSOID_DESCRIPTOR: EllipsoidDescriptor = EllipsoidDescriptor{ mjsa_m: 6378136.49, ifltn: 298.25645 };
pub const KRSY_ELLIPSOID_DESCRIPTOR: EllipsoidDescriptor = EllipsoidDescriptor{ mjsa_m: 6378245.0, ifltn: 298.3 };


pub fn wrap(value: f64, bound: f64) -> f64 {
    if bound <= 0.0 {
        panic!("Specified 'bound' value should be greater than zero");
    }

    let mut vl = value.abs();
    let sign = value.signum();

    while vl > bound {
        vl -= bound;
    }

    (vl * sign)
}

pub fn wrap_2pi(value: f64) -> f64 {
    wrap(value, PI2)
}

pub fn lat_1deg_length(lat_rad: f64, el: &Ellipsoid) -> f64 {
    if el.mjsa_m <= 0.0 {
        panic!("Specified ellipsoid's major semiaxis (mjsa_m) should be greater than zero");
    }
    if (el.ecnt_sq < 0.0) || (el.ecnt_sq >= 1.0) {
        panic!("Specified ellipsoid's eccentrisity squared (ecnt_sq) should be in the range (0.0 .. 1.0) exclusively");
    }

    (PI_DBY_180 * el.mjsa_m * (1.0 - el.ecnt_sq) / (1.0 - el.ecnt_sq * lat_rad.sin().powi(2)).powf(1.5)).abs()
}
           
pub fn lon_1deg_length(lat_rad: f64, el: &Ellipsoid) -> f64 {
    if el.mjsa_m <= 0.0 {
        panic!("Specified ellipsoid's major semiaxis (mjsa_m) should be greater than zero");
    }
    if (el.ecnt_sq <= 0.0) || (el.ecnt_sq >= 1.0) {
        panic!("Specified ellipsoid's eccentrisity squared (ecnt_sq) should be in the range (0.0 .. 1.0) exclusively");
    }    

    (PI_DBY_180 * el.mjsa_m * lat_rad.cos() / (1.0 - el.ecnt_sq * lat_rad.sin().powi(2)).sqrt()).abs()
} 


pub fn geopoint_offset_by_deltas(lat_rad: f64, lon_rad: f64, lat_offset_m: f64, lon_offset_m: f64, el: &Ellipsoid) -> (f64, f64) {
    let m_per_deg_lat = lat_1deg_length(lat_rad, el);
    let m_per_deg_lon = lon_1deg_length(lat_rad, el);
    (lat_rad - PI_DBY_180 * lat_offset_m / m_per_deg_lat, lon_rad - PI_DBY_180 * lon_offset_m / m_per_deg_lon)    
}

pub fn geopoint_offset_by_deltas_wgs84(lat_rad: f64, lon_rad: f64, lat_offset_m: f64, lon_offset_m: f64) -> (f64, f64) {
    let m_per_deg_lat = 111132.92 - 559.82 * (2.0 * lat_rad).cos() + 1.175 * (4.0 * lat_rad).cos();
    let m_per_deg_lon = 111412.84 * (lat_rad).cos() - 93.5 * (3.0 * lat_rad).cos();
    (lat_rad - PI_DBY_180 * lat_offset_m / m_per_deg_lat, lon_rad - PI_DBY_180 * lon_offset_m / m_per_deg_lon)    
}

pub fn get_deltas_by_geopoints(sp_lat_rad: f64, sp_lon_rad: f64, ep_lat_rad: f64, ep_lon_rad: f64, el: &Ellipsoid) -> (f64, f64) {
    let m_lat_rad = (sp_lat_rad + ep_lat_rad) / 2.0;
    let m_per_deg_lat = lat_1deg_length(m_lat_rad, el);
    let m_per_deg_lon = lon_1deg_length(m_lat_rad, el);
    ((sp_lat_rad - ep_lat_rad) * m_per_deg_lat * D180_DBY_PI, (sp_lon_rad - ep_lon_rad) * m_per_deg_lon * D180_DBY_PI)
}

pub fn get_deltas_by_geopoints_wgs84(sp_lat_rad: f64, sp_lon_rad: f64, ep_lat_rad: f64, ep_lon_rad: f64) -> (f64, f64) {    
    let m_lat_rad = (sp_lat_rad + ep_lat_rad) / 2.0;
    let m_per_deg_lat = 111132.92 - 559.82 * (2.0 * m_lat_rad).cos() + 1.175 * (4.0 * m_lat_rad).cos();
    let m_per_deg_lon = 111412.84 * (m_lat_rad).cos() - 93.5 * (3.0 * m_lat_rad).cos();
    ((sp_lat_rad - ep_lat_rad) * m_per_deg_lat * D180_DBY_PI, (sp_lon_rad - ep_lon_rad) * m_per_deg_lon * D180_DBY_PI)    
}
 

pub fn haversine_inverse(sp_lat_rad: f64, sp_lon_rad: f64, ep_lat_rad: f64, ep_lon_rad: f64, e_radius_m: f64) -> f64 {
    let a = ((ep_lat_rad - sp_lat_rad) / 2.0).sin().powi(2) + sp_lat_rad.cos() * ep_lat_rad.cos() *
            ((ep_lon_rad - sp_lon_rad) / 2.0).sin().powi(2);
    (e_radius_m * 2.0 * a.sqrt().atan2((1.0 - a).sqrt()))
}
     
pub fn haversine_direct(sp_lat_rad: f64, sp_lon_rad: f64, dst_m: f64, fwd_az_rad: f64, e_radius_m: f64) -> (f64, f64) {
    let delta = dst_m / e_radius_m;
    let ep_lat_rad = wrap_2pi((sp_lat_rad.sin() * delta.cos() + sp_lat_rad.cos() * delta.sin() * fwd_az_rad.cos()).asin());
    let ep_lon_rad = wrap_2pi(PI2 + consts::PI + 
        (sp_lon_rad + (fwd_az_rad.sin() * delta.sin() * sp_lat_rad.cos()).atan2(delta.cos() - sp_lat_rad.sin().powi(2)))) - consts::PI;
    (ep_lat_rad, ep_lon_rad)
}
        
pub fn haversine_initial_bearing(sp_lat_rad: f64, sp_lon_rad: f64, ep_lat_rad: f64, ep_lon_rad: f64) -> f64 {
    let y = (ep_lon_rad - sp_lon_rad).sin() * ep_lat_rad.cos();
    let x = sp_lat_rad.cos() * ep_lat_rad.sin() - sp_lat_rad.sin() * ep_lat_rad.cos() * (ep_lon_rad - sp_lon_rad).cos();
    wrap_2pi(y.atan2(x))
}
        
pub fn haversine_final_bearing(sp_lat_rad: f64, sp_lon_rad: f64, ep_lat_rad: f64, ep_lon_rad: f64) -> f64 {
    wrap_2pi(consts::PI + haversine_initial_bearing(ep_lat_rad, ep_lon_rad, sp_lat_rad, sp_lon_rad))
} 


pub fn vincenty_inverse(sp_lat_rad: f64, sp_lon_rad: f64, ep_lat_rad: f64, ep_lon_rad: f64, el: &Ellipsoid, eps: f64, it_limit: i32) -> (f64, f64, f64, i32, bool) {
    
    let l_ = ep_lon_rad - sp_lon_rad;
    let tan_u_1 = (1.0 - el.fltn) * sp_lat_rad.tan();
    let cos_u_1 = 1.0 / (1.0 + tan_u_1.powi(2)).sqrt();
    let sin_u_1 = tan_u_1 * cos_u_1;

    let tan_u_2 = (1.0 - el.fltn) * ep_lat_rad.tan();
    let cos_u_2 = 1.0 / (1.0 + tan_u_2 * tan_u_2).sqrt();
    let sin_u_2 = tan_u_2 * cos_u_2;
            
    let mut sin_lambda;
    let mut cos_lambda;
    let mut sin_sigma = 0.0;
    let mut cos_sigma = 0.0;
    let mut sin_alpha;
    let mut sin_sq_sigma;
    let mut cos_sq_alpha = 0.0;
    let mut cos_2_sigma_m = 0.0;
    let mut sigma = 0.0;
    let mut c_;

    let mut lambda = l_;
    let mut lambda_ = 0.0;
    let mut its = 0;

    let mut it_check = 0.0;    
    let antimeridian : bool = l_.abs() > consts::PI;

    while {

        sin_lambda = lambda.sin();
        cos_lambda = lambda.cos();

        sin_sq_sigma = (cos_u_2 * sin_lambda) * (cos_u_2 * sin_lambda) +
                       (cos_u_1 * sin_u_2 - sin_u_1 * cos_u_2 * cos_lambda) * (cos_u_1 * sin_u_2 - sin_u_1 * cos_u_2 * cos_lambda);

        if sin_sq_sigma.abs() > f64::EPSILON
        {
            sin_sigma = sin_sq_sigma.sqrt();
            cos_sigma = sin_u_1 * sin_u_2 + cos_u_1 * cos_u_2 * cos_lambda;
            sigma = sin_sigma.atan2(cos_sigma);
            sin_alpha = cos_u_1 * cos_u_2 * sin_lambda / sin_sigma;

            cos_sq_alpha = 1.0 - sin_alpha * sin_alpha;
            
            if cos_sq_alpha != 0.0 {
                cos_2_sigma_m = cos_sigma - 2.0 * sin_u_1 * sin_u_2 / cos_sq_alpha;
            }
            else {
                cos_sq_alpha = 0.0;
            }

            c_ = el.fltn / 16.0 * cos_sq_alpha * (4.0 + el.fltn * (4.0 - 3.0 * cos_sq_alpha));
            lambda_ = lambda;
            lambda = l_ + (1.0 - c_) * el.fltn * sin_alpha *
                     (sigma + c_ * sin_sigma * (cos_2_sigma_m + c_ * cos_sigma * (-1.0 + 2.0 * cos_2_sigma_m * cos_2_sigma_m)));
        
            if antimeridian {
                it_check = lambda.abs() - consts::PI;
            }
            else
            {
                it_check = lambda.abs();
            }        
        }
            
        its += 1;
        (((lambda - lambda_).abs() > eps) && (its < it_limit) && (it_check < consts::PI))
    } { }


    let u_sq = cos_sq_alpha * (el.mjsa_m.powi(2) - el.mnsa_m.powi(2)) / el.mnsa_m.powi(2);
    let a_ = 1.0 + u_sq/16384.0 * (4096.0 + u_sq * (-768.0 + u_sq * (320.0 - 175.0 * u_sq)));
    let b_ = u_sq / 1024.0 * (256.0 + u_sq * (-128.0 + u_sq * (74.0 - 47.0 * u_sq)));
    let delta_sigma = b_ * sin_sigma * (cos_2_sigma_m + b_/4.0 * (cos_sigma * (-1.0 + 2.0 * cos_2_sigma_m * cos_2_sigma_m) -
        b_/6.0 * cos_2_sigma_m * (-3.0 + 4.0 * sin_sigma * sin_sigma) * (-3.0 + 4.0 * cos_2_sigma_m * cos_2_sigma_m)));

    let dst_m = el.mnsa_m * a_ * (sigma - delta_sigma);

    let fwd_az_rad = (cos_u_2 * sin_lambda).atan2(cos_u_1 * sin_u_2 - sin_u_1 * cos_u_2 * cos_lambda);
    let rev_az_rad = (cos_u_1 * sin_lambda).atan2(- sin_u_1 * cos_u_2 + cos_u_1 * sin_u_2 * cos_lambda);

    (dst_m, wrap_2pi(fwd_az_rad), wrap_2pi(rev_az_rad), its, ((its < it_limit) && (it_check < consts::PI)))
}

pub fn vincenty_direct(sp_lat_rad: f64, sp_lon_rad: f64, fwd_az_rad: f64, dst_m: f64, el: &Ellipsoid, eps: f64, it_limit: i32) -> (f64, f64, f64, i32) {

    let sin_alpha_1 = fwd_az_rad.sin();
    let cos_alpha_1 = fwd_az_rad.cos();
    let tan_u_1 = (1.0 - el.fltn) * sp_lat_rad.tan();
    let cos_u_1 = 1.0 / (1.0 + tan_u_1 * tan_u_1).sqrt();
    let sin_u_1 = tan_u_1 * cos_u_1;

    let sigma_1 = tan_u_1.atan2(cos_alpha_1);
    let sin_alpha = cos_u_1 * sin_alpha_1;
    let cos_sq_alpha = 1.0 - sin_alpha * sin_alpha;
    let u_sq = cos_sq_alpha * (el.mjsa_m.powi(2) - el.mnsa_m.powi(2)) / el.mnsa_m.powi(2);
    let a_ = 1.0 + u_sq / 16384.0 * (4096.0 + u_sq * (-768.0 + u_sq * (320.0 - 175.0 * u_sq)));
    let b_ = u_sq / 1024.0 * (256.0 + u_sq * (-128.0 + u_sq * (74.0 - 47.0 * u_sq)));

    let mut cos_2_sigma_m;
    let mut sin_sigma;
    let mut cos_sigma;
    let mut delta_sigma;

    let mut sigma = dst_m / (el.mnsa_m * a_);
    let mut sigma_;
    let mut its = 0;
            
    while {
        cos_2_sigma_m = (2.0 * sigma_1 + sigma).cos();
        sin_sigma = sigma.sin();
        cos_sigma = sigma.cos();

        delta_sigma = b_ * sin_sigma * (cos_2_sigma_m + b_ / 4.0 * (cos_sigma * (-1.0 + 2.0 * cos_2_sigma_m * cos_2_sigma_m) -
                      b_ / 6.0 * cos_2_sigma_m * (-3.0 + 4.0 * sin_sigma * sin_sigma) * (-3.0 + 4.0 * cos_2_sigma_m * cos_2_sigma_m)));

        sigma_ = sigma;
        sigma = dst_m / (el.mnsa_m * a_) + delta_sigma;

        its += 1;
        (((sigma - sigma_).abs() > eps) && (its < it_limit))
    } { }
    
    let x = sin_u_1 * sin_sigma - cos_u_1 * cos_sigma * cos_alpha_1;
    let ep_lat_rad = (sin_u_1 * cos_sigma + cos_u_1 * sin_sigma * cos_alpha_1).atan2((1.0 - el.fltn) * (sin_alpha * sin_alpha + x * x).sqrt());
    
    let lambda = (sin_sigma * sin_alpha_1).atan2(cos_u_1 * cos_sigma - sin_u_1 * sin_sigma * cos_alpha_1);
    let c_ = el.fltn / 16.0 * cos_sq_alpha * (4.0 + el.fltn * (4.0 - 3.0 * cos_sq_alpha));
    
    let l_ = lambda - (1.0 - c_) * el.fltn * sin_alpha * (sigma + c_ * sin_sigma * 
             (cos_2_sigma_m + c_ * cos_sigma * (-1.0 + 2.0 * cos_2_sigma_m * cos_2_sigma_m)));
    
    let ep_lon_rad = sp_lon_rad + l_;
    let rev_az_rad = sin_alpha.atan2(-x);
    (wrap_2pi(ep_lat_rad), wrap_2pi(ep_lon_rad), wrap_2pi(rev_az_rad), its)
}



pub struct TOABasePoint {
    x: f64,
    y: f64,
    z: f64,
    d: f64,
}   

pub struct TDOABaseline {
    x1: f64,
    y1: f64,
    z1: f64,
    x2: f64,
    y2: f64,
    z2: f64,
    prd: f64, // Pseudorange difference    
}

pub type Eps3dFunc<T> = fn(&Vec<T>, f64, f64, f64) -> f64;

pub fn eps_toa3d(base_points: &Vec<TOABasePoint>, x: f64, y: f64, z: f64) -> f64 {
    
    let mut result: f64 = 0.0;
    
    for base_point in base_points {
        result += (((base_point.x - x).powi(2) +
                    (base_point.y - y).powi(2) +
                    (base_point.z - z).powi(2)).sqrt() - base_point.d).powi(2);
    }

    result
}
 
pub fn eps_tdoa3d(base_lines: &Vec<TDOABaseline>, x: f64, y: f64, z: f64) -> f64 {

    let mut result: f64 = 0.0;
    
    for base_line in base_lines {
        result += (((base_line.x1 - x).powi(2) +
                    (base_line.y1 - y).powi(2) +
                    (base_line.z1 - z).powi(2)).sqrt() -                  
                   ((base_line.x2 - x).powi(2) +
                    (base_line.y2 - y).powi(2) +
                    (base_line.z2 - z).powi(2)).sqrt() - base_line.prd).powi(2);        
    }

    result
}


pub fn nlm_2d_solve<T>(eps: Eps3dFunc<T>, base_elements: &Vec<T>, x_prev: f64, y_prev: f64, z: f64,
                  max_iterations: i32, precision_threshold: f64, simplex_size: f64) -> (f64, f64, f64, i32) {

    let mut is_finished: bool = false;
                  
    let mut tmp: f64;
    let mut tmp1: f64;
    let mut xix: [f64; 3] = [0.0; 3];
    let mut xiy: [f64; 3] = [0.0; 3];
    let mut fxi: [f64; 3] = [0.0; 3];
    
    let mut fl;
    let mut fg;
    let mut fh;
    let mut fr; 
    let mut fe;
    let mut fs;
    let mut xcx;
    let mut xcy;
    let mut xrx;
    let mut xry;
    let mut xex; 
    let mut xey;
    let mut xsx;
    let mut xsy;
      
    let mut it_cnt: i32 = 0;
      
    xix[0] = x_prev;
    xiy[0] = y_prev;
    xix[1] = xix[0] + simplex_size;
    xiy[1] = xiy[0] + simplex_size;
    xix[2] = xix[0] - simplex_size / 2.0;
    xiy[2] = xiy[0] + simplex_size / 2.0;

    while !is_finished {

        fxi[0] = eps(base_elements, xix[0], xiy[0], z);
        fxi[1] = eps(base_elements, xix[1], xiy[1], z);
        fxi[2] = eps(base_elements, xix[2], xiy[2], z);

        if fxi[0] > fxi[1] {
            tmp = fxi[0]; fxi[0] = fxi[1]; fxi[1] = tmp;
            tmp = xix[0]; xix[0] = xix[1]; xix[1] = tmp;
            tmp = xiy[0]; xiy[0] = xiy[1]; xiy[1] = tmp;
        }

        if fxi[0] > fxi[2] {
            tmp = fxi[0]; fxi[0] = fxi[2]; fxi[2] = tmp;
            tmp = xix[0]; xix[0] = xix[2]; xix[2] = tmp;
            tmp = xiy[0]; xiy[0] = xiy[2]; xiy[2] = tmp;
        }

        if fxi[1] > fxi[2] {
            tmp = fxi[1]; fxi[1] = fxi[2]; fxi[2] = tmp;
            tmp = xix[1]; xix[1] = xix[2]; xix[2] = tmp;
            tmp = xiy[1]; xiy[1] = xiy[2]; xiy[2] = tmp;
        }

        fl = fxi[0];
        fg = fxi[1];
        fh = fxi[2];

        xcx = (xix[0] + xix[1]) / 2.0;
        xcy = (xiy[0] + xiy[1]) / 2.0;

        xrx = (1.0 + NLM_A) * xcx - NLM_A * xix[2];
        xry = (1.0 + NLM_A) * xcy - NLM_A * xiy[2];

        fr = eps(base_elements, xrx, xry, z);

        if fr < fl {
            xex = (1.0 - NLM_G) * xcx + NLM_G * xrx;
            xey = (1.0 - NLM_G) * xcy + NLM_G * xry;

            fe = eps(base_elements, xex, xey, z);

            if fe < fr {
                xix[2] = xex;
                xiy[2] = xey;
            }
            else {
                xix[2] = xrx;
                xiy[2] = xry;
            }
        }
        else
        {
            if (fr > fl) && (fr < fg) {
                xix[2] = xrx;
                xiy[2] = xry;
            }
            else {
                if (fr > fg) && (fr < fh) {
                    xix[2] = xrx; //tmp = xix[2]; xix[2] = xrx; xrx = tmp;
                    xiy[2] = xry; //tmp = xiy[2]; xiy[2] = xry; xry = tmp;
                    fxi[2] = fr;  //tmp = fxi[2]; fxi[2] = fr; fr = tmp;
                }
                else {
                    if fh < fr {
                        //
                    }
                }

                xsx = NLM_B * xix[2] + (1.0 - NLM_B) * xcx;
                xsy = NLM_B * xiy[2] + (1.0 - NLM_B) * xcy;
                fs = eps(base_elements, xsx, xsy, z);

                if fs < fh {
                    xix[2] = xsx;
                    xiy[2] = xsy;
                }
                else {
                    xix[1] = (xix[1] - xix[0]) / 2.0;
                    xiy[1] = (xiy[1] - xiy[0]) / 2.0;
                    xix[2] = (xix[2] - xix[0]) / 2.0;
                    xiy[2] = (xiy[2] - xiy[0]) / 2.0;
                }
            }
        }

        tmp = (fxi[0] + fxi[1] + fxi[2]) / 3.0;
        tmp1 = ((fxi[0] - tmp) * (fxi[0] - tmp) +
                (fxi[1] - tmp) * (fxi[1] - tmp) +
                (fxi[2] - tmp) * (fxi[2] - tmp)) / 3.0;

        it_cnt += 1;
        is_finished = (it_cnt < max_iterations) && (tmp1.sqrt() <= precision_threshold);
    }

    //
    (xix[0], xiy[0], eps(base_elements, xix[0], xiy[0], z).sqrt(), it_cnt)
}

pub fn nlm_3d_solve<T>(eps: Eps3dFunc<T>, base_elements: &Vec<T>, x_prev: f64, y_prev: f64, z_prev: f64,
    max_iterations: i32, precision_threshold: f64, simplex_size: f64) -> (f64, f64, f64, f64, i32) {    
   
    let mut xix: [f64; 4] = [0.0; 4];
    let mut xiy: [f64; 4] = [0.0; 4];
    let mut xia: [f64; 4] = [0.0; 4];
    let mut fxi: [f64; 4] = [0.0; 4];

    xix[0] = x_prev;
    xiy[0] = y_prev;
    xia[0] = z_prev;
    xix[1] = xix[0] + simplex_size;
    xiy[1] = xiy[0] + simplex_size;
    xia[1] = xia[0] - simplex_size;
    xix[2] = xix[0] - simplex_size / 2.0;
    xiy[2] = xiy[0] + simplex_size / 2.0;
    xia[2] = xia[0] + simplex_size / 2.0;
    xix[3] = xix[0] + simplex_size / 2.0;
    xiy[3] = xiy[0] - simplex_size / 2.0;
    xia[3] = xia[0] - simplex_size / 2.0;

    let mut is_finished: bool = false;
    let mut sigma;
    let mut mean;
    let mut tmp;
    let mut xrx;
    let mut xry;
    let mut xra;
    let mut x0x;
    let mut x0y;
    let mut x0a;
    let mut xex;
    let mut xey;
    let mut xea;
    let mut xcx;
    let mut xcy;
    let mut xca;
    let mut fc;
    let mut fr;
    let mut fe;
    let mut it_cnt: i32 = 0;

    while !is_finished {
        
        for i in 0..4 {
            fxi[i] = eps(base_elements, xix[i], xiy[i], xia[i]);
        }

        mean = (fxi[0] + fxi[1] + fxi[2] + fxi[3]) / 4.0;
        sigma = (((fxi[0] - mean).powi(2) +
                  (fxi[1] - mean).powi(2) +
                  (fxi[2] - mean).powi(2) +
                  (fxi[3] - mean).powi(2)) / 4.0).sqrt();

        it_cnt += 1;
        if (it_cnt > max_iterations) || (sigma < precision_threshold) {
           is_finished = true;
        }
        else
        {
            // Sort vertices
            if fxi[0] > fxi[1] {
                tmp = fxi[0]; fxi[0] = fxi[1]; fxi[1] = tmp;
                tmp = xix[0]; xix[0] = xix[1]; xix[1] = tmp;
                tmp = xiy[0]; xiy[0] = xiy[1]; xiy[1] = tmp;
                tmp = xia[0]; xia[0] = xia[1]; xia[1] = tmp;
            }
            if fxi[0] > fxi[2] {
                tmp = fxi[0]; fxi[0] = fxi[2]; fxi[2] = tmp;
                tmp = xix[0]; xix[0] = xix[2]; xix[2] = tmp;
                tmp = xiy[0]; xiy[0] = xiy[2]; xiy[2] = tmp;
                tmp = xia[0]; xia[0] = xia[2]; xia[2] = tmp;
            }
            if fxi[0] > fxi[3] {
                tmp = fxi[0]; fxi[0] = fxi[3]; fxi[3] = tmp;
                tmp = xix[0]; xix[0] = xix[3]; xix[3] = tmp;
                tmp = xiy[0]; xiy[0] = xiy[3]; xiy[3] = tmp;
                tmp = xia[0]; xia[0] = xia[3]; xia[3] = tmp;
            }
            if fxi[1] > fxi[2] {
                tmp = fxi[1]; fxi[1] = fxi[2]; fxi[2] = tmp;
                tmp = xix[1]; xix[1] = xix[2]; xix[2] = tmp;
                tmp = xiy[1]; xiy[1] = xiy[2]; xiy[2] = tmp;
                tmp = xia[1]; xia[1] = xia[2]; xia[2] = tmp;
            }
            if fxi[1] > fxi[3] {
                tmp = fxi[1]; fxi[1] = fxi[3]; fxi[3] = tmp;
                tmp = xix[1]; xix[1] = xix[3]; xix[3] = tmp;
                tmp = xiy[1]; xiy[1] = xiy[3]; xiy[3] = tmp;
                tmp = xia[1]; xia[1] = xia[3]; xia[3] = tmp;
            }
            if fxi[2] > fxi[3] {
                tmp = fxi[2]; fxi[2] = fxi[3]; fxi[3] = tmp;
                tmp = xix[2]; xix[2] = xix[3]; xix[3] = tmp;
                tmp = xiy[2]; xiy[2] = xiy[3]; xiy[3] = tmp;
                tmp = xia[2]; xia[2] = xia[3]; xia[3] = tmp;
            }

            // (2) Calculate x0, the centroid of all points except xn+1
            x0x = (xix[0] + xix[1] + xix[2]) / 3.0;
            x0y = (xiy[0] + xiy[1] + xiy[2]) / 3.0;
            x0a = (xia[0] + xia[1] + xia[2]) / 3.0;

            // (3) reflect the xh (xi(4)) point to xr
            xrx = x0x + NLM_A * (x0x - xix[3]);
            xry = x0y + NLM_A * (x0y - xiy[3]);
            xra = x0a + NLM_A * (x0a - xia[3]);

            // function value in the xr point                    
            fr = eps(base_elements, xrx, xry, xra);

            // if fx1 <= fxr <= fxn replace worst point xn+1 with xr
            if (fr >= fxi[0]) && (fxi[2] >= fr) {
                xix[3] = xrx;
                xiy[3] = xry;
                xia[3] = xra;
                fxi[3] = fr;
            }
            else {
                // (4) expansion
                if fr < fxi[0] {
                    xex = x0x + NLM_G * (xrx - x0x);
                    xey = x0x + NLM_G * (xry - x0x);
                    xea = x0x + NLM_G * (xra - x0x);
                    fe = eps(base_elements, xex, xey, xea);

                    if fe < fr {
                        xix[3] = xex;
                        xiy[3] = xey;
                        xia[3] = xea;
                        fxi[3] = fe;
                        // % to step 1
                    }
                    else {
                        xix[3] = xrx;
                        xiy[3] = xry;
                        xia[3] = xra;
                        fxi[3] = fr;
                        // % to step 1
                    }
                }
                else {
                    xcx = x0x + NLM_R * (xix[3] - x0x);
                    xcy = x0y + NLM_R * (xiy[3] - x0y);
                    xca = x0a + NLM_R * (xia[3] - x0a);
                    fc = eps(base_elements, xcx, xcy, xca);

                    if fc < fxi[3] {
                        xix[3] = xcx;
                        xiy[3] = xcy;
                        xia[3] = xca;
                        fxi[3] = fc;
                        // % to step 1
                    }
                    else {
                        xix[1] = xix[0] + NLM_Q * (xix[1] - xix[0]);
                        xiy[1] = xiy[0] + NLM_Q * (xiy[1] - xiy[0]);
                        xia[1] = xia[0] + NLM_Q * (xia[1] - xia[0]);

                        xix[2] = xix[0] + NLM_Q * (xix[2] - xix[0]);
                        xiy[2] = xiy[0] + NLM_Q * (xiy[2] - xiy[0]);
                        xia[2] = xia[0] + NLM_Q * (xia[2] - xia[0]);

                        xix[3] = xix[0] + NLM_Q * (xix[3] - xix[0]);
                        xiy[3] = xiy[0] + NLM_Q * (xiy[3] - xiy[0]);
                        xia[3] = xia[0] + NLM_Q * (xia[3] - xia[0]);
                    }
                }
            }
        }
    }
        
    (xix[0], xiy[0], xia[0], eps(base_elements, xix[0], xiy[0], xia[0]).sqrt(), it_cnt)
}

                                                            
pub fn toa_nlm_2d_solve(base_points: &Vec<TOABasePoint>, x_prev: f64, y_prev: f64, z: f64,
    max_iterations: i32, precision_threshold: f64, simplex_size: f64) -> (f64, f64, f64, i32) {
    nlm_2d_solve::<TOABasePoint>(eps_toa3d, base_points, x_prev, y_prev, z, max_iterations, precision_threshold, simplex_size)    
}


// TODO: tests for routines below

pub fn tdoa_nlm_2d_solve(base_lines: &Vec<TDOABaseline>,  x_prev: f64, y_prev: f64, z: f64,
    max_iterations: i32, precision_threshold: f64, simplex_size: f64) -> (f64, f64, f64, i32) {
    nlm_2d_solve::<TDOABaseline>(eps_tdoa3d, base_lines, x_prev, y_prev, z, max_iterations, precision_threshold, simplex_size)       
}
       
pub fn toa_nlm_3d_solve(base_points: &Vec<TOABasePoint>, x_prev: f64, y_prev: f64, z_prev: f64,
    max_iterations: i32, precision_threshold: f64, simplex_size: f64) -> (f64, f64, f64, f64, i32) {

    nlm_3d_solve::<TOABasePoint>(eps_toa3d, base_points, x_prev, y_prev, z_prev, max_iterations, precision_threshold, simplex_size)
}

pub fn tdoa_nlm_3d_solve(base_lines: &Vec<TDOABaseline>,  x_prev: f64, y_prev: f64, z_prev: f64,
    max_iterations: i32, precision_threshold: f64, simplex_size: f64) -> (f64, f64, f64, f64, i32) {

    nlm_3d_solve::<TDOABaseline>(eps_tdoa3d, base_lines, x_prev, y_prev, z_prev, max_iterations, precision_threshold, simplex_size)
}

pub fn get_nearest_item_index(base_points: &Vec<TOABasePoint>) -> usize {

    let mut nrst_idx = 0;
    let mut min_dst: f64 = f64::MIN;

    for (idx, base_point) in base_points.iter().enumerate() {
        if base_point.d < min_dst {
            min_dst = base_point.d;
            nrst_idx = idx;
        }
    }

    nrst_idx
}

pub fn toa_circles_intersection_solve(base_points: &Vec<TOABasePoint>, anchor_x: f64, anchor_y: f64, radius: f64, z: f64,
                            arc_mid_rad: f64, arc_angle_rad: f64, steps: i32) -> f64 {

    let mut a: f64 = arc_mid_rad - arc_angle_rad / 2.0;
    let a_end: f64 = arc_mid_rad + arc_angle_rad / 2.0;
    let mut a_best: f64 = a;
    let step_rad: f64 = arc_angle_rad / steps as f64;
    let mut eps_best: f64 = f64::MAX;
    let mut x;
    let mut y;
    let mut eps;

    while a < a_end {
        x = anchor_x + radius * a.cos();
        y = anchor_y + radius * a.sin();
        eps = eps_toa3d(base_points, x, y, z);

        if eps < eps_best {
            eps_best = eps;
            a_best = a;
        }

        a += step_rad;
    }

    a_best
}

pub fn toa_circles_1d_solve(base_points: &Vec<TOABasePoint>, z: f64, end_arc_angle_rad: f64, steps: i32, arc_angle_decrease_factor: f64) -> (f64, f64, f64) {                
    let nrst_idx = get_nearest_item_index(base_points);
    let d_z: f64 = (base_points[nrst_idx].z - z).abs();
    let radius: f64 = if base_points[nrst_idx].d < d_z { 0.0 } else { (base_points[nrst_idx].d.powi(2) - d_z.powi(2)).sqrt() };
    
    let anchor_x: f64 = base_points[nrst_idx].x;
    let anchor_y: f64 = base_points[nrst_idx].y;
    let mut arc_angle: f64 = PI2;
    let mut alpha: f64 = 0.0;    

    while arc_angle > end_arc_angle_rad {    
        alpha = toa_circles_intersection_solve(base_points, anchor_x, anchor_y, radius, z, alpha, arc_angle, steps);
        arc_angle *= arc_angle_decrease_factor;
    }

    let x_best = anchor_x + radius * alpha.cos();
    let y_best = anchor_y + radius * alpha.sin();
    (x_best, y_best, eps_toa3d(base_points, x_best, y_best, z).sqrt())
} 
 
pub fn dist_3d(x1: f64, y1: f64, z1: f64, x2: f64, y2: f64, z2: f64) -> f64 {
    ((x1 - x2).powi(2) + (y1 - y2).powi(2) + (z1 - z2).powi(2)).sqrt()
}


#[macro_export]
macro_rules! assert_approx_eq {
    ($a:expr, $b:expr, $eps:expr) => {{
        let (a, b) = (&$a, &$b);
        let eps = $eps;
        assert!(
            (*a - *b).abs() < eps,
            "assertion failed \
             (l_op: `{:?}`, r_op: `{:?}`, eps: `{:?}`, real diff: `{:?}`)",
            *a,
            *b,
            eps,
            (*a - *b).abs()
        );
    }};
}

#[cfg(test)]
mod tests {
    use super::*;
    extern crate rand;
    use rand::{Rng, thread_rng};

    #[test]
    fn test_const_ellipsoid_descriptors() {
        assert_approx_eq!(WGS72_ELLIPSOID_DESCRIPTOR.mjsa_m, 6378135.0, 1E-6);
        assert_approx_eq!(WGS72_ELLIPSOID_DESCRIPTOR.ifltn, 298.26, 1E-6);
            
        assert_approx_eq!(WGS84_ELLIPSOID_DESCRIPTOR.mjsa_m, 6378137.0, 1E-6);
        assert_approx_eq!(WGS84_ELLIPSOID_DESCRIPTOR.ifltn, 298.257223563, 1E-6);

        assert_approx_eq!(GRS80_ELLIPSOID_DESCRIPTOR.mjsa_m, 6378137.0, 1E-6); 
        assert_approx_eq!(GRS80_ELLIPSOID_DESCRIPTOR.ifltn, 298.257222100882711, 1E-6);

        assert_approx_eq!(PZ90_ELLIPSOID_DESCRIPTOR.mjsa_m, 6378136.0, 1E-6);
        assert_approx_eq!(PZ90_ELLIPSOID_DESCRIPTOR.ifltn, 298.257839303, 1E-6);

        assert_approx_eq!(IERS_ELLIPSOID_DESCRIPTOR.mjsa_m, 6378136.49, 1E-6);
        assert_approx_eq!(IERS_ELLIPSOID_DESCRIPTOR.ifltn, 298.25645, 1E-6);

        assert_approx_eq!(KRSY_ELLIPSOID_DESCRIPTOR.mjsa_m, 6378245.0, 1E-6);
        assert_approx_eq!(KRSY_ELLIPSOID_DESCRIPTOR.ifltn, 298.3, 1E-6);        
    }    

    #[test]
    #[should_panic(expected = "Specified major semi-axis 'mjsa_m' should be greater than zero")]
    fn test_ellipsoid_new_non_positive_major_semiaxis_panic() {
        let _el: Ellipsoid = Ellipsoid::new(-1.0, 0.0);
    }

    #[test]
    #[should_panic(expected = "Specifed inverse flattening 'ifltn' should be greater than zero")]
    fn test_ellipsoid_new_inverse_flattening_out_of_range_panic() {
        let _el: Ellipsoid = Ellipsoid::new(1.0, 0.0);
    }

    #[test]
    #[should_panic(expected = "Specified major semi-axis 'mjsa_m' should be greater than zero")]
    fn test_ellipsoid_from_descriptor_non_positive_major_semiaxis_panic() {
        let ed: EllipsoidDescriptor = EllipsoidDescriptor { mjsa_m: -1.0, ifltn: 0.0 };
        let _el: Ellipsoid = Ellipsoid::from_descriptor(&ed);
    }

    #[test]
    #[should_panic(expected = "Specifed inverse flattening 'ifltn' should be greater than zero")]
    fn test_ellipsoid_from_descriptor_inverse_flattening_out_of_range_panic() {
        let ed: EllipsoidDescriptor = EllipsoidDescriptor { mjsa_m: 1.0, ifltn: 0.0 };
        let _el: Ellipsoid = Ellipsoid::from_descriptor(&ed);
    }

    #[test]
    fn test_ellipsoid_new() {
        let el: Ellipsoid = Ellipsoid::new(WGS84_ELLIPSOID_DESCRIPTOR.mjsa_m, WGS84_ELLIPSOID_DESCRIPTOR.ifltn);        
        assert_approx_eq!(el.mjsa_m, WGS84_ELLIPSOID_DESCRIPTOR.mjsa_m, 1E-6);
        assert_approx_eq!(el.ifltn, WGS84_ELLIPSOID_DESCRIPTOR.ifltn, 1E-6);
        assert_approx_eq!(el.fltn, 1.0 / el.ifltn, 1E-6);
        assert_approx_eq!(el.mnsa_m, el.mjsa_m * (1.0 - el.fltn), 1E-6);
        assert_approx_eq!(el.ecnt, (el.mjsa_m.powi(2) - el.mnsa_m.powi(2)) / el.mjsa_m.powi(2), 1E-6);        
    }

    #[test]
    fn test_ellipsoid_from_descriptor() {
        let el: Ellipsoid = Ellipsoid::from_descriptor(&WGS84_ELLIPSOID_DESCRIPTOR);
        assert_approx_eq!(el.mjsa_m, WGS84_ELLIPSOID_DESCRIPTOR.mjsa_m, 1E-6);
        assert_approx_eq!(el.ifltn, WGS84_ELLIPSOID_DESCRIPTOR.ifltn, 1E-6);
        assert_approx_eq!(el.fltn, 1.0 / el.ifltn, 1E-6);
        assert_approx_eq!(el.mnsa_m, el.mjsa_m * (1.0 - el.fltn), 1E-6);
        assert_approx_eq!(el.ecnt, (el.mjsa_m.powi(2) - el.mnsa_m.powi(2)) / el.mjsa_m.powi(2), 1E-6);
    }

    #[test]
    fn test_wrap() {
        assert_eq!(wrap(10.0, 5.0), 5.0);
        assert_eq!(wrap(-10.0, 5.0), -5.0);
        assert_eq!(wrap(0.0, 5.0), 0.0);
        assert_eq!(wrap(5.0, 10.0), 5.0);
        assert_eq!(wrap(-5.0, 10.0), -5.0);
    }

    #[test]
    #[should_panic(expected = "Specified 'bound' value should be greater than zero")]
    fn test_wrap_non_positive_bound_panic() {
        wrap(0.0, -1.0);
    }

    #[test]
    fn test_wrap_2pi() {
        assert_eq!(wrap_2pi(PI2), PI2);
        assert_eq!(wrap_2pi(consts::PI * -3.0), -consts::PI);
        assert_eq!(wrap_2pi(consts::PI / 3.0), consts::PI / 3.0);
    }

    #[test]
    fn test_lat_1deg_length() {
        let wgs84_ellipsoid: Ellipsoid = Ellipsoid::from_descriptor(&WGS84_ELLIPSOID_DESCRIPTOR);        
        assert_approx_eq!(lat_1deg_length(-3.141592654 / 2.0, &wgs84_ellipsoid), lat_1deg_length(3.141592654 / 2.0, &wgs84_ellipsoid), 1E-6);
        
        assert_approx_eq!(lat_1deg_length(-3.141592654, &wgs84_ellipsoid), 111314.502041079, 1E-6);
        assert_approx_eq!(lat_1deg_length(3.141592654, &wgs84_ellipsoid), 111314.502041079, 1E-6);
        assert_approx_eq!(lat_1deg_length(-3.141592654, &wgs84_ellipsoid), lat_1deg_length(3.141592654, &wgs84_ellipsoid), 1E-6);
        assert_approx_eq!(lat_1deg_length(0.000000000, &wgs84_ellipsoid), 111314.502041079, 1E-6);
        assert_approx_eq!(lat_1deg_length(0.000000000, &wgs84_ellipsoid), 111314.502041079, 1E-6);
        assert_approx_eq!(lat_1deg_length(0.174532925, &wgs84_ellipsoid), 111314.727675276, 1E-6);
        assert_approx_eq!(lat_1deg_length(-0.174532925, &wgs84_ellipsoid), 111314.727675276, 1E-6);
        assert_approx_eq!(lat_1deg_length(0.349065850, &wgs84_ellipsoid), 111315.377367309, 1E-6);
        assert_approx_eq!(lat_1deg_length(-0.349065850, &wgs84_ellipsoid), 111315.377367309, 1E-6);
        assert_approx_eq!(lat_1deg_length(0.523598776, &wgs84_ellipsoid), 111316.372765512, 1E-6);
        assert_approx_eq!(lat_1deg_length(-0.523598776, &wgs84_ellipsoid), 111316.372765512, 1E-6);
        assert_approx_eq!(lat_1deg_length(0.698131701, &wgs84_ellipsoid), 111317.593822430, 1E-6);
        assert_approx_eq!(lat_1deg_length(-0.698131701, &wgs84_ellipsoid), 111317.593822430, 1E-6);
        assert_approx_eq!(lat_1deg_length(0.872664626, &wgs84_ellipsoid), 111318.893268579, 1E-6);
        assert_approx_eq!(lat_1deg_length(-0.872664626, &wgs84_ellipsoid), 111318.893268579, 1E-6);
        assert_approx_eq!(lat_1deg_length(1.047197551, &wgs84_ellipsoid), 111320.114371577, 1E-6);
        assert_approx_eq!(lat_1deg_length(-1.047197551, &wgs84_ellipsoid), 111320.114371577, 1E-6);
        assert_approx_eq!(lat_1deg_length(1.221730476, &wgs84_ellipsoid), 111321.109840379, 1E-6);
        assert_approx_eq!(lat_1deg_length(-1.221730476, &wgs84_ellipsoid), 111321.109840379, 1E-6);
        assert_approx_eq!(lat_1deg_length(1.396263402, &wgs84_ellipsoid), 111321.759594497, 1E-6);
        assert_approx_eq!(lat_1deg_length(-1.396263402, &wgs84_ellipsoid), 111321.759594497, 1E-6);
        assert_approx_eq!(lat_1deg_length(1.570796327, &wgs84_ellipsoid), 111321.985253213, 1E-6);
        assert_approx_eq!(lat_1deg_length(-1.570796327, &wgs84_ellipsoid), 111321.985253213, 1E-6);

        for i in -90..90 {
            let i = (i as f64).to_radians();            
            assert_approx_eq!(lat_1deg_length(i, &wgs84_ellipsoid), lat_1deg_length(-i, &wgs84_ellipsoid), 1E-3);
        }
    }

    #[test]
    #[should_panic(expected = "Specified ellipsoid's major semiaxis (mjsa_m) should be greater than zero")]
    fn test_lat_1deg_length_non_positive_mjsa_panic() {
        let mut bad_ellipsoid: Ellipsoid = Ellipsoid::from_descriptor(&WGS84_ELLIPSOID_DESCRIPTOR);
        bad_ellipsoid.mjsa_m = -1.0;
        lat_1deg_length(-1.570796, &bad_ellipsoid);
    }

    #[test]
    #[should_panic(expected = "Specified ellipsoid's eccentrisity squared (ecnt_sq) should be in the range (0.0 .. 1.0) exclusively")]
    fn test_lat_1deg_length_ecnt_sq_out_of_range_panic() {
        let mut bad_ellipsoid: Ellipsoid = Ellipsoid::from_descriptor(&WGS84_ELLIPSOID_DESCRIPTOR);
        bad_ellipsoid.ecnt_sq = 1.0;
        lat_1deg_length(-1.570796, &bad_ellipsoid);
    }

    #[test]
    fn test_lon_1deg_length() {
        let wgs84_ellipsoid: Ellipsoid = Ellipsoid::from_descriptor(&WGS84_ELLIPSOID_DESCRIPTOR);

        assert_approx_eq!(lon_1deg_length(-3.141592654, &wgs84_ellipsoid), 111319.490793274, 1E-5);
        assert_approx_eq!(lon_1deg_length(3.141592654, &wgs84_ellipsoid), 111319.490793274, 1E-5);
        assert_approx_eq!(lon_1deg_length(-3.141592654, &wgs84_ellipsoid), lon_1deg_length(3.141592654, &wgs84_ellipsoid), 1E-5);
        
        assert_approx_eq!(lon_1deg_length(0.000000000, &wgs84_ellipsoid), 111319.490793274, 1E-4);
        assert_approx_eq!(lon_1deg_length(0.174532925, &wgs84_ellipsoid), 109628.371666625, 1E-4);
        assert_approx_eq!(lon_1deg_length(-0.174532925, &wgs84_ellipsoid), 109628.371666625, 1E-4);
        assert_approx_eq!(lon_1deg_length(0.349065850, &wgs84_ellipsoid), 104606.378238853, 1E-4);
        assert_approx_eq!(lon_1deg_length(-0.349065850, &wgs84_ellipsoid), 104606.378238853, 1E-4);
        assert_approx_eq!(lon_1deg_length(0.523598776, &wgs84_ellipsoid), 96406.047016128, 1E-4);
        assert_approx_eq!(lon_1deg_length(-0.523598776, &wgs84_ellipsoid), 96406.047016128, 1E-4);
        assert_approx_eq!(lon_1deg_length(0.698131701, &wgs84_ellipsoid), 85276.466841735, 1E-4);
        assert_approx_eq!(lon_1deg_length(-0.698131701, &wgs84_ellipsoid), 85276.466841735, 1E-4);
        assert_approx_eq!(lon_1deg_length(0.872664626, &wgs84_ellipsoid), 71555.730303868, 1E-4);
        assert_approx_eq!(lon_1deg_length(-0.872664626, &wgs84_ellipsoid), 71555.730303868, 1E-4);
        assert_approx_eq!(lon_1deg_length(1.047197551, &wgs84_ellipsoid), 55660.680811254, 1E-4);
        assert_approx_eq!(lon_1deg_length(-1.047197551, &wgs84_ellipsoid), 55660.680811254, 1E-4);
        assert_approx_eq!(lon_1deg_length(1.221730476, &wgs84_ellipsoid), 38074.261548399, 1E-4);
        assert_approx_eq!(lon_1deg_length(-1.221730476, &wgs84_ellipsoid), 38074.261548399, 1E-4);
        assert_approx_eq!(lon_1deg_length(1.396263402, &wgs84_ellipsoid), 19330.846811735, 1E-4);
        assert_approx_eq!(lon_1deg_length(-1.396263402, &wgs84_ellipsoid), 19330.846811735, 1E-4);
        assert_approx_eq!(lon_1deg_length(1.570796327, &wgs84_ellipsoid), 0.000000000, 1E-4);
        assert_approx_eq!(lon_1deg_length(-1.570796327, &wgs84_ellipsoid), 0.000000000, 1E-4);
        assert_approx_eq!(lon_1deg_length(1.745329252, &wgs84_ellipsoid), 19330.846811735, 1E-4);
        assert_approx_eq!(lon_1deg_length(-1.745329252, &wgs84_ellipsoid), 19330.846811735, 1E-4);
        assert_approx_eq!(lon_1deg_length(1.919862177, &wgs84_ellipsoid), 38074.261548399, 1E-4);
        assert_approx_eq!(lon_1deg_length(-1.919862177, &wgs84_ellipsoid), 38074.261548399, 1E-4);
        assert_approx_eq!(lon_1deg_length(2.094395102, &wgs84_ellipsoid), 55660.680811254, 1E-4);
        assert_approx_eq!(lon_1deg_length(-2.094395102, &wgs84_ellipsoid), 55660.680811254, 1E-4);
        assert_approx_eq!(lon_1deg_length(2.268928028, &wgs84_ellipsoid), 71555.730303868, 1E-4);
        assert_approx_eq!(lon_1deg_length(-2.268928028, &wgs84_ellipsoid), 71555.730303868, 1E-4);
        assert_approx_eq!(lon_1deg_length(2.443460953, &wgs84_ellipsoid), 85276.466841735, 1E-4);
        assert_approx_eq!(lon_1deg_length(-2.443460953, &wgs84_ellipsoid), 85276.466841735, 1E-4);
        assert_approx_eq!(lon_1deg_length(2.617993878, &wgs84_ellipsoid), 96406.047016128, 1E-4);
        assert_approx_eq!(lon_1deg_length(-2.617993878, &wgs84_ellipsoid), 96406.047016128, 1E-4);
        assert_approx_eq!(lon_1deg_length(2.792526803, &wgs84_ellipsoid), 104606.378238853, 1E-4);
        assert_approx_eq!(lon_1deg_length(-2.792526803, &wgs84_ellipsoid), 104606.378238853, 1E-4);
        assert_approx_eq!(lon_1deg_length(2.967059728, &wgs84_ellipsoid), 109628.371666625, 1E-4);
        assert_approx_eq!(lon_1deg_length(-2.967059728, &wgs84_ellipsoid), 109628.371666625, 1E-4);
        assert_approx_eq!(lon_1deg_length(3.141592654, &wgs84_ellipsoid), 111319.490793274, 1E-4);
        assert_approx_eq!(lon_1deg_length(-3.141592654, &wgs84_ellipsoid), 111319.490793274, 1E-4);

        for i in -180..180 {
            let i = (i as f64).to_radians();            
            assert_approx_eq!(lon_1deg_length(i, &wgs84_ellipsoid), lon_1deg_length(-i, &wgs84_ellipsoid), 1E-5);
        }
    }

    #[test]
    #[should_panic(expected = "Specified ellipsoid's major semiaxis (mjsa_m) should be greater than zero")]
    fn test_lon_1deg_length_non_positive_mjsa_panic() {
        let mut bad_ellipsoid: Ellipsoid = Ellipsoid::from_descriptor(&WGS84_ELLIPSOID_DESCRIPTOR);
        bad_ellipsoid.mjsa_m = -1.0;
        lon_1deg_length(-1.570796, &bad_ellipsoid);
    }

    #[test]
    #[should_panic(expected = "Specified ellipsoid's eccentrisity squared (ecnt_sq) should be in the range (0.0 .. 1.0) exclusively")]
    fn test_lon_1deg_length_ecnt_sq_out_of_range_panic() {
        let mut bad_ellipsoid: Ellipsoid = Ellipsoid::from_descriptor(&WGS84_ELLIPSOID_DESCRIPTOR);
        bad_ellipsoid.ecnt_sq = 1.0;
        lon_1deg_length(-1.570796, &bad_ellipsoid);
    }

    #[test]
    fn test_geopoint_offset_by_deltas_and_get_deltas_by_geopoints() {
        
        let el: Ellipsoid = Ellipsoid::from_descriptor(&WGS84_ELLIPSOID_DESCRIPTOR);
        let mut lat_rad: f64;
        let mut lon_rad: f64;

        let mut gcs2lcs: (f64, f64);
        let mut lcs2gcs: (f64, f64);

        for lat_deg in -8..8 {
            lat_rad = (lat_deg as f64 * 10.0).to_radians();
            for lon_deg in -17..17 {               
               lon_rad = (lon_deg as f64 * 10.0).to_radians();

               for delta_m in -1..3 {
                   let delta_m = (10.0 as f64).powi(delta_m);
               
                   gcs2lcs = geopoint_offset_by_deltas(lat_rad, lon_rad, delta_m, delta_m, &el);
                   lcs2gcs = get_deltas_by_geopoints(lat_rad, lon_rad, gcs2lcs.0, gcs2lcs.1, &el);

                   assert_approx_eq!(lcs2gcs.0, delta_m, 1E-4 * delta_m);
                   assert_approx_eq!(lcs2gcs.1, delta_m, 1E-4 * delta_m);
               }
            }           
        }
    }

    #[test]
    fn test_geopoint_offset_by_deltas_and_get_deltas_by_geopoints_wgs84() {
        let mut lat_rad: f64;
        let mut lon_rad: f64;
        let mut gcs2lcs: (f64, f64);
        let mut lcs2gcs: (f64, f64);

        for lat_deg in -8..8 {
            lat_rad = (lat_deg as f64 * 10.0).to_radians();
            for lon_deg in -17..17 {               
               lon_rad = (lon_deg as f64 * 10.0).to_radians();
               for delta_m in -1..3 {
                   let delta_m = (10.0 as f64).powi(delta_m);               
                   gcs2lcs = geopoint_offset_by_deltas_wgs84(lat_rad, lon_rad, delta_m, delta_m);
                   lcs2gcs = get_deltas_by_geopoints_wgs84(lat_rad, lon_rad, gcs2lcs.0, gcs2lcs.1);
                   assert_approx_eq!(lcs2gcs.0, delta_m, 1E-4 * delta_m);
                   assert_approx_eq!(lcs2gcs.1, delta_m, 1E-4 * delta_m);
               }
            }           
        }
    }

    #[test]
    fn test_haversine_inverse_and_direct_and_initial_bearing() {

        let mut sp_lat_rad: f64;
        let mut sp_lon_rad: f64;
        let mut ep_lat_rad: f64;
        let mut ep_lon_rad: f64;
        let mut fwd_az_rad: f64;
        let mut a_az_rad: f64;
        let mut d_point: (f64, f64);
        let mut adist_m: f64;
    
        for lat_deg in -8..8 {
            sp_lat_rad = (lat_deg as f64 * 10.0).to_radians();
            for lon_deg in -17..17 {
                sp_lon_rad = (lon_deg as f64 * 10.0).to_radians();
                for az_deg in 0..35 {
                    fwd_az_rad = wrap_2pi((az_deg as f64 * 10.0).to_radians());
                    for dist_m in -1..3 {
                        let dst_m = (10.0 as f64).powi(dist_m);
                        d_point = haversine_direct(sp_lat_rad, sp_lon_rad, dst_m, fwd_az_rad, WGS84_ELLIPSOID_DESCRIPTOR.mjsa_m);
                        ep_lat_rad = d_point.0;
                        ep_lon_rad = d_point.1;
                        adist_m = haversine_inverse(sp_lat_rad, sp_lon_rad, ep_lat_rad, ep_lon_rad, WGS84_ELLIPSOID_DESCRIPTOR.mjsa_m);
                        a_az_rad = haversine_initial_bearing(sp_lat_rad, sp_lon_rad, ep_lat_rad, ep_lon_rad);

                        if ((a_az_rad - fwd_az_rad).abs() - PI2).abs() < 10E-3 {
                            assert_approx_eq!((a_az_rad - fwd_az_rad).abs(), PI2, 10E-3);                            
                        } else {
                            assert_approx_eq!(fwd_az_rad, a_az_rad, 10E-3);
                        }

                        assert_approx_eq!(dst_m, adist_m, 10E-2);                        
                    }
                }
            }
        }        
    }

    #[test]
    fn test_vincenty_equations() {

        let el: Ellipsoid = Ellipsoid::from_descriptor(&WGS84_ELLIPSOID_DESCRIPTOR);
        let mut sp_lat_rad: f64;
        let mut sp_lon_rad: f64;
        let mut ep_lat_rad: f64;
        let mut ep_lon_rad: f64;
        let mut fwd_az_rad: f64;
        let mut rev_az_rad: f64;
        let mut a_az_rad: f64;
        let mut a_raz_rad: f64;
        let mut d_point: (f64, f64, f64, i32);
        let mut dd_point: (f64, f64, f64, i32, bool);
        let mut adist_m: f64;        
    
        for lat_deg in -9..9 {
            sp_lat_rad = (lat_deg as f64 * 10.0).to_radians();
            for lon_deg in -18..18 {
                sp_lon_rad = (lon_deg as f64 * 10.0).to_radians();
                for az_deg in 0..35 {
                    fwd_az_rad = wrap_2pi((az_deg as f64 * 10.0).to_radians());
                    for dist_m in -1..3 {
                        let dst_m = (10.0 as f64).powi(dist_m);
                                                                        
                        d_point = vincenty_direct(sp_lat_rad, sp_lon_rad, fwd_az_rad, dst_m, &el, VNC_DEF_EPSILON, VNC_DEF_IT_LIMIT);
                                                
                        ep_lat_rad = d_point.0;
                        ep_lon_rad = d_point.1;
                        rev_az_rad = d_point.2;

                        dd_point = vincenty_inverse(sp_lat_rad, sp_lon_rad, ep_lat_rad, ep_lon_rad, &el, VNC_DEF_EPSILON, VNC_DEF_IT_LIMIT);
                        adist_m = dd_point.0;
                        a_az_rad = dd_point.1;
                        a_raz_rad = dd_point.2;

                        if ((a_az_rad - fwd_az_rad).abs() - PI2).abs() < 10E-3 {
                            assert_approx_eq!((a_az_rad - fwd_az_rad).abs(), PI2, 10E-6);                            
                        } else {
                            assert_approx_eq!(fwd_az_rad, a_az_rad, 10E-6);
                        }

                        assert_approx_eq!(rev_az_rad, a_raz_rad, 10E-6);                      
                        assert_approx_eq!(dst_m, adist_m, 10E-6);                        
                    }
                }
            }
        }        




    }

    #[test]
    fn test_toa_nlm_2d_solve_no_noise() {

        let mut base_points = Vec::new();
    
        // random relevant values for actual point position
        let actual_x = 10.5;
        let actual_y = -11.0;
        let actual_z = 18.3;

        let mut x;
        let mut y;
        let mut z;
        let mut d;
            
        x = -100.0;
        y = 200.0;
        z = 11.0;
        d = dist_3d(x, y, z, actual_x, actual_y, actual_z);
        base_points.push(TOABasePoint {x: x, y: y, z: z, d: d});

        x = 115.0;
        y = 267.0;
        z = 17.3;
        d = dist_3d(x, y, z, actual_x, actual_y, actual_z);
        base_points.push(TOABasePoint {x: x, y: y, z: z, d: d});

        x = 315.4;
        y = -118.4;
        z = 31.2;
        d = dist_3d(x, y, z, actual_x, actual_y, actual_z);
        base_points.push(TOABasePoint {x: x, y: y, z: z, d: d});

        x = -470.1;
        y = -216.7;
        z = 12.5;
        d = dist_3d(x, y, z, actual_x, actual_y, actual_z);
        base_points.push(TOABasePoint {x: x, y: y, z: z, d: d});

        
        let result = toa_nlm_2d_solve(&base_points, 0.0, 0.0, actual_z, NLM_DEF_IT_LIMIT, NLM_DEF_PREC_THRLD, 1.0);

        assert_approx_eq!(result.0, actual_x, 10E-4);
        assert_approx_eq!(result.1, actual_y, 10E-4);
        
        assert!(result.2 < 1.0, "Residual function greater than limit: {}", result.2);
        assert!(result.3 < NLM_DEF_IT_LIMIT, "Method did not converge: iterations limit exeeded {}", result.3);
    }

    #[test]
    fn test_toa_nlm_2d_solve_noise() {

        let mut base_points = Vec::new();
        let err_h_amp = 1.5; // distance measurement error in meters
            
        // random relevant values for actual point position
        let actual_x = 10.5;
        let actual_y = -11.0;
        let actual_z = 18.3;

        let mut x;
        let mut y;
        let mut z;
        let mut d;

        let mut d_err: f64;
            
        x = -100.0;
        y = 200.0;
        z = 11.0;
        d_err = thread_rng().gen();
        d = dist_3d(x, y, z, actual_x, actual_y, actual_z) + d_err * err_h_amp * 2.0 - err_h_amp;
        base_points.push(TOABasePoint {x: x, y: y, z: z, d: d});

        x = 115.0;
        y = 267.0;
        z = 17.3;
        d_err = thread_rng().gen();
        d = dist_3d(x, y, z, actual_x, actual_y, actual_z) + d_err * err_h_amp * 2.0 - err_h_amp;
        base_points.push(TOABasePoint {x: x, y: y, z: z, d: d});

        x = 315.4;
        y = -118.4;
        z = 31.2;
        d_err = thread_rng().gen();
        d = dist_3d(x, y, z, actual_x, actual_y, actual_z) + d_err * err_h_amp * 2.0 - err_h_amp;
        base_points.push(TOABasePoint {x: x, y: y, z: z, d: d});

        x = -470.1;
        y = -216.7;
        z = 12.5;
        d_err = thread_rng().gen();
        d = dist_3d(x, y, z, actual_x, actual_y, actual_z) + d_err * err_h_amp * 2.0 - err_h_amp;
        base_points.push(TOABasePoint {x: x, y: y, z: z, d: d});
        
        let result = toa_nlm_2d_solve(&base_points, 0.0, 0.0, actual_z, NLM_DEF_IT_LIMIT, NLM_DEF_PREC_THRLD, 1.0);

        assert_approx_eq!(result.0, actual_x, err_h_amp * 2.0);
        assert_approx_eq!(result.1, actual_y, err_h_amp * 2.0);
        
        assert!(result.2 < err_h_amp * 2.0, "Residual function greater than limit: {}", result.2);
        assert!(result.3 < NLM_DEF_IT_LIMIT, "Method did not converge: iterations limit exeeded {}", result.3);
    }

    #[test]
    fn test_tdoa_nlm_2d_solve_no_noise() {

        let mut base_lines = Vec::new();
                    
        // random relevant values for actual point position
        let actual_x = 10.5;
        let actual_y = -11.0;
        let actual_z = 18.3;               
            
        let x1 = -100.0;
        let y1 = 200.0;
        let z1 = 11.0;        
        let d1 = dist_3d(x1, y1, z1, actual_x, actual_y, actual_z);
        let x2 = 115.0;
        let y2 = 267.0;
        let z2 = 17.3;        
        let d2 = dist_3d(x2, y2, z2, actual_x, actual_y, actual_z);
        let x3 = 315.4;
        let y3 = -118.4;
        let z3 = 31.2;        
        let d3 = dist_3d(x3, y3, z3, actual_x, actual_y, actual_z);                
        let x4 = -470.1;
        let y4 = -216.7;
        let z4 = 12.5;        
        let d4 = dist_3d(x4, y4, z4, actual_x, actual_y, actual_z);
        

        base_lines.push(TDOABaseline { x1: x1, y1: y1, z1: z1, x2: x2, y2: y2, z2: z2, prd: d1 - d2 });
        base_lines.push(TDOABaseline { x1: x1, y1: y1, z1: z1, x2: x3, y2: y3, z2: z3, prd: d1 - d3 });
        base_lines.push(TDOABaseline { x1: x1, y1: y1, z1: z1, x2: x4, y2: y4, z2: z4, prd: d1 - d4 });
        base_lines.push(TDOABaseline { x1: x2, y1: y2, z1: z2, x2: x3, y2: y3, z2: z3, prd: d2 - d3 });
        base_lines.push(TDOABaseline { x1: x2, y1: y2, z1: z2, x2: x4, y2: y4, z2: z4, prd: d2 - d4 });        
        base_lines.push(TDOABaseline { x1: x3, y1: y3, z1: z3, x2: x4, y2: y4, z2: z4, prd: d3 - d4 });
        
        let result = tdoa_nlm_2d_solve(&base_lines, 0.0, 0.0, actual_z, NLM_DEF_IT_LIMIT, NLM_DEF_PREC_THRLD, 1.0);

        assert_approx_eq!(result.0, actual_x, 10E-4);
        assert_approx_eq!(result.1, actual_y, 10E-4);
        
        assert!(result.2 < 1.0, "Residual function greater than limit: {}", result.2);
        assert!(result.3 < NLM_DEF_IT_LIMIT, "Method did not converge: iterations limit exeeded {}", result.3);
    }   

    // 3D solutions - not completed!
    #[test]
    fn test_toa_nlm_3d_solve_no_noise() {

        let mut base_points = Vec::new();
    
        // random relevant values for actual point position
        let actual_x = 10.5;
        let actual_y = -11.0;
        let actual_z = 18.3;

        let mut x;
        let mut y;
        let mut z;
        let mut d;        

        let x_prev = 0.0;
        let y_prev = 0.0;
        let z_prev = 0.0;
            
        x = -100.0;
        y = 200.0;
        z = 10.0;
        d = dist_3d(x, y, z, actual_x, actual_y, actual_z);
        base_points.push(TOABasePoint {x: x, y: y, z: z, d: d});            

        x = 115.0;
        y = 267.0;
        z = 17.3;
        d = dist_3d(x, y, z, actual_x, actual_y, actual_z);
        base_points.push(TOABasePoint {x: x, y: y, z: z, d: d});
        
        x = 315.4;
        y = -118.4;
        z = 31.2;
        d = dist_3d(x, y, z, actual_x, actual_y, actual_z);
        base_points.push(TOABasePoint {x: x, y: y, z: z, d: d});
        
        x = -470.1;
        y = -216.7;
        z = 12.5;
        d = dist_3d(x, y, z, actual_x, actual_y, actual_z);
        base_points.push(TOABasePoint {x: x, y: y, z: z, d: d});
                
        let result = toa_nlm_3d_solve(&base_points, x_prev, y_prev, z_prev, NLM_DEF_IT_LIMIT, NLM_DEF_PREC_THRLD, 10.0);

        assert_approx_eq!(result.0, actual_x, 10.0);
        assert_approx_eq!(result.1, actual_y, 10.0);
        assert_approx_eq!(result.2, actual_z, 10.0);
        
        assert!(result.3 < 10.0, "Residual function greater than limit: {}", result.3);
        assert!(result.4 < NLM_DEF_IT_LIMIT, "Method did not converge: iterations limit exeeded {}", result.4);
    }

    #[test]
    fn test_tdoa_nlm_3d_solve_no_noise() {

        let mut base_lines = Vec::new();
                    
        // random relevant values for actual point position
        let actual_x = 10.5;
        let actual_y = -11.0;
        let actual_z = 18.3;           
        
        let x_prev = 0.0;
        let y_prev = 0.0;
        let z_prev = 0.0;
            
        let x1 = -100.0;
        let y1 = 200.0;
        let z1 = 11.0;        
        let d1 = dist_3d(x1, y1, z1, actual_x, actual_y, actual_z);
        let x2 = 115.0;
        let y2 = 267.0;
        let z2 = 17.3;        
        let d2 = dist_3d(x2, y2, z2, actual_x, actual_y, actual_z);
        let x3 = 315.4;
        let y3 = -118.4;
        let z3 = 31.2;        
        let d3 = dist_3d(x3, y3, z3, actual_x, actual_y, actual_z);                
        let x4 = -470.1;
        let y4 = -216.7;
        let z4 = 12.5;        
        let d4 = dist_3d(x4, y4, z4, actual_x, actual_y, actual_z);        

        base_lines.push(TDOABaseline { x1: x1, y1: y1, z1: z1, x2: x2, y2: y2, z2: z2, prd: d1 - d2 });
        base_lines.push(TDOABaseline { x1: x1, y1: y1, z1: z1, x2: x3, y2: y3, z2: z3, prd: d1 - d3 });
        base_lines.push(TDOABaseline { x1: x1, y1: y1, z1: z1, x2: x4, y2: y4, z2: z4, prd: d1 - d4 });
        base_lines.push(TDOABaseline { x1: x2, y1: y2, z1: z2, x2: x3, y2: y3, z2: z3, prd: d2 - d3 });
        base_lines.push(TDOABaseline { x1: x2, y1: y2, z1: z2, x2: x4, y2: y4, z2: z4, prd: d2 - d4 });        
        base_lines.push(TDOABaseline { x1: x3, y1: y3, z1: z3, x2: x4, y2: y4, z2: z4, prd: d3 - d4 });
        
        let result = tdoa_nlm_3d_solve(&base_lines, x_prev, y_prev, z_prev, NLM_DEF_IT_LIMIT, NLM_DEF_PREC_THRLD, 10.0);

        assert_approx_eq!(result.0, actual_x, 20.0);
        assert_approx_eq!(result.1, actual_y, 20.0);
        assert_approx_eq!(result.2, actual_z, 20.0);
        
        assert!(result.3 < 10.0, "Residual function greater than limit: {}", result.3);
        assert!(result.4 < NLM_DEF_IT_LIMIT, "Method did not converge: iterations limit exeeded {}", result.4);
    }   
}
