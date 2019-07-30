using System;

namespace UCNLNav
{
    public class Ellipsoid
    {
        public double MajorSemiAxis_m { get; private set; }
        public double MinorSemiAxis_m { get; private set; }
        public double Flattening { get; private set; }
        public double InverseFlattening { get; private set; }
        public double Eccentricity { get; private set; }
        public double EccentricitySq { get; private set; }

        public Ellipsoid(double mjsa, double inverseFlattening)
        {
            MajorSemiAxis_m = mjsa;
            InverseFlattening = inverseFlattening;

            Flattening = 1 / inverseFlattening;
            MinorSemiAxis_m = MajorSemiAxis_m * (1 - Flattening);
            Eccentricity = ((MajorSemiAxis_m * MajorSemiAxis_m) - (MinorSemiAxis_m * MinorSemiAxis_m)) / (MajorSemiAxis_m * MajorSemiAxis_m);
            EccentricitySq = Eccentricity * Eccentricity;
        }
    }

    public class Algorithms
    {
        public static readonly double PI2 = Math.PI * 2;
        public static readonly double PI_DBY_180 = Math.PI / 180.0;
        public static readonly double D180_DBY_PI = 180.0 / Math.PI;

        public static readonly Ellipsoid WGS84Ellipsoid = new Ellipsoid(6378137, 298.257223563);

        public static double Deg2Rad(double angle_rad)
        {
            return angle_rad * PI_DBY_180;
        }

        public static double Rad2Deg(double angle_deg)
        {
            return angle_deg * D180_DBY_PI;
        }
        
        public static double Wrap(double val, double lim)
        {
            double vl = Math.Abs(val);
            double sign = Math.Sign(val);

            while (vl > lim)
                vl -= lim;

            return vl * sign;
        }

        public static double Wrap2PI(double angle_rad)
        {
            return Wrap(angle_rad, PI2);
        }

        public static double Wrap360(double angle_deg)
        {
            return Wrap(angle_deg, 360);
        }



        public static double HaversineInverse(double sp_lat_rad, double sp_lon_rad, double ep_lat_rad, double ep_lon_rad, double e_radius_m)
        {
            double dLat = ep_lat_rad - sp_lat_rad;
            double dLon = ep_lon_rad - sp_lon_rad;
            double a = Math.Pow(Math.Sin(dLat / 2), 2) +
                      Math.Cos(sp_lat_rad) * Math.Cos(ep_lat_rad) *
                      Math.Sin(dLon / 2) * Math.Sin(dLon / 2);
            return e_radius_m * 2 * Math.Atan2(Math.Sqrt(a), Math.Sqrt(1 - a));
        }

        public static void HaversineDirect(double sp_lat_rad, double sp_lon_rad, double dst_m, double fwd_az_rad, double e_radius_m,
            out double ep_lat_rad, out double ep_lon_rad)
        {
            double delta = dst_m / e_radius_m;
            ep_lat_rad = Math.Asin(Math.Sin(sp_lat_rad) * Math.Cos(delta) + Math.Cos(sp_lat_rad) * Math.Sin(delta) * Math.Cos(fwd_az_rad));
            ep_lon_rad = Wrap2PI(3 * Math.PI + (sp_lon_rad + Math.Atan2(Math.Sin(fwd_az_rad) * Math.Sin(delta) * Math.Cos(sp_lat_rad),
                Math.Cos(delta) - Math.Sin(sp_lat_rad) * Math.Sin(ep_lat_rad)))) - Math.PI;
        }

        public static double HaversineInitialBearing(double sp_lat_rad, double sp_lon_rad, double ep_lat_rad, double ep_lon_rad)
        {
            double y = Math.Sin(ep_lon_rad - sp_lon_rad) * Math.Cos(ep_lat_rad);
            double x = Math.Cos(sp_lat_rad) * Math.Sin(ep_lat_rad) - Math.Sin(sp_lat_rad) * Math.Cos(ep_lat_rad) * Math.Cos(ep_lon_rad - sp_lon_rad);
            return Wrap2PI(Math.PI + Math.Atan2(y, x));
        }

        public static double HaversineFinalBearing(double sp_lat_rad, double sp_lon_rad, double ep_lat_rad, double ep_lon_rad)
        {
            return Wrap2PI(HaversineInitialBearing(ep_lat_rad, ep_lon_rad, sp_lat_rad, sp_lon_rad) + Math.PI);
        }



        public static double Lat1DegLength(double lat_rad, Ellipsoid el)
        {
            return PI_DBY_180 * (1 - el.EccentricitySq) / Math.Pow(1 - el.EccentricitySq * Math.Sin(lat_rad) * Math.Sin(lat_rad), 1.5);
        }

        public static double Lon1DegLength(double lat_rad, Ellipsoid el)
        {
            return PI_DBY_180 * el.MajorSemiAxis_m * Math.Cos(lat_rad) / Math.Sqrt(1 - el.EccentricitySq * Math.Sin(lat_rad) * Math.Sin(lat_rad));
        }



        public static void GeopointOffsetByDeltas(double lat_rad, double lon_rad, double lat_offset_m, double lon_offset_m,
            Ellipsoid el,
            out double e_lat_rad, out double e_lon_rad)
        {
            double m_per_deg_lat = Lat1DegLength(lat_rad, el);
            double m_per_deg_lon = Lon1DegLength(lat_rad, el);
            e_lat_rad = lat_rad - PI_DBY_180 * lat_offset_m / m_per_deg_lat;
            e_lon_rad = lon_rad - PI_DBY_180 * lon_offset_m / m_per_deg_lon;
        }

        public static void GeopointOffsetByDeltas_WGS84(double lat_rad, double lon_rad, double lat_offset_m, double lon_offset_m,
            out double e_lat_rad, out double e_lon_rad)
        {
            double m_per_deg_lat = 111132.92 - 559.82 * Math.Cos(2.0 * lat_rad) + 1.175 * Math.Cos(4.0 * lat_rad);
            double m_per_deg_lon = 111412.84 * Math.Cos(lat_rad) - 93.5 * Math.Cos(3.0 * lat_rad);
            e_lat_rad = lat_rad - PI_DBY_180 * lat_offset_m / m_per_deg_lat;
            e_lon_rad = lon_rad - PI_DBY_180 * lon_offset_m / m_per_deg_lon;
        }

        public static void GetDeltasByGeopoints(double sp_lat_rad, double sp_lon_rad, double ep_lat_rad, double ep_lon_rad,
            Ellipsoid el,
            out double delta_lat_m, out double delta_lon_m)
        {
            double m_lat_rad = (sp_lat_rad + ep_lat_rad) / 2.0;
            double m_per_deg_lat = Lat1DegLength(m_lat_rad, el);
            double m_per_deg_lon = Lon1DegLength(m_lat_rad, el);
            delta_lat_m = (sp_lat_rad - ep_lat_rad) * m_per_deg_lat * D180_DBY_PI;
            delta_lon_m = (sp_lon_rad - ep_lon_rad) * m_per_deg_lon * D180_DBY_PI;
        }

        public static void GetDeltasByGeopoints_WGS84(double sp_lat_rad, double sp_lon_rad, double ep_lat_rad, double ep_lon_rad,
            out double delta_lat_m, out double delta_lon_m)
        {
            double m_lat_rad = (sp_lat_rad + ep_lat_rad) / 2.0;
            double m_per_deg_lat = 111132.92 - 559.82 * Math.Cos(2.0 * m_lat_rad) + 1.175 * Math.Cos(4.0 * m_lat_rad);
            double m_per_deg_lon = 111412.84 * Math.Cos(m_lat_rad) - 93.5 * Math.Cos(3.0 * m_lat_rad);
            delta_lat_m = (sp_lat_rad - ep_lat_rad) * m_per_deg_lat * D180_DBY_PI;
            delta_lon_m = (sp_lon_rad - ep_lon_rad) * m_per_deg_lon * D180_DBY_PI;
        }



        public static bool VincentyInverse(double sp_lat_rad, double sp_lon_rad, double ep_lat_rad, double ep_lon_rad,
            Ellipsoid el, double eps, int it_limit,
            out double dst_m, out double fwd_az_rad, out double rev_az_rad,
            out int its)
        {
            double l_ = ep_lon_rad - sp_lon_rad;
            double tan_u_1 = (1.0 - el.Flattening) * Math.Tan(sp_lat_rad);
            double cos_u_1 = 1.0 / Math.Sqrt((1.0 + tan_u_1 * tan_u_1));
            double sin_u_1 = tan_u_1 * cos_u_1;

            double tan_u_2 = (1.0 - el.Flattening) * Math.Tan(ep_lat_rad);
            double cos_u_2 = 1.0 / Math.Sqrt((1.0 + tan_u_2 * tan_u_2));
            double sin_u_2 = tan_u_2 * cos_u_2;

            double sin_lambda = 0;
            double cos_lambda = 0;
            double sin_sigma = 0;
            double cos_sigma = 0;
            double sin_alpha = 0;
            double sin_sq_sigma = 0;
            double cos_sq_alpha = 0;
            double cos_2_sigma_m = 0;
            double sigma = 0;
            double c_ = 0;

            double lambda = l_;
            double lambda_ = 0;
            its = 0;
            double it_check = 0;

            bool antimeridian = Math.Abs(l_) > Math.PI;

            do
            {
                sin_lambda = Math.Sin(lambda);
                cos_lambda = Math.Cos(lambda);

                sin_sq_sigma = (cos_u_2 * sin_lambda) * (cos_u_2 * sin_lambda) + (cos_u_1 * sin_u_2 - sin_u_1 * cos_u_2 * cos_lambda) * (cos_u_1 * sin_u_2 - sin_u_1 * cos_u_2 * cos_lambda);

                if (Math.Abs(sin_sq_sigma) > double.Epsilon)
                {
                    sin_sigma = Math.Sqrt(sin_sq_sigma);
                    cos_sigma = sin_u_1 * sin_u_2 + cos_u_1 * cos_u_2 * cos_lambda;
                    sigma = Math.Atan2(sin_sigma, cos_sigma);
                    sin_alpha = cos_u_1 * cos_u_2 * sin_lambda / sin_sigma;

                    cos_sq_alpha = 1 - sin_alpha * sin_alpha;

                    if (cos_sq_alpha != 0)
                        cos_2_sigma_m = cos_sigma - 2 * sin_u_1 * sin_u_2 / cos_sq_alpha;
                    else
                        cos_sq_alpha = 0;

                    c_ = el.Flattening / 16.0 * cos_sq_alpha * (4 + el.Flattening * (4 - 3 * cos_sq_alpha));
                    lambda_ = lambda;
                    lambda = l_ + (1 - c_) * el.Flattening * sin_alpha * (sigma + c_ * sin_sigma * (cos_2_sigma_m + c_ * cos_sigma * (-1 + 2 * cos_2_sigma_m * cos_2_sigma_m)));

                    if (antimeridian)
                        it_check = Math.Abs(lambda) - Math.PI;
                    else
                        it_check = Math.Abs(lambda);
                }
                
            } while ((Math.Abs(lambda - lambda_) > eps) && (++its < it_limit) && (it_check < Math.PI));


            double u_sq = cos_sq_alpha * (el.MajorSemiAxis_m * el.MajorSemiAxis_m - el.MinorSemiAxis_m * el.MajorSemiAxis_m) / (el.MinorSemiAxis_m * el.MinorSemiAxis_m);
            double a_ = 1 + u_sq / 16384 * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)));
            double b_ = u_sq / 1024 * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)));
            double delta_sigma = b_ * sin_sigma * (cos_2_sigma_m + b_ / 4 * (cos_sigma * (-1 + 2 * cos_2_sigma_m * cos_2_sigma_m) -
                b_ / 6 * cos_2_sigma_m * (-3 + 4 * sin_sigma * sin_sigma) * (-3 + 4 * cos_2_sigma_m * cos_2_sigma_m)));


            dst_m = el.MinorSemiAxis_m * a_ * (sigma - delta_sigma);

            fwd_az_rad = Math.Atan2(cos_u_2 * sin_lambda, cos_u_1 * sin_u_2 - sin_u_1 * cos_u_2 * cos_lambda);
            rev_az_rad = Math.Atan2(cos_u_1 * sin_lambda, -sin_u_1 * cos_u_2 + cos_u_1 * sin_u_2 * cos_lambda);

            return (its < it_limit) && (it_check < Math.PI);
        }

        public static bool VincentyDirect(double sp_lat_rad, double sp_lon_rad, double fwd_az_rad, double dst_m,
            Ellipsoid el, double eps, int it_limit,
            out double ep_lat_rad, out double ep_lon_rad, out double rev_az_rad,
            out int its)
        {
            double sin_alpha_1 = Math.Sin(fwd_az_rad);
            double cos_alpha_1 = Math.Cos(fwd_az_rad);
            double tan_u_1 = (1.0 - el.Flattening) * Math.Tan(sp_lat_rad);
            double cos_u_1 = 1.0 / Math.Sqrt(1.0 + tan_u_1 * tan_u_1);
            double sin_u_1 = tan_u_1 * cos_u_1;

            double sigma_1 = Math.Atan2(tan_u_1, cos_alpha_1);
            double sin_alpha = cos_u_1 * sin_alpha_1;
            double cos_sq_alpha = 1.0 - sin_alpha * sin_alpha;
            double u_sq = cos_sq_alpha * (el.MajorSemiAxis_m * el.MajorSemiAxis_m - el.MinorSemiAxis_m * el.MinorSemiAxis_m) / (el.MinorSemiAxis_m * el.MinorSemiAxis_m);
            double a_ = 1.0 + u_sq / 16384.0 * (4096.0 + u_sq * (-768.0 + u_sq * (320.0 - 175.0 * u_sq)));
            double b_ = u_sq / 1024.0 * (256.0 + u_sq * (-128.0 + u_sq * (74.0 - 47.0 * u_sq)));

            double cos_2_sigma_m;
            double sin_sigma;
            double cos_sigma;
            double delta_sigma;

            double sigma = dst_m / (el.MinorSemiAxis_m * a_);
            double sigma_;
            its = 0;

            do
            {
                cos_2_sigma_m = Math.Cos(2.0 * sigma_1 + sigma);
                sin_sigma = Math.Sin(sigma);
                cos_sigma = Math.Cos(sigma);

                delta_sigma = b_ * sin_sigma * (cos_2_sigma_m + b_ / 4.0 * (cos_sigma * (-1.0 + 2.0 * cos_2_sigma_m * cos_2_sigma_m) -
                      b_ / 6.0 * cos_2_sigma_m * (-3.0 + 4.0 * sin_sigma * sin_sigma) * (-3.0 + 4.0 * cos_2_sigma_m * cos_2_sigma_m)));

                sigma_ = sigma;
                sigma = dst_m / (el.MinorSemiAxis_m * a_) + delta_sigma;

                its++;
            } while ((Math.Abs(sigma - sigma_) > eps) && (its < it_limit));

            double x = sin_u_1 * sin_sigma - cos_u_1 * cos_sigma * cos_alpha_1;
            ep_lat_rad = Math.Atan2(sin_u_1 * cos_sigma + cos_u_1 * sin_sigma * cos_alpha_1, (1.0 - el.Flattening) * Math.Sqrt(sin_alpha * sin_alpha + x * x));

            double lambda = Math.Atan2(sin_sigma * sin_alpha_1, cos_u_1 * cos_sigma - sin_u_1 * sin_sigma * cos_alpha_1);
            double c_ = el.Flattening / 16.0 * cos_sq_alpha * (4.0 + el.Flattening * (4.0 - 3.0 * cos_sq_alpha));

            double l_ = lambda - (1.0 - c_) * el.Flattening * sin_alpha * (sigma + c_ * sin_sigma * (cos_2_sigma_m + c_ * cos_sigma * (-1.0 + 2.0 * cos_2_sigma_m * cos_2_sigma_m)));

            ep_lon_rad = sp_lon_rad + l_;
            rev_az_rad = Math.Atan2(sin_alpha, -x);

            return (its < it_limit);
        }
    }
}