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
        #region Nelder-Mead constants

        static readonly double NLM_A = 1.0;
        static readonly double NLM_B = 0.5;
        static readonly double NLM_R = 0.5;
        static readonly double NLM_Q = 0.5;
        static readonly double NLM_G = 2.0;

        public static readonly int    NLM_DEF_IT_LIMIT = 600;
        public static readonly double NLM_DEF_PREC_THRLD = 1E-8;

        #endregion

        #region Commonly used reference ellipsoids

        public static readonly Ellipsoid WGS72Ellipsoid = new Ellipsoid(6378135, 298.26);
        public static readonly Ellipsoid WGS84Ellipsoid = new Ellipsoid(6378137, 298.257223563);
        public static readonly Ellipsoid GRS80Ellipsoid = new Ellipsoid(6378137, 298.257222100882711);
        public static readonly Ellipsoid PZ90Ellipsoid = new Ellipsoid(6378136, 298.257839303);
        public static readonly Ellipsoid IERSEllipsoid = new Ellipsoid(6378136.49, 298.25645);
        public static readonly Ellipsoid KRSYEllipsoid = new Ellipsoid(6378245, 298.3);

        #endregion

        public static readonly double PI2 = Math.PI * 2.0;
        public static readonly double PI_DBY_180 = Math.PI / 180.0;
        public static readonly double D180_DBY_PI = 180.0 / Math.PI;

        public static readonly int    VNC_DEF_IT_LIMIT = 200;
        public static readonly double VNC_DEF_EPSILON = 1E-12;        

        #region Methods

        #region Utilities

        /// <summary>
        /// Converts degrees to radians
        /// </summary>
        /// <param name="angle_deg">Angle in degrees</param>
        /// <returns>Angle in radians</returns>
        public static double Deg2Rad(double angle_deg)
        {
            return angle_deg * PI_DBY_180;
        }

        /// <summary>
        /// Converts radians to degrees
        /// </summary>
        /// <param name="angle_rad">Angle in radians</param>
        /// <returns>Angle in degrees</returns>
        public static double Rad2Deg(double angle_rad)
        {
            return angle_rad * D180_DBY_PI;
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

        #endregion

        #region Degrees and meters conversion routines

        /// <summary>
        /// Calculates length in meter of 1 degree latitude according to specified ellipsoid
        /// </summary>
        /// <param name="lat_rad">Latitude, signed degrees</param>
        /// <param name="el">Ellipsoid</param>
        /// <returns>Length in meters</returns>
        public static double Lat1DegLength(double lat_rad, Ellipsoid el)
        {
            return PI_DBY_180 * el.MajorSemiAxis_m * (1 - el.EccentricitySq) / Math.Pow(1 - el.EccentricitySq * Math.Sin(lat_rad) * Math.Sin(lat_rad), 1.5);
        }

        /// <summary>
        /// Calculates length in meter of 1 degree longitude according to specified ellipsoid and latitude
        /// </summary>
        /// <param name="lat_rad">Latitude, signed degrees</param>
        /// <param name="el">Ellipsoid</param>
        /// <returns>Length in meters</returns>
        public static double Lon1DegLength(double lat_rad, Ellipsoid el)
        {
            return PI_DBY_180 * el.MajorSemiAxis_m * Math.Cos(lat_rad) / Math.Sqrt(1 - el.EccentricitySq * Math.Sin(lat_rad) * Math.Sin(lat_rad));
        }

        /// <summary>
        /// Calculates location of a point by base point and latitudal and longitudal projections on ellipsoid
        /// </summary>
        /// <param name="lat_rad">base point latitude, ragians</param>
        /// <param name="lon_rad">base point longitude, ragians</param>
        /// <param name="lat_offset_m">latitudal offset in meters</param>
        /// <param name="lon_offset_m">longitudal offset in meters</param>
        /// <param name="el">ellipsoid</param>
        /// <param name="e_lat_rad">new point latitude, ragians</param>
        /// <param name="e_lon_rad">new point longitude, ragians</param>
        public static void GeopointOffsetByDeltas(double lat_rad, double lon_rad, double lat_offset_m, double lon_offset_m,
            Ellipsoid el,
            out double e_lat_rad, out double e_lon_rad)
        {
            double m_per_deg_lat = Lat1DegLength(lat_rad, el);
            double m_per_deg_lon = Lon1DegLength(lat_rad, el);
            e_lat_rad = lat_rad - PI_DBY_180 * lat_offset_m / m_per_deg_lat;
            e_lon_rad = lon_rad - PI_DBY_180 * lon_offset_m / m_per_deg_lon;
        }

        /// <summary>
        /// Calculates location of a point by base point and latitudal and longitudal projections on WGS84 ellipsoid
        /// </summary>
        /// <param name="lat_rad">base point latitude, ragians</param>
        /// <param name="lon_rad">base point longitude, ragians</param>
        /// <param name="lat_offset_m">latitudal offset in meters</param>
        /// <param name="lon_offset_m">longitudal offset in meters</param>        
        /// <param name="e_lat_rad">new point latitude, ragians</param>
        /// <param name="e_lon_rad">new point longitude, ragians</param>
        public static void GeopointOffsetByDeltas_WGS84(double lat_rad, double lon_rad, double lat_offset_m, double lon_offset_m,
            out double e_lat_rad, out double e_lon_rad)
        {
            double m_per_deg_lat = 111132.92 - 559.82 * Math.Cos(2.0 * lat_rad) + 1.175 * Math.Cos(4.0 * lat_rad);
            double m_per_deg_lon = 111412.84 * Math.Cos(lat_rad) - 93.5 * Math.Cos(3.0 * lat_rad);
            e_lat_rad = lat_rad - PI_DBY_180 * lat_offset_m / m_per_deg_lat;
            e_lon_rad = lon_rad - PI_DBY_180 * lon_offset_m / m_per_deg_lon;
        }

        /// <summary>
        /// Calculates latitudal and longitudal projections of a line on ellipsoid between specified points
        /// </summary>
        /// <param name="sp_lat_rad">start point latitude, ragians</param>
        /// <param name="sp_lon_rad">start point longitude, ragians</param>
        /// <param name="ep_lat_rad">end point latitude, ragians</param>
        /// <param name="ep_lon_rad">end point longitude, ragians</param>
        /// <param name="el">ellipsoid</param>
        /// <param name="delta_lat_m">latitudal projection, meters</param>
        /// <param name="delta_lon_m">longitudal projection, meters</param>
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

        /// <summary>
        /// Calculates latitudal and longitudal projections of a line on WGS84 ellipsoid between specified points
        /// </summary>
        /// <param name="sp_lat_rad">start point latitude, ragians</param>
        /// <param name="sp_lon_rad">start point longitude, ragians</param>
        /// <param name="ep_lat_rad">end point latitude, ragians</param>
        /// <param name="ep_lon_rad">end point longitude, ragians</param>        
        /// <param name="delta_lat_m">latitudal projection, meters</param>
        /// <param name="delta_lon_m">longitudal projection, meters</param>
        public static void GetDeltasByGeopoints_WGS84(double sp_lat_rad, double sp_lon_rad, double ep_lat_rad, double ep_lon_rad,
            out double delta_lat_m, out double delta_lon_m)
        {
            double m_lat_rad = (sp_lat_rad + ep_lat_rad) / 2.0;
            double m_per_deg_lat = 111132.92 - 559.82 * Math.Cos(2.0 * m_lat_rad) + 1.175 * Math.Cos(4.0 * m_lat_rad);
            double m_per_deg_lon = 111412.84 * Math.Cos(m_lat_rad) - 93.5 * Math.Cos(3.0 * m_lat_rad);
            delta_lat_m = (sp_lat_rad - ep_lat_rad) * m_per_deg_lat * D180_DBY_PI;
            delta_lon_m = (sp_lon_rad - ep_lon_rad) * m_per_deg_lon * D180_DBY_PI;
        }

        #endregion
        
        #region Haversine equations

        /// <summary>
        /// Solves inverse geodetic problem according to haversine equation
        /// </summary>
        /// <param name="sp_lat_rad">start point latitude, radians</param>
        /// <param name="sp_lon_rad">start point longitude, radians</param>
        /// <param name="ep_lat_rad">end point latitude, radians</param>
        /// <param name="ep_lon_rad">end point longitude, radians</param>
        /// <param name="e_radius_m">Earth radius</param>
        /// <returns>distance between specified points</returns>
        public static double HaversineInverse(double sp_lat_rad, double sp_lon_rad, double ep_lat_rad, double ep_lon_rad, double e_radius_m)
        {
            double dLat = ep_lat_rad - sp_lat_rad;
            double dLon = ep_lon_rad - sp_lon_rad;
            double a = Math.Pow(Math.Sin(dLat / 2), 2) +
                      Math.Cos(sp_lat_rad) * Math.Cos(ep_lat_rad) *
                      Math.Sin(dLon / 2) * Math.Sin(dLon / 2);
            return e_radius_m * 2 * Math.Atan2(Math.Sqrt(a), Math.Sqrt(1 - a));
        }

        /// <summary>
        /// Solves direct geodetic problem according to haversine equation
        /// </summary>
        /// <param name="sp_lat_rad">start point latitude, radians</param>
        /// <param name="sp_lon_rad">start point longitude, radians</param>
        /// <param name="dst_m">distance on sphere, meters</param>
        /// <param name="fwd_az_rad">forward azimuth, radians</param>
        /// <param name="e_radius_m">earth radius, meters</param>
        /// <param name="ep_lat_rad">end point latitude, radians</param>
        /// <param name="ep_lon_rad">end point longitude, radians</param>
        public static void HaversineDirect(double sp_lat_rad, double sp_lon_rad, double dst_m, double fwd_az_rad, double e_radius_m,
            out double ep_lat_rad, out double ep_lon_rad)
        {
            double delta = dst_m / e_radius_m;
            ep_lat_rad = Math.Asin(Math.Sin(sp_lat_rad) * Math.Cos(delta) + Math.Cos(sp_lat_rad) * Math.Sin(delta) * Math.Cos(fwd_az_rad));
            ep_lon_rad = Wrap2PI(3 * Math.PI + (sp_lon_rad + Math.Atan2(Math.Sin(fwd_az_rad) * Math.Sin(delta) * Math.Cos(sp_lat_rad),
                Math.Cos(delta) - Math.Sin(sp_lat_rad) * Math.Sin(ep_lat_rad)))) - Math.PI;
        }

        /// <summary>
        /// Calculates initial bearing (forward azimuth) to a point
        /// </summary>
        /// <param name="sp_lat_rad">start point latitude, radians</param>
        /// <param name="sp_lon_rad">start point longitude, radians</param>
        /// <param name="ep_lat_rad">end point latitude, radians</param>
        /// <param name="ep_lon_rad">end point longitude, radians</param>
        /// <returns>initial bearing from start to end point</returns>
        public static double HaversineInitialBearing(double sp_lat_rad, double sp_lon_rad, double ep_lat_rad, double ep_lon_rad)
        {
            double y = Math.Sin(ep_lon_rad - sp_lon_rad) * Math.Cos(ep_lat_rad);
            double x = Math.Cos(sp_lat_rad) * Math.Sin(ep_lat_rad) - Math.Sin(sp_lat_rad) * Math.Cos(ep_lat_rad) * Math.Cos(ep_lon_rad - sp_lon_rad);
            return Wrap2PI(Math.PI + Math.Atan2(y, x));
        }

        /// <summary>
        /// Calculates final bearing (reverse azimuth) to a point
        /// </summary>
        /// <param name="sp_lat_rad">start point latitude, radians</param>
        /// <param name="sp_lon_rad">start point longitude, radians</param>
        /// <param name="ep_lat_rad">end point latitude, radians</param>
        /// <param name="ep_lon_rad">end point longitude, radians</param>
        /// <returns>final bearing from end to start point</returns>
        public static double HaversineFinalBearing(double sp_lat_rad, double sp_lon_rad, double ep_lat_rad, double ep_lon_rad)
        {
            return Wrap2PI(HaversineInitialBearing(ep_lat_rad, ep_lon_rad, sp_lat_rad, sp_lon_rad) + Math.PI);
        }

        #endregion

        #region Vincenty equations

        /// <summary>
        /// Solves inverse geodetic problem according to vincenty equation
        /// </summary>
        /// <param name="sp_lat_rad">start point latitude, radians</param>
        /// <param name="sp_lon_rad">start point longitude, radians</param>
        /// <param name="ep_lat_rad">end point latitude, radians</param>
        /// <param name="ep_lon_rad">end point longitude, radians</param>
        /// <param name="el">ellipsoid</param>
        /// <param name="eps">epsilon</param>
        /// <param name="it_limit">iterations limit</param>
        /// <param name="dst_m">distance, meters</param>
        /// <param name="fwd_az_rad">forward azimuth, radians</param>
        /// <param name="rev_az_rad">reverse azimuth, radians</param>
        /// <param name="its">number of iterations taken</param>
        /// <returns>true if algorithm converges, false otherwise</returns>
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


            double u_sq = cos_sq_alpha * (el.MajorSemiAxis_m * el.MajorSemiAxis_m - el.MinorSemiAxis_m * el.MinorSemiAxis_m) / (el.MinorSemiAxis_m * el.MinorSemiAxis_m);
            double a_ = 1 + u_sq / 16384 * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)));
            double b_ = u_sq / 1024 * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)));
            double delta_sigma = b_ * sin_sigma * (cos_2_sigma_m + b_ / 4 * (cos_sigma * (-1 + 2 * cos_2_sigma_m * cos_2_sigma_m) -
                b_ / 6 * cos_2_sigma_m * (-3 + 4 * sin_sigma * sin_sigma) * (-3 + 4 * cos_2_sigma_m * cos_2_sigma_m)));


            dst_m = el.MinorSemiAxis_m * a_ * (sigma - delta_sigma);

            fwd_az_rad = Math.Atan2(cos_u_2 * sin_lambda, cos_u_1 * sin_u_2 - sin_u_1 * cos_u_2 * cos_lambda);
            rev_az_rad = Math.Atan2(cos_u_1 * sin_lambda, -sin_u_1 * cos_u_2 + cos_u_1 * sin_u_2 * cos_lambda);

            while (fwd_az_rad < 0)
                fwd_az_rad += PI2;

            while (rev_az_rad < 0)
                rev_az_rad += PI2;

            return (its < it_limit) && (it_check < Math.PI);
        }

        /// <summary>
        /// Solves direct geodetic problem according to vincenty equation
        /// </summary>
        /// <param name="sp_lat_rad">start point latitude, rad</param>
        /// <param name="sp_lon_rad">start point longitude, rad</param>
        /// <param name="fwd_az_rad">forward azimuth, rad</param>
        /// <param name="dst_m">distance, m</param>
        /// <param name="el">ellipsoid</param>
        /// <param name="eps">epsilon</param>
        /// <param name="it_limit">iterations limit</param>
        /// <param name="ep_lat_rad">end point latitude, rad</param>
        /// <param name="ep_lon_rad">end point longitude, rad</param>
        /// <param name="rev_az_rad">reverse azimuth, rad</param>
        /// <param name="its">number of iterations taken</param>
        /// <returns>true is algorithm converges, false - otherwise</returns>
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

        #endregion
        
        #region TOA & TDOA solvers
       
        /// <summary>
        /// TOA problem residual function
        /// </summary>
        /// <param name="basePoints">base points with known locations and distances to them</param>
        /// <param name="x">current x coordinate</param>
        /// <param name="y">current y coordinate</param>
        /// <param name="z">current z coordinate</param>
        /// <returns>value of residual function in specified location</returns>
        public static double Eps_TOA3D(TOABasePoint[] basePoints, double x, double y, double z)
        {
            double result = 0;
            double eps = 0;

            for (int i = 0; i < basePoints.Length; i++)
            {
                eps =  Math.Sqrt((basePoints[i].X - x) * (basePoints[i].X - x) +
                                 (basePoints[i].Y - y) * (basePoints[i].Y - y) +
                                 (basePoints[i].Z - z) * (basePoints[i].Z - z)) - basePoints[i].D;
                result += eps * eps;
            }

            return result;
        }

        /// <summary>
        /// TDOA problem residual function
        /// </summary>
        /// <param name="baseLines">base lines, each represented by two base points with known locations and times of arrival</param>
        /// <param name="x">current x coordinate</param>
        /// <param name="y">current y coordinate</param>
        /// <param name="z">current z coordinate</param>
        /// <returns>value of residual function in specified location</returns>
        public static double Eps_TDOA3D(TDOABaseline[] baseLines, double x, double y, double z)
        {
            double result = 0;
            double eps;

            for (int i = 0; i < baseLines.Length; i++)
            {
                eps = Math.Sqrt((baseLines[i].X1 - x) * (baseLines[i].X1 - x) +
                                (baseLines[i].Y1 - y) * (baseLines[i].Y1 - y) +
                                (baseLines[i].Z1 - z) * (baseLines[i].Z1 - z)) -
                      Math.Sqrt((baseLines[i].X2 - x) * (baseLines[i].X2 - x) +
                                (baseLines[i].Y2 - y) * (baseLines[i].Y2 - y) +
                                (baseLines[i].Z2 - z) * (baseLines[i].Z2 - z)) - baseLines[i].PRD;
                result += eps * eps;
            }

            return result;
        }

        /// <summary>
        /// Finds minimum of specifief residual function by Nelder-Mead (simplex) method (2D - x and y)
        /// </summary>
        /// <typeparam name="T">Can be TOABasePoint or TDOABaseLine</typeparam>
        /// <param name="eps">residual function</param>
        /// <param name="baseElements">a set of base elements of type T</param>
        /// <param name="xPrev">previous solution x coordinate</param>
        /// <param name="yPrev">previous solution y coordinate</param>
        /// <param name="z">z coordinate (depth), suppose to be a constant (known by direct measurement)</param>
        /// <param name="maxIterations">iterations limit</param>
        /// <param name="precisionThreshold">precision threshold</param>
        /// <param name="simplexSize">initial size of the simplex</param>
        /// <param name="xBest">x coordinate of minimum</param>
        /// <param name="yBest">y coordinate of minimum</param>
        /// <param name="radialError">radial error, m</param>
        /// <param name="itCnt">number of iterations taken</param>
        public static void NLM2D_Solve<T>(Func<T[], double, double, double, double> eps,
                                          T[] baseElements, double xPrev, double yPrev, double z,
                                          int maxIterations, double precisionThreshold, double simplexSize,
                                          out double xBest, out double yBest, out double radialError, out int itCnt)
        {
            #region Nelder-Mead 2D optimization

            bool isFinished = false;
            
            double tmp, tmp1;
            double[] xix = new double[3];
            double[] xiy = new double[3];
            double[] fxi = new double[3];
            double fl, fg, fh, fr, fe, fs;
            double xcx, xcy, xrx, xry, xex, xey, xsx, xsy;

            itCnt = 0;

            xix[0] = xPrev;
            xiy[0] = yPrev;
            xix[1] = xix[0] + simplexSize;
            xiy[1] = xiy[0] + simplexSize;
            xix[2] = xix[0] - simplexSize / 2;
            xiy[2] = xiy[0] + simplexSize / 2;

            while (!isFinished)
            {
                fxi[0] = eps(baseElements, xix[0], xiy[0], z);
                fxi[1] = eps(baseElements, xix[1], xiy[1], z);
                fxi[2] = eps(baseElements, xix[2], xiy[2], z);

                if (fxi[0] > fxi[1])
                {
                    tmp = fxi[0]; fxi[0] = fxi[1]; fxi[1] = tmp;
                    tmp = xix[0]; xix[0] = xix[1]; xix[1] = tmp;
                    tmp = xiy[0]; xiy[0] = xiy[1]; xiy[1] = tmp;
                }

                if (fxi[0] > fxi[2])
                {
                    tmp = fxi[0]; fxi[0] = fxi[2]; fxi[2] = tmp;
                    tmp = xix[0]; xix[0] = xix[2]; xix[2] = tmp;
                    tmp = xiy[0]; xiy[0] = xiy[2]; xiy[2] = tmp;
                }

                if (fxi[1] > fxi[2])
                {
                    tmp = fxi[1]; fxi[1] = fxi[2]; fxi[2] = tmp;
                    tmp = xix[1]; xix[1] = xix[2]; xix[2] = tmp;
                    tmp = xiy[1]; xiy[1] = xiy[2]; xiy[2] = tmp;
                }

                fl = fxi[0];
                fg = fxi[1];
                fh = fxi[2];

                xcx = (xix[0] + xix[1]) / 2.0f;
                xcy = (xiy[0] + xiy[1]) / 2.0f;

                xrx = (1.0f + NLM_A) * xcx - NLM_A * xix[2];
                xry = (1.0f + NLM_A) * xcy - NLM_A * xiy[2];

                fr = eps(baseElements, xrx, xry, z);

                if (fr < fl)
                {
                    xex = (1.0f - NLM_G) * xcx + NLM_G * xrx;
                    xey = (1.0f - NLM_G) * xcy + NLM_G * xry;

                    fe = eps(baseElements, xex, xey, z);

                    if (fe < fr)
                    {
                        xix[2] = xex;
                        xiy[2] = xey;
                    }
                    else
                    {
                        xix[2] = xrx;
                        xiy[2] = xry;
                    }
                }
                else
                {
                    if ((fr > fl) && (fr < fg))
                    {
                        xix[2] = xrx;
                        xiy[2] = xry;
                    }
                    else
                    {
                        if ((fr > fg) && (fr < fh))
                        {
                            tmp = xix[2]; xix[2] = xrx; xrx = tmp;
                            tmp = xiy[2]; xiy[2] = xry; xry = tmp;
                            tmp = fxi[2]; fxi[2] = fr; fr = tmp;
                        }
                        else
                        {
                            if (fh < fr)
                            {
                                //
                            }
                        }

                        xsx = NLM_B * xix[2] + (1.0f - NLM_B) * xcx;
                        xsy = NLM_B * xiy[2] + (1.0f - NLM_B) * xcy;
                        fs = eps(baseElements, xsx, xsy, z);

                        if (fs < fh)
                        {
                            xix[2] = xsx;
                            xiy[2] = xsy;
                        }
                        else
                        {
                            xix[1] = (xix[1] - xix[0]) / 2.0f;
                            xiy[1] = (xiy[1] - xiy[0]) / 2.0f;
                            xix[2] = (xix[2] - xix[0]) / 2.0f;
                            xiy[2] = (xiy[2] - xiy[0]) / 2.0f;
                        }
                    }
                }

                tmp = (fxi[0] + fxi[1] + fxi[2]) / 3.0f;
                tmp1 = ((fxi[0] - tmp) * (fxi[0] - tmp) +
                        (fxi[1] - tmp) * (fxi[1] - tmp) +
                        (fxi[2] - tmp) * (fxi[2] - tmp)) / 3.0f;

                isFinished = (++itCnt < maxIterations) && (Math.Sqrt(tmp1) <= precisionThreshold);
            }

            #endregion

            xBest = xix[0];
            yBest = xiy[0];
            radialError = Math.Sqrt(eps(baseElements, xBest, yBest, z));
        }

        /// <summary>
        /// Finds minimum of specifief residual function by Nelder-Mead (simplex) method (3D - x, y and z)
        /// </summary>
        /// <typeparam name="T">Can be TOABasePoint or TDOABaseLine</typeparam>
        /// <param name="eps">residual function</param>
        /// <param name="baseElements">a set of base elements of type T</param>
        /// <param name="xPrev">previous solution x coordinate</param>
        /// <param name="yPrev">previous solution y coordinate</param>
        /// <param name="zPrev">previous solution z coordinate (depth)</param>
        /// <param name="maxIterations">iterations limit</param>
        /// <param name="precisionThreshold">precision threshold</param>
        /// <param name="simplexSize">initial size of the simplex</param>
        /// <param name="xBest">x coordinate of minimum</param>
        /// <param name="yBest">y coordinate of minimum</param>
        /// <param name="zBest">z coordinate of minimum</param>
        /// <param name="radialError">radial error, m</param>
        /// <param name="itCnt">number of iterations taken</param>
        public static void NLM3D_Solve<T>(Func<T[], double, double, double, double> eps,
                                          T[] baseElements, double xPrev, double yPrev, double zPrev,
                                          int maxIterations, double precisionThreshold, double simplexSize,
                                          out double xBest, out double yBest, out double zBest, out double radialError, out int itCnt)
        {
            #region Nelder-Mead 3D optimization

            int v_num = 4; // vertices number is 4 for 3D-task

            double[] xix = new double[v_num]; 
            double[] xiy = new double[v_num];
            double[] xia = new double[v_num];
            double[] fxi = new double[v_num];

            xix[0] = xPrev;
            xiy[0] = yPrev;
            xia[0] = zPrev;
            xix[1] = xix[0] + simplexSize;
            xiy[1] = xiy[0] + simplexSize;
            xia[1] = xia[0] - simplexSize;
            xix[2] = xix[0] - simplexSize / 2;
            xiy[2] = xiy[0] + simplexSize / 2;
            xia[2] = xia[0] + simplexSize / 2;
            xix[3] = xix[0] + simplexSize / 2;
            xiy[3] = xiy[0] - simplexSize / 2;
            xia[3] = xia[0] - simplexSize / 2;

            bool isFinished = false;
            double sigma, mean, tmp;
            double xrx, xry, xra, x0x, x0y, x0a, xex, xey, xea, xcx, xcy, xca, fc;
            double fr, fe;
            itCnt = 0;

            while (!isFinished)
            {
                for (int i = 0; i < v_num; i++)
                {
                    fxi[i] = eps(baseElements, xix[i], xiy[i], xia[i]);
                }

                mean = (fxi[0] + fxi[1] + fxi[2] + fxi[3]) / v_num;
                sigma = Math.Sqrt(((fxi[0] - mean) * (fxi[0] - mean) +
                                   (fxi[1] - mean) * (fxi[1] - mean) +
                                   (fxi[2] - mean) * (fxi[2] - mean) +
                                   (fxi[3] - mean) * (fxi[3] - mean)) / v_num);

                if ((++itCnt > maxIterations) || (sigma < precisionThreshold))
                {
                    isFinished = true;
                }
                else
                {
                    // Sort vertices
                    if (fxi[0] > fxi[1])
                    {
                        tmp = fxi[0]; fxi[0] = fxi[1]; fxi[1] = tmp;
                        tmp = xix[0]; xix[0] = xix[1]; xix[1] = tmp;
                        tmp = xiy[0]; xiy[0] = xiy[1]; xiy[1] = tmp;
                        tmp = xia[0]; xia[0] = xia[1]; xia[1] = tmp;
                    }
                    if (fxi[0] > fxi[2])
                    {
                        tmp = fxi[0]; fxi[0] = fxi[2]; fxi[2] = tmp;
                        tmp = xix[0]; xix[0] = xix[2]; xix[2] = tmp;
                        tmp = xiy[0]; xiy[0] = xiy[2]; xiy[2] = tmp;
                        tmp = xia[0]; xia[0] = xia[2]; xia[2] = tmp;
                    }
                    if (fxi[0] > fxi[3])
                    {
                        tmp = fxi[0]; fxi[0] = fxi[3]; fxi[3] = tmp;
                        tmp = xix[0]; xix[0] = xix[3]; xix[3] = tmp;
                        tmp = xiy[0]; xiy[0] = xiy[3]; xiy[3] = tmp;
                        tmp = xia[0]; xia[0] = xia[3]; xia[3] = tmp;
                    }
                    if (fxi[1] > fxi[2])
                    {
                        tmp = fxi[1]; fxi[1] = fxi[2]; fxi[2] = tmp;
                        tmp = xix[1]; xix[1] = xix[2]; xix[2] = tmp;
                        tmp = xiy[1]; xiy[1] = xiy[2]; xiy[2] = tmp;
                        tmp = xia[1]; xia[1] = xia[2]; xia[2] = tmp;
                    }
                    if (fxi[1] > fxi[3])
                    {
                        tmp = fxi[1]; fxi[1] = fxi[3]; fxi[3] = tmp;
                        tmp = xix[1]; xix[1] = xix[3]; xix[3] = tmp;
                        tmp = xiy[1]; xiy[1] = xiy[3]; xiy[3] = tmp;
                        tmp = xia[1]; xia[1] = xia[3]; xia[3] = tmp;
                    }
                    if (fxi[2] > fxi[3])
                    {
                        tmp = fxi[2]; fxi[2] = fxi[3]; fxi[3] = tmp;
                        tmp = xix[2]; xix[2] = xix[3]; xix[3] = tmp;
                        tmp = xiy[2]; xiy[2] = xiy[3]; xiy[3] = tmp;
                        tmp = xia[2]; xia[2] = xia[3]; xia[3] = tmp;
                    }

                    // (2) Calculate x0, the centroid of all points except xn+1

                    x0x = (xix[0] + xix[1] + xix[2]) / 3;
                    x0y = (xiy[0] + xiy[1] + xiy[2]) / 3;
                    x0a = (xia[0] + xia[1] + xia[2]) / 3;

                    // (3) reflect the xh (xi(4)) point to xr
                    xrx = x0x + NLM_A * (x0x - xix[3]);
                    xry = x0y + NLM_A * (x0y - xiy[3]);
                    xra = x0a + NLM_A * (x0a - xia[3]);

                    // function value in the xr point                    
                    fr = eps(baseElements, xrx, xry, xra);

                    // if fx1 <= fxr <= fxn replace worst point xn+1 with xr
                    if ((fr >= fxi[0]) && (fxi[2] >= fr))
                    {
                        xix[3] = xrx;
                        xiy[3] = xry;
                        xia[3] = xra;
                        fxi[3] = fr;
                    }
                    else
                    {
                        // (4) expansion
                        if (fr < fxi[0])
                        {
                            xex = x0x + NLM_G * (xrx - x0x);
                            xey = x0x + NLM_G * (xry - x0x);
                            xea = x0x + NLM_G * (xra - x0x);
                            fe = eps(baseElements, xex, xey, xea);

                            if (fe < fr)
                            {
                                xix[3] = xex;
                                xiy[3] = xey;
                                xia[3] = xea;
                                fxi[3] = fe;
                                // % to step 1
                            }
                            else
                            {
                                xix[3] = xrx;
                                xiy[3] = xry;
                                xia[3] = xra;
                                fxi[3] = fr;
                                // % to step 1
                            }
                        }
                        else
                        {
                            xcx = x0x + NLM_R * (xix[3] - x0x);
                            xcy = x0y + NLM_R * (xiy[3] - x0y);
                            xca = x0a + NLM_R * (xia[3] - x0a);
                            fc = eps(baseElements, xcx, xcy, xca);

                            if (fc < fxi[3])
                            {
                                xix[3] = xcx;
                                xiy[3] = xcy;
                                xia[3] = xca;
                                fxi[3] = fc;
                                // % to step 1
                            }
                            else
                            {
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

            #endregion

            xBest = xix[0];
            yBest = xiy[0];
            zBest = xia[0];
            radialError = Math.Sqrt(eps(baseElements, xBest, yBest, zBest));
        }

        /// <summary>
        /// Solves navigation problem: finds target location by TOA (by base points with known location and distances to them)
        /// with Nelder-Mead (simplex) algorithm. 2D problem (x and y are variables, z is supposed to be known by direct measurement)
        /// </summary>
        /// <param name="basePoints">Base points with known locations and distances to them from target</param>
        /// <param name="xPrev">previous solution x coordinate</param>
        /// <param name="yPrev">previous solution x coordinate</param>
        /// <param name="z">z coordinate</param>
        /// <param name="maxIterations">iterations limit</param>
        /// <param name="precisionThreshold">precision threshold</param>
        /// <param name="simplexSize">initial size of simplex</param>
        /// <param name="xBest">x coordinate of target</param>
        /// <param name="yBest">y coordinate of target</param>
        /// <param name="radialError">radial error</param>
        /// <param name="itCnt">number of iterations taken</param>                                                                                
        public static void TOA_NLM2D_Solve(TOABasePoint[] basePoints,
                                           double xPrev, double yPrev, double z,
                                           int maxIterations, double precisionThreshold, double simplexSize,
                                           out double xBest, out double yBest, out double radialError, out int itCnt)
        {
            NLM2D_Solve<TOABasePoint>(Eps_TOA3D,
                                      basePoints, xPrev, yPrev, z,
                                      maxIterations, precisionThreshold, simplexSize,
                                      out xBest, out yBest, out radialError, out itCnt);
        }

        /// <summary>
        /// Solves navigation problem: finds target location by TDOA (by base points with known location and times of arrival)
        /// with Nelder-Mead (simplex) algorithm. 2D problem (x and y are variables, z is supposed to be known by direct measurement)
        /// </summary>
        /// <param name="baseLines">A set of base lines represented each by two base points with known locations and time of arrival</param>
        /// <param name="xPrev">previous solution x coordinate</param>
        /// <param name="yPrev">previous solution x coordinate</param>
        /// <param name="z">z coordinate</param>
        /// <param name="maxIterations">iterations limit</param>
        /// <param name="precisionThreshold">precision threshold</param>
        /// <param name="simplexSize">initial size of simplex</param>
        /// <param name="xBest">x coordinate of target</param>
        /// <param name="yBest">y coordinate of target</param>
        /// <param name="radialError">radial error</param>
        /// <param name="itCnt">number of iterations taken</param>   
        public static void TDOA_NLM2D_Solve(TDOABaseline[] baseLines,
                                            double xPrev, double yPrev, double z,
                                            int maxIterations, double precisionThreshold, double simplexSize,
                                            out double xBest, out double yBest, out double radialError, out int itCnt)
        {
            NLM2D_Solve<TDOABaseline>(Eps_TDOA3D,
                                      baseLines, xPrev, yPrev, z,
                                      maxIterations, precisionThreshold, simplexSize,
                                      out xBest, out yBest, out radialError, out itCnt);
        }

        /// <summary>
        /// Solves navigation problem: finds target location by TOA (by base points with known location and distances to them)
        /// with Nelder-Mead (simplex) algorithm. 3D problem (x, y and z are variables)
        /// </summary>
        /// <param name="basePoints">Base points with known locations and distances to them from target</param>
        /// <param name="xPrev">previous solution x coordinate</param>
        /// <param name="yPrev">previous solution x coordinate</param>
        /// <param name="zPrev">previous solution z coordinate</param>
        /// <param name="maxIterations">iterations limit</param>
        /// <param name="precisionThreshold">precision threshold</param>
        /// <param name="simplexSize">initial size of simplex</param>
        /// <param name="xBest">x coordinate of target</param>
        /// <param name="yBest">y coordinate of target</param>
        /// <param name="zBest">z coordinate of target</param>
        /// <param name="radialError">radial error</param>
        /// <param name="itCnt">number of iterations taken</param>         
        public static void TOA_NLM3D_Solve(TOABasePoint[] basePoints,
                                           double xPrev, double yPrev, double zPrev,
                                           int maxIterations, double precisionThreshold, double simplexSize,
                                           out double xBest, out double yBest, out double zBest, out double radialError, out int itCnt)
        {
            NLM3D_Solve<TOABasePoint>(Eps_TOA3D,
                                      basePoints, xPrev, yPrev, zPrev,
                                      maxIterations, precisionThreshold, simplexSize,
                                      out xBest, out yBest, out zBest, out radialError, out itCnt);
        }

        /// <summary>
        /// Solves navigation problem: finds target location by TDOA (by base points with known location and times of arrival)
        /// with Nelder-Mead (simplex) algorithm. 3D problem (x, y and z are variables)
        /// </summary>
        /// <param name="baseLines">A set of base lines represented each by two base points with known locations and time of arrival</param>
        /// <param name="xPrev">previous solution x coordinate</param>
        /// <param name="yPrev">previous solution x coordinate</param>
        /// <param name="zPrev">previous solution z coordinate</param>
        /// <param name="maxIterations">iterations limit</param>
        /// <param name="precisionThreshold">precision threshold</param>
        /// <param name="simplexSize">initial size of simplex</param>
        /// <param name="xBest">x coordinate of target</param>
        /// <param name="yBest">y coordinate of target</param>
        /// <param name="zBest">z coordinate of target</param>
        /// <param name="radialError">radial error</param>
        /// <param name="itCnt">number of iterations taken</param>   
        public static void TDOA_NLM3D_Solve(TDOABaseline[] baseLines,
                                            double xPrev, double yPrev, double zPrev,
                                            int maxIterations, double precisionThreshold, double simplexSize,
                                            out double xBest, out double yBest, out double zBest, out double radialError, out int itCnt)
        {
            NLM3D_Solve<TDOABaseline>(Eps_TDOA3D,
                                      baseLines, xPrev, yPrev, zPrev,
                                      maxIterations, precisionThreshold, simplexSize,
                                      out xBest, out yBest, out zBest, out radialError, out itCnt);
        }
        
        /// <summary>
        /// Finds the index of the TOABasePoint with smallest D property
        /// </summary>
        /// <param name="basePoints">Source base points</param>
        /// <returns>Index of an item with smallest D value</returns>
        public static int GetNearestItemIndex(TOABasePoint[] basePoints)
        {
            int nrstIdx = 0;
            double minDst = double.MaxValue;

            for (int i = 0; i < basePoints.Length; i++)
            {
                if (basePoints[i].D < minDst)
                {
                    minDst = basePoints[i].D;
                    nrstIdx = i;
                }
            }

            return nrstIdx;
        }

        /// <summary>
        /// Searches minimum of residual function as a point on a circle
        /// </summary>
        /// <param name="basePoints">Base points with known locations and distances to target</param>
        /// <param name="anchorx">x coordinate of the center of circle (e.g. nearest base point)</param>
        /// <param name="anchory">x coordinate of the center of circle (e.g. nearest base point)</param>
        /// <param name="radius">radius of circle (distance from base point to target) in meters</param>
        /// <param name="z">z coordinate of the target (supposed to be known by direct measurement)</param>
        /// <param name="arcMidRad">start angle on arc, in radians</param>
        /// <param name="arcAngleRad">angular dimension of arc to search in it, in radians</param>
        /// <param name="steps">number of steps in which arcAngleRag will be divided</param>
        /// <returns>angle, where the value of residual function is minimal</returns>
        public static double TOA_CirclesIntersection_Solve(TOABasePoint[] basePoints,                                                           
                                                           double anchorx, double anchory, double radius,
                                                           double z,
                                                           double arcMidRad, double arcAngleRad, int steps)
        {
            double a = arcMidRad - arcAngleRad / 2;
            double aEnd = arcMidRad + arcAngleRad / 2;
            double aBest = a;
            double stepRad = arcAngleRad / steps;
            double epsBest = double.MaxValue;
            double x, y, eps;

            while (a < aEnd)
            {
                x = anchorx + radius * Math.Cos(a);
                y = anchory + radius * Math.Sin(a);
                eps = Eps_TOA3D(basePoints, x, y, z);

                if (eps < epsBest)
                {
                    epsBest = eps;
                    aBest = a;
                }

                a += stepRad;
            }

            return aBest;
        }

        /// <summary>
        /// Solves TOA problem by searching minimum of the residual function on a circle, whose center is a base point and radius is measured distance from 
        /// this base point to target
        /// </summary>
        /// <param name="basePoints">Base points with known location and distances to target</param>
        /// <param name="z">z coordinate of the target, supposed to be known by direct measurement, in meters</param>
        /// <param name="endArcAngleRad">minimal angular size of an arc, when algorithm should stop, in radians</param>
        /// <param name="steps">number of part to divide the arc</param>
        /// <param name="arcAngleDecreaseFactor">multiplier of arc angular size on each step</param>
        /// <param name="xBest">x coordinate of the target, in meters</param>
        /// <param name="yBest">y coordinate of the target, in meters</param>
        /// <param name="radialError">radial error, in meters</param>
        public static void TOA_Circles1D_Solve(TOABasePoint[] basePoints,
                                               double z,
                                               double endArcAngleRad, int steps,
                                               double arcAngleDecreaseFactor,
                                               out double xBest, out double yBest, out double radialError)                                                              
        {            
            int nrstIdx = GetNearestItemIndex(basePoints);
            double dZ = Math.Abs(basePoints[nrstIdx].Z - z);
            double radius = basePoints[nrstIdx].D < dZ ? 0.0 : Math.Sqrt(basePoints[nrstIdx].D * basePoints[nrstIdx].D - dZ * dZ);
            double anchorx = basePoints[nrstIdx].X;
            double anchory = basePoints[nrstIdx].Y;
            double arcAngle = Math.PI * 2;
            double alpha = 0;

            while (arcAngle > endArcAngleRad)
            {
                alpha = TOA_CirclesIntersection_Solve(basePoints,
                                                      anchorx, anchory, radius,
                                                      z,
                                                      alpha, arcAngle, steps);

                arcAngle *= arcAngleDecreaseFactor;                
            }
                        
            xBest = anchorx + radius * Math.Cos(alpha);
            yBest = anchory + radius * Math.Sin(alpha);
            radialError = Math.Sqrt(Eps_TOA3D(basePoints, xBest, yBest, z));
        }
                
        #endregion

        #endregion
    }
}