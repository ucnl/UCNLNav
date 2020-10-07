using System;
using System.Collections.Generic;

namespace UCNLNav
{
    public static class Navigation
    {
        
        #region Methods

        /// <summary>
        /// Calculates the centroid of a cloud of points      
        /// </summary>
        /// <param name="points">Cloud of points</param>
        /// <returns>The centroid of the specified cloud of points</returns>
        public static GeoPoint GetPointsCentroid2D(IEnumerable<GeoPoint> points)
        {            
            double cLat = 0.0;
            double cLon = 0.0;
            int cnt = 0;

            foreach (var point in points)
            {
                cLat += point.Latitude;
                cLon += point.Longitude;
                cnt++;
            }

            if (cnt > 0)
            {
                return new GeoPoint(cLat / cnt, cLon / cnt);
            }
            else
            {
                return new GeoPoint();
            }
        }

        public static MPoint GetPointsCentroid2D(IEnumerable<MPoint> points)
        {
            double cx = 0.0;
            double cy = 0.0;
            int cnt = 0;

            foreach (var point in points)
            {
                cx += point.X;
                cy += point.Y;
                cnt++;
            }

            if (cnt > 0)
            {
                return new MPoint(cx / cnt, cy / cnt);
            }
            else
            {
                return new MPoint();
            }
        }

        public static MPoint3D GetPointsCentroid3D(IEnumerable<MPoint3D> points)
        {
            double cx = 0.0;
            double cy = 0.0;
            double cz = 0.0;
            int cnt = 0;

            foreach (var point in points)
            {
                cx += point.X;
                cy += point.Y;
                cz += point.Z;
                cnt++;
            }

            if (cnt > 0)
            {
                return new MPoint3D(cx / cnt, cy / cnt, cz / cnt);
            }
            else
            {
                return new MPoint3D(0, 0, 0);
            }
        }        

        public static void GetPointsSTD2D(IEnumerable<MPoint> points, out double sigmax, out double sigmay)
        {
            sigmax = 0.0;
            sigmay = 0.0;
            int cnt = 0;

            var centroid = GetPointsCentroid2D(points);

            foreach (var point in points)
            {
                sigmax += Math.Pow(point.X - centroid.X, 2);
                sigmay += Math.Pow(point.Y - centroid.Y, 2);
                cnt++;
            }

            if (cnt > 0)
            {
                sigmax /= cnt;
                sigmay /= cnt;
            }            
        }

        public static void GetPointsSTD3D(IEnumerable<MPoint3D> points, out double sigmax, out double sigmay, out double sigmaz)
        {
            sigmax = 0.0;
            sigmay = 0.0;
            sigmaz = 0.0;
            int cnt = 0;

            var centroid = GetPointsCentroid3D(points);

            foreach (var point in points)
            {
                sigmax += Math.Pow(point.X - centroid.X, 2);
                sigmay += Math.Pow(point.Y - centroid.Y, 2);
                sigmaz += Math.Pow(point.Z - centroid.Z, 2);
                cnt++;
            }

            if (cnt > 0)
            {
                sigmax /= cnt;
                sigmay /= cnt;
                sigmaz /= cnt;
            }
        }


        /// <summary>
        /// Converts geographic coordinates to local cartesian coordinates
        /// </summary>
        /// <param name="points">Geographic points</param>
        /// <param name="el">Reference ellipsoid</param>
        /// <returns>List of MPoints in cartesian coordinates system, which center is in the centroid on specifief points</returns>
        public static List<MPoint> GCSToLCS(IEnumerable<GeoPoint> points, Ellipsoid el)
        {
            List<MPoint> result = new List<MPoint>();
            var centroid = GetPointsCentroid2D(points);            

            double cLat = Algorithms.Deg2Rad(centroid.Latitude);
            double cLon = Algorithms.Deg2Rad(centroid.Longitude);
            double d_lat_m, d_lon_m;

            foreach (var point in points)
            {
                Algorithms.GetDeltasByGeopoints(cLat, cLon, 
                    Algorithms.Deg2Rad(point.Latitude), Algorithms.Deg2Rad(point.Longitude),
                    el, 
                    out d_lat_m, out d_lon_m);

                result.Add(new MPoint(d_lon_m, d_lat_m));
            }

            return result;
        }

        public static List<GeoPoint> LCSToGCS(IEnumerable<MPoint> points, GeoPoint centroid, Ellipsoid el)
        {
            List<GeoPoint> result = new List<GeoPoint>();

            double cLat = Algorithms.Deg2Rad(centroid.Latitude);
            double cLon = Algorithms.Deg2Rad(centroid.Longitude);
            double eLat = 0;
            double eLon = 0;

            foreach (var point in points)
            {
                Algorithms.GeopointOffsetByDeltas(cLat, cLon, point.Y, point.X, el, out eLat, out eLon);
                result.Add(new GeoPoint(eLat, eLon));
            }

            return result;
        }


        /// <summary>
        /// Calculates the Distance Root Mean Squared (DRMS) of the specified
        /// standard errors.
        /// DRMS -> 65% probability
        /// 2DRMS -> 95% probability
        /// 3DRMS -> 98% probability
        /// </summary>
        /// <param name="sigmax">Standard error of estimates X coordinate</param>
        /// <param name="sigmay">Standard error of estimates Y coordinate</param>
        /// <returns>The Distance Root Mean Squared of specified standard errors</returns>
        public static double DRMS(double sigmax, double sigmay)
        {
            return Math.Sqrt(sigmax * sigmax + sigmay * sigmay);
        }

        /// <summary>
        /// Calculates the the Circular Error Probability - the radius of circle centered at the true position,
        /// containing the position estimate with probability of 50%
        /// </summary>
        /// <param name="sigmax">Standard error of estimates X coordinate</param>
        /// <param name="sigmay">Standard error of estimates Y coordinate</param>
        /// <returns>The Circular Error Probability of specified standard errors</returns>
        public static double CEP(double sigmax, double sigmay)
        {
            return 0.62 * sigmay + 0.56 * sigmax;            
        }

        /// <summary>
        /// Calculates Spherical Error Probability - the radius of the sphere
        /// centered at the true position, containing the position estimate in 3D
        /// with probability of 50%
        /// </summary>
        /// <param name="sigmax">Standard error of estimates X coordinate</param>
        /// <param name="sigmay">Standard error of estimates Y coordinate</param>
        /// <param name="sigmaz">Standard error of estimates Z coordinate</param>
        /// <returns>The Spherical Error Probability of specified standard errors</returns>
        public static double SEP(double sigmax, double sigmay, double sigmaz)
        {
            return 0.51 * (sigmay + sigmax + sigmaz);
        }

        /// <summary>
        /// Calculates the Mean Radial Spherical Error - the radius of sphere
        /// centered at the true position, containing the position estimate in 3D
        /// with probability of 61%
        /// </summary>
        /// <param name="sigmax">Standard error of estimates X coordinate</param>
        /// <param name="sigmay">Standard error of estimates Y coordinate</param>
        /// <param name="sigmaz">Standard error of estimates Z coordinate</param>
        /// <returns>The Mean Radial Spherical Error of specified standard errors</returns>
        public static double MRSE(double sigmax, double sigmay, double sigmaz)
        {
            return Math.Sqrt(sigmax * sigmax + sigmay * sigmay + sigmaz * sigmaz);
        }



               
        
        /// <summary>
        /// Converts array of GeoPoint3DD to an array of TOABasePoints, converts lat/lon to metric local CS whose center is in specified centroid
        /// </summary>
        /// <param name="bases">An array of base points with measured distances</param>
        /// <param name="centroid">Start of new local coordinate system</param>
        /// <param name="el">ellipsoid</param>
        /// <returns>An array of TOABasePoints</returns>
        public static TOABasePoint[] ConvertToLCS(GeoPoint3DD[] bases, GeoPoint centroid, Ellipsoid el)
        {
            TOABasePoint[] result = new TOABasePoint[bases.Length];

            double cLat = Algorithms.Deg2Rad(centroid.Latitude);
            double cLon = Algorithms.Deg2Rad(centroid.Longitude);
            double d_lat_m, d_lon_m;

            for (int i = 0; i < bases.Length; i++)
            {
                Algorithms.GetDeltasByGeopoints(cLat, cLon,
                    Algorithms.Deg2Rad(bases[i].Latitude), Algorithms.Deg2Rad(bases[i].Longitude),
                    el,
                    out d_lat_m, out d_lon_m);

                result[i].X = d_lon_m;
                result[i].Y = d_lat_m;
                result[i].Z = bases[i].Depth;
                result[i].D = bases[i].SlantRange;
            }

            return result;
        }

        /// <summary>
        /// Converts array of GeoPoint3DT to an array of TOABasePoints, converts lat/lon to metric local CS whose center is in specified centroid
        /// </summary>
        /// <param name="bases">An array of base points with measured time of arrival</param>
        /// <param name="centroid">Start of new local coordinate system</param>
        /// <param name="el">ellipsoid</param>
        /// <returns>An array of TOABasePoints</returns>
        public static TOABasePoint[] ConvertToLCS(GeoPoint3DT[] bases, GeoPoint centroid, Ellipsoid el)
        {
            TOABasePoint[] result = new TOABasePoint[bases.Length];

            double cLat = Algorithms.Deg2Rad(centroid.Latitude);
            double cLon = Algorithms.Deg2Rad(centroid.Longitude);
            double d_lat_m, d_lon_m;

            for (int i = 0; i < bases.Length; i++)
            {
                Algorithms.GetDeltasByGeopoints(cLat, cLon,
                    Algorithms.Deg2Rad(bases[i].Latitude), Algorithms.Deg2Rad(bases[i].Longitude),
                    el,
                    out d_lat_m, out d_lon_m);

                result[i].X = d_lon_m;
                result[i].Y = d_lat_m;
                result[i].Z = bases[i].Depth;
                result[i].D = bases[i].TOASec;
            }

            return result;
        }

        /// <summary>
        /// Build baselines on collection of base points (combinatoric)
        /// </summary>
        /// <param name="bases">Base points</param>
        /// <param name="velocity">Signal propagation velocity in m/s</param>
        /// <returns></returns>
        public static TDOABaseline[] BuildBaseLines(TOABasePoint[] bases, double velocity)
        {
            List<TDOABaseline> result = new List<TDOABaseline>();

            for (int i = 0; i < bases.Length - 1; i++)
                for (int j = i + 1; j < bases.Length; j++)
                {
                    result.Add(new TDOABaseline(bases[i].X, bases[i].Y, bases[i].Z,
                                                bases[j].X, bases[j].Y, bases[j].Z,
                                                (bases[i].D - bases[j].D) * velocity));
                }

            return result.ToArray();
        }
        
        /// <summary>
        /// Solves TOA navigation problem: searches for a target location by base points - points with known locations
        /// and measured slant ranges to the target. Nelder-Mead (simplex) method is used with preliminary 1D optimization
        /// </summary>
        /// <param name="bases">A collection of base points</param>
        /// <param name="lat_prev_deg">Previous location latitude, signed degrees</param>
        /// <param name="lon_prev_deg">Previous location longitude, signed degrees</param>
        /// <param name="z_m">Target's depth, known by direct measurement</param>
        /// <param name="maxIterations">Nelder-Mead iterations limit</param>
        /// <param name="precisionThreshold">Precision threshold</param>
        /// <param name="simplexSize">Initial size of simplex, meters</param>
        /// <param name="el">Reference ellipsoid</param>
        /// <param name="lat_deg">Found latitude of the target</param>
        /// <param name="lon_deg">Found longitude of the target</param>
        /// <param name="radialError">Radial error. Square root from the final value of residual function</param>
        /// <param name="itCnt">Number of iterations taken</param>
        public static void TOA_Locate2D(GeoPoint3DD[] bases,
                                        double lat_prev_deg, double lon_prev_deg, double z_m,
                                        int maxIterations, double precisionThreshold, double simplexSize,
                                        Ellipsoid el,
                                        out double lat_deg, out double lon_deg, out double radialError, out int itCnt)
        {
            double xPrev, yPrev;
            double xBest = 0;
            double yBest = 0;
            var basesCentroid = GetPointsCentroid2D(bases);
            var basePoints = ConvertToLCS(bases, basesCentroid, el);

            if (double.IsNaN(lat_prev_deg) || double.IsNaN(lon_prev_deg))
            {
                Algorithms.TOA_Circles1D_Solve(basePoints, z_m, Algorithms.PI_DBY_180, 10, 0.1, out xPrev, out yPrev, out radialError);                
            }
            else
            {
                Algorithms.GetDeltasByGeopoints(Algorithms.Deg2Rad(basesCentroid.Latitude), Algorithms.Deg2Rad(basesCentroid.Longitude),
                                                Algorithms.Deg2Rad(lat_prev_deg), Algorithms.Deg2Rad(lon_prev_deg), 
                                                el, 
                                                out yPrev, out xPrev);
            }

            Algorithms.TOA_NLM2D_Solve(basePoints, xPrev, yPrev, z_m,
                                       maxIterations, precisionThreshold, simplexSize,
                                       out xBest, out yBest, out radialError, out itCnt);

            Algorithms.GeopointOffsetByDeltas(Algorithms.Deg2Rad(basesCentroid.Latitude), Algorithms.Deg2Rad(basesCentroid.Longitude),
                                              yBest, xBest,
                                              el, 
                                              out yPrev, out xPrev);

            lat_deg = Algorithms.Rad2Deg(yPrev);
            lon_deg = Algorithms.Rad2Deg(xPrev);
        }

        /// <summary>
        /// Solves TOA navigation problem: searches for a target location by base points - points with known locations
        /// and measured slant ranges to the target. Nelder-Mead (simplex) method is used with preliminary 1D optimization
        /// </summary>
        /// <param name="bases">A collection of base points</param>
        /// <param name="lat_prev_deg">Previous location latitude, signed degrees</param>
        /// <param name="lon_prev_deg">Previous location longitude, signed degrees</param>
        /// <param name="z_m">Target's depth, known by direct measurement</param>
        /// <param name="maxIterations">Nelder-Mead iterations limit</param>
        /// <param name="precisionThreshold">Precision threshold</param>
        /// <param name="simplexSize">Initial size of simplex, meters</param>
        /// <param name="el">Reference ellipsoid</param>
        /// <param name="lat_deg">Found latitude of the target</param>
        /// <param name="lon_deg">Found longitude of the target</param>
        /// <param name="radialError">Radial error. Square root from the final value of residual function</param>
        /// <param name="itCnt">Number of iterations taken</param>
        public static void TOA_Locate3D(GeoPoint3DD[] bases,
                                        double lat_prev_deg, double lon_prev_deg, double prev_z_m,
                                        int maxIterations, double precisionThreshold, double simplexSize,
                                        Ellipsoid el,
                                        out double lat_deg, out double lon_deg, out double z_m, out double radialError, out int itCnt)
        {
            double xPrev, yPrev;
            double xBest = 0;
            double yBest = 0;
            double zBest = 0;
            var basesCentroid = GetPointsCentroid2D(bases);
            var basePoints = ConvertToLCS(bases, basesCentroid, el);

            if (double.IsNaN(prev_z_m))
                prev_z_m = bases[0].Depth;

            if (double.IsNaN(lat_prev_deg) || double.IsNaN(lon_prev_deg))
            {
                lat_prev_deg = basesCentroid.Latitude;
                lon_prev_deg = basesCentroid.Longitude;
            }
                        
            Algorithms.GetDeltasByGeopoints(Algorithms.Deg2Rad(basesCentroid.Latitude), Algorithms.Deg2Rad(basesCentroid.Longitude),
                                            Algorithms.Deg2Rad(lat_prev_deg), Algorithms.Deg2Rad(lon_prev_deg),
                                            el,
                                            out yPrev, out xPrev);            

            Algorithms.TOA_NLM3D_Solve(basePoints, xPrev, yPrev, prev_z_m,
                                       maxIterations, precisionThreshold, simplexSize,
                                       out xBest, out yBest, out zBest, out radialError, out itCnt);

            Algorithms.GeopointOffsetByDeltas(Algorithms.Deg2Rad(basesCentroid.Latitude), Algorithms.Deg2Rad(basesCentroid.Longitude),
                                              yBest, xBest,
                                              el,
                                              out yPrev, out xPrev);

            lat_deg = Algorithms.Rad2Deg(yPrev);
            lon_deg = Algorithms.Rad2Deg(xPrev);
            z_m = zBest;
        }

        /// <summary>
        /// Solves TDOA navigation problem: searches for a target location by base points - points with known locations
        /// and measured times of arrival. Nelder-Mead (simplex) method is used.
        /// </summary>
        /// <param name="bases">A collection of base points</param>
        /// <param name="lat_prev_deg">Previous location latitude, signed degrees. Set it to NaN if previous location is unknown</param>
        /// <param name="lon_prev_deg">Previous location longitude, signed degrees. Set it to NaN if previous location is unknown</param>
        /// <param name="z_m">Target's depth, known by direct measurement</param>
        /// <param name="maxIterations">Nelder-Mead iterations limit</param>
        /// <param name="precisionThreshold">Precision threshold</param>
        /// <param name="simplexSize">Initial size of simplex, meters</param>
        /// <param name="el">Reference ellipsoid</param>
        /// <param name="velocity">Signal propagation velocity, e.g. speed of sound in m/s</param>
        /// <param name="lat_deg">Found latitude of the target</param>
        /// <param name="lon_deg">Found longitude of the target</param>
        /// <param name="radialError">Radial error. Square root from the final value of residual function</param>
        /// <param name="itCnt">Number of iterations taken</param>
        public static void TDOA_Locate2D(GeoPoint3DT[] bases,
                                         double lat_prev_deg, double lon_prev_deg, double z_m,
                                         int maxIterations, double precisionThreshold, double simplexSize,
                                         Ellipsoid el,
                                         double velocity,
                                         out double lat_deg, out double lon_deg, out double radialError, out int itCnt)
        {
            double xPrev = 0;
            double yPrev = 0;
            double xBest = 0;
            double yBest = 0;

            var basesCentroid = GetPointsCentroid2D(bases);
            var basePoints = ConvertToLCS(bases, basesCentroid, el);
            var baseLines = BuildBaseLines(basePoints, velocity);

            if (!double.IsNaN(lat_prev_deg) && !double.IsNaN(lon_prev_deg))
            {
                Algorithms.GetDeltasByGeopoints(Algorithms.Deg2Rad(basesCentroid.Latitude), Algorithms.Deg2Rad(basesCentroid.Longitude),
                    Algorithms.Deg2Rad(lat_prev_deg), Algorithms.Deg2Rad(lon_prev_deg), el, out yPrev, out xPrev);
            }

            Algorithms.TDOA_NLM2D_Solve(baseLines, xPrev, yPrev, z_m,
                                        maxIterations, precisionThreshold, simplexSize,
                                        out xBest, out yBest, out radialError, out itCnt);
            
            Algorithms.GeopointOffsetByDeltas(Algorithms.Deg2Rad(basesCentroid.Latitude), Algorithms.Deg2Rad(basesCentroid.Longitude), 
                                              yBest, xBest,                                             
                                              el,
                                              out yPrev, out xPrev);

            lat_deg = Algorithms.Rad2Deg(yPrev);
            lon_deg = Algorithms.Rad2Deg(xPrev);
        }

        /// <summary>
        /// Solves TDOA navigation problem: searches for a target location by base points - points with known locations
        /// and measured times of arrival. Nelder-Mead (simplex) method is used.
        /// </summary>
        /// <param name="bases">A collection of base points</param>
        /// <param name="lat_prev_deg">Previous location latitude, signed degrees. Set it to NaN if previous location is unknown</param>
        /// <param name="lon_prev_deg">Previous location longitude, signed degrees. Set it to NaN if previous location is unknown</param>
        /// <param name="z_m">Previous estimation of Target's depth. Set it to NaN of previos depth is unknown</param>
        /// <param name="maxIterations">Nelder-Mead iterations limit</param>
        /// <param name="precisionThreshold">Precision threshold</param>
        /// <param name="simplexSize">Initial size of simplex, meters</param>
        /// <param name="el">Reference ellipsoid</param>
        /// <param name="velocity">Signal propagation velocity, e.g. speed of sound in m/s</param>
        /// <param name="lat_deg">Found latitude of the target</param>
        /// <param name="lon_deg">Found longitude of the target</param>
        /// <param name="radialError">Radial error. Square root from the final value of residual function</param>
        /// <param name="itCnt">Number of iterations taken</param>
        public static void TDOA_Locate3D(GeoPoint3DT[] bases,
                                         double lat_prev_deg, double lon_prev_deg, double prev_z_m,
                                         int maxIterations, double precisionThreshold, double simplexSize,
                                         Ellipsoid el,
                                         double velocity,
                                         out double lat_deg, out double lon_deg, out double z_m, out double radialError, out int itCnt)
        {
            double xPrev = 0;
            double yPrev = 0;
            double xBest = 0;
            double yBest = 0;

            var basesCentroid = GetPointsCentroid2D(bases);
            var basePoints = ConvertToLCS(bases, basesCentroid, el);
            var baseLines = BuildBaseLines(basePoints, velocity);

            if (!double.IsNaN(lat_prev_deg) && !double.IsNaN(lon_prev_deg))
            {
                Algorithms.GetDeltasByGeopoints(Algorithms.Deg2Rad(basesCentroid.Latitude), Algorithms.Deg2Rad(basesCentroid.Longitude),
                    Algorithms.Deg2Rad(lat_prev_deg), Algorithms.Deg2Rad(lon_prev_deg), el, out yPrev, out xPrev);
            }
            
            if (double.IsNaN(prev_z_m))
                prev_z_m = bases[0].Depth;

            Algorithms.TDOA_NLM3D_Solve(baseLines, xPrev, yPrev, prev_z_m,
                                        maxIterations, precisionThreshold, simplexSize,
                                        out xBest, out yBest, out z_m, out radialError, out itCnt);
            

            
            Algorithms.GeopointOffsetByDeltas(Algorithms.Deg2Rad(basesCentroid.Latitude), Algorithms.Deg2Rad(basesCentroid.Longitude),
                                              yBest, xBest,
                                              el,
                                              out yPrev, out xPrev);

            lat_deg = Algorithms.Rad2Deg(yPrev);
            lon_deg = Algorithms.Rad2Deg(xPrev);
        }

        #endregion
    }
}
