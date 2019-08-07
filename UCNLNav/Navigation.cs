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

            return new GeoPoint(cLat / cnt, cLon / cnt);
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

        #endregion
    }
}
