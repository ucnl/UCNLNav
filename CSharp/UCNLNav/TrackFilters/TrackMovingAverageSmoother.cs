using System;
using System.Collections.Generic;

namespace UCNLNav.TrackFilters
{

    public class TrackMovingAverageSmootherCartesian : ITrackCartesianFilter
    {
        #region Properties

        public int FilterSize { get; private set; }

        public bool IsEmpty
        {
            get { return points.Count == 0; }
        }

        List<MPoint> points;        

        double max_dist_m_to_reset = 3000;

        double prev_x = double.NaN;
        double prev_y = double.NaN;

        #endregion

        #region Constructor

        public TrackMovingAverageSmootherCartesian(int filterSize)
            : this(filterSize, 3000)
        {
        }

        public TrackMovingAverageSmootherCartesian(int filterSize, double reset_delta_m)
        {
            if (reset_delta_m <= 0)
                throw new ArgumentOutOfRangeException("reset_delta_m should be greater than zero");

            if (filterSize <= 0)
                throw new ArgumentOutOfRangeException("filterSize should be greater than zero");

            max_dist_m_to_reset = reset_delta_m;
            FilterSize = filterSize;
            points = new List<MPoint>(filterSize);
        }

        #endregion

        #region Methods

        private double Dist2D(double x1, double y1, double x2, double y2)
        {
            double dx = x2 - x1;
            double dy = y2 - y1;

            return Math.Sqrt(dx*dx + dy*dy);
        }

        public void Reset()
        {            
            points.Clear();
        }

        public bool Process(double x, double y, double z, DateTime inTS,
            out double xflt, out double yflt, out double zflt, out DateTime outTS)
        {
            zflt = z;
            outTS = inTS;

            if (points.Count == 0)
            {
                points.Add(new MPoint(x, y));
                xflt = x;
                yflt = y;
                prev_x = x;
                prev_y = y;
            }
            else
            {
                if (Dist2D(prev_x, prev_y, x, y) >= max_dist_m_to_reset)
                {
                    Reset();
                    xflt = x;
                    yflt = y;
                }
                else
                {
                    if (points.Count >= FilterSize)
                        points.RemoveAt(0);

                    points.Add(new MPoint(x, y));

                    double meanx = 0.0, meany = 0.0;
                    for (int i = 0; i < points.Count; i++)
                    {
                        meanx += points[i].X * (i + 1);
                        meany += points[i].Y * (i + 1);
                    }

                    double fWeight = (points.Count + points.Count * points.Count) / 2.0;
                    meanx /= fWeight;
                    meany /= fWeight;

                    xflt = meanx; yflt = meany;
                    prev_x = xflt; prev_y = yflt;
                }
            }

            return true;
        }

        #endregion
    }

    public class TrackMovingAverageSmoother : ITrackFilter
    {
        #region Properties

        public int FilterSize { get; private set; }

        public bool IsEmpty
        {
            get { return points.Count == 0; }
        }

        List<MPoint> points;

        double anchor_lat_rad = double.NaN;
        double anchor_lon_rad = double.NaN;

        double max_dist_m_to_reset = 3000;

        double prev_lat_rad = double.NaN;
        double prev_lon_rad = double.NaN;

        #endregion

        #region Constructor

        public TrackMovingAverageSmoother(int filterSize)
            : this(filterSize, 3000)
        {
        }

        public TrackMovingAverageSmoother(int filterSize, double reset_delta_m)
        {
            if (reset_delta_m <= 0)
                throw new ArgumentOutOfRangeException("reset_delta_m should be greater than zero");

            if (filterSize <= 0)
                throw new ArgumentOutOfRangeException("filterSize should be greater than zero");

            max_dist_m_to_reset = reset_delta_m;
            FilterSize = filterSize;
            points = new List<MPoint>(filterSize);
        }

        #endregion

        #region Methods

        public void Reset()
        {
            anchor_lat_rad = double.NaN;
            anchor_lon_rad = double.NaN;
            points.Clear();
        }

        public bool Process(double inLat_rad, double inLon_rad, double inDpt_m, DateTime inTS,
            out double outLat_rad, out double outLon_rad, out double outDpt_m, out DateTime outTS)
        {
            outDpt_m = inDpt_m;
            outTS = inTS;

            if (points.Count == 0)
            {
                anchor_lat_rad = inLat_rad;
                anchor_lon_rad = inLon_rad;
                points.Add(new MPoint(0, 0));
                prev_lat_rad = anchor_lat_rad;
                prev_lon_rad = anchor_lon_rad;

                outLat_rad = inLat_rad;
                outLon_rad = inLon_rad;
            }
            else
            {
                if (Algorithms.HaversineInverse(prev_lat_rad, prev_lon_rad, inLat_rad, inLon_rad, Algorithms.WGS84Ellipsoid.MajorSemiAxis_m) >= max_dist_m_to_reset)
                {
                    Reset();

                    outLat_rad = inLat_rad;
                    outLon_rad = inLon_rad;
                }
                else
                {

                    if (points.Count >= FilterSize)
                        points.RemoveAt(0);

                    double deltay = 0, deltax = 0;
                    Algorithms.GetDeltasByGeopoints_WGS84(anchor_lat_rad, anchor_lon_rad,
                        inLat_rad, inLon_rad,
                        out deltay, out deltax);

                    points.Add(new MPoint(deltax, deltay));

                    double meanx = 0.0, meany = 0.0;
                    for (int i = 0; i < points.Count; i++)
                    {
                        meanx += points[i].X * (i + 1);
                        meany += points[i].Y * (i + 1);
                    }

                    double fWeight = (points.Count + points.Count * points.Count) / 2.0;
                    meanx /= fWeight;
                    meany /= fWeight;

                    Algorithms.GeopointOffsetByDeltas_WGS84(anchor_lat_rad, anchor_lon_rad, meany, meanx, out outLat_rad, out outLon_rad);

                    prev_lat_rad = outLat_rad;
                    prev_lon_rad = outLon_rad;
                }
            }

            return true;
        }

        #endregion
    }
}
