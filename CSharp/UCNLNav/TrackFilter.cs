using System;
using System.Collections.Generic;

namespace UCNLNav
{
    public class TrackFilter
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

        public TrackFilter(int filterSize)
            : this(filterSize, 3000)
        {
        }

        public TrackFilter(int filterSize, double reset_delta_m)
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

        public GeoPoint Filter(double lat_deg, double lon_deg)
        {
            if (points.Count == 0)
            {
                anchor_lat_rad = Algorithms.Deg2Rad(lat_deg);
                anchor_lon_rad = Algorithms.Deg2Rad(lon_deg);
                points.Add(new MPoint(0, 0));
                prev_lat_rad = anchor_lat_rad;
                prev_lon_rad = anchor_lon_rad;
                return new GeoPoint(lat_deg, lon_deg);
            }
            else
            {                
                double lat_rad = Algorithms.Deg2Rad(lat_deg);
                double lon_rad = Algorithms.Deg2Rad(lon_deg);

                if (Algorithms.HaversineInverse(prev_lat_rad, prev_lon_rad, lat_rad, lon_rad, Algorithms.WGS84Ellipsoid.MajorSemiAxis_m) >= max_dist_m_to_reset)
                {
                    Reset();
                    return new GeoPoint(lat_deg, lon_deg);
                }

                if (points.Count >= FilterSize)
                    points.RemoveAt(0);

                double deltay = 0, deltax = 0;
                Algorithms.GetDeltasByGeopoints_WGS84(anchor_lat_rad, anchor_lon_rad,
                    lat_rad, lon_rad,
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

                
                Algorithms.GeopointOffsetByDeltas_WGS84(anchor_lat_rad, anchor_lon_rad, meany, meanx, out lat_rad, out lon_rad);

                prev_lat_rad = lat_rad;
                prev_lon_rad = lon_rad;

                return new GeoPoint(Algorithms.Rad2Deg(lat_rad), Algorithms.Rad2Deg(lon_rad));
            }
        }

        #endregion
    }   
}
