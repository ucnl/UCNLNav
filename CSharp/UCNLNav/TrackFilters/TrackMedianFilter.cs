using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace UCNLNav.TrackFilters
{
    public class TrackMedianFilter : ITrackFilter
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

        public TrackMedianFilter(int filterSize)
            : this(filterSize, 3000)
        {
        }

        public TrackMedianFilter(int filterSize, double reset_delta_m)
        {
            if (reset_delta_m <= 0)
                throw new ArgumentOutOfRangeException("reset_delta_m should be greater than zero");

            if ((filterSize <= 0) || (filterSize % 2 == 0))
                throw new ArgumentOutOfRangeException("filterSize should be an odd value greater than zero");

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

                    double[] xs = new double[points.Count];
                    double[] ys = new double[points.Count];

                    for (int i = 0; i < points.Count; i++)
                    {
                        xs[i] = points[i].X;
                        ys[i] = points[i].Y;
                    }

                    Array.Sort(xs);
                    Array.Sort(ys);

                    double vx = xs[xs.Length / 2];
                    double vy = ys[ys.Length / 2];

                    Algorithms.GeopointOffsetByDeltas_WGS84(anchor_lat_rad, anchor_lon_rad, vy, vx, out outLat_rad, out outLon_rad);

                    prev_lat_rad = outLat_rad;
                    prev_lon_rad = outLon_rad;
                }
            }

            return true;
        }

        #endregion
    }
}
