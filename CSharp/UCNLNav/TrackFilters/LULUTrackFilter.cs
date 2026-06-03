using System;
using System.Collections.Generic;

namespace UCNLNav.TrackFilters
{
    public class LULUTrackFilter : ITrackFilter
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

        int wincenter_idx;

        #endregion

        #region Constructor

        public LULUTrackFilter(int filterSize)
        {
            if (filterSize % 2 == 0)
                throw new ArgumentOutOfRangeException("filterSize", "should be a positive odd value");

            FilterSize = filterSize;
            wincenter_idx = filterSize / 2;
            points = new List<MPoint>(filterSize);
        }


        #endregion

        #region Methods

        private void FindMins(out double eval_x, out double eval_y)
        {
            double minx = points[0].X;
            double miny = points[0].Y;
            
            for (int i = 1; i < FilterSize; i++)
            {
               
                if (points[i].X < minx)
                    minx = points[i].X;

                if (points[i].Y < miny)
                    miny = points[i].Y;
            }

            eval_x = minx;
            eval_y = miny;
        }

        private void FindMins(double cval_x, double cval_y, out double eval_x, out double eval_y)
        {
            double minx = points[0].X;
            double miny = points[0].Y;

            double cvx = 0;
            double cvy = 0;

            for (int i = 1; i < FilterSize; i++)
            {
                if (i == wincenter_idx)
                {
                    cvx = cval_x;
                    cvy = cval_y;
                }
                else
                {
                    cvx = points[i].X;
                    cvy = points[i].Y;
                }

                if (cvx < minx)
                    minx = cvx;

                if (cvy < miny)
                    miny = cvy;
            }

            eval_x = minx;
            eval_y = miny;
        }

        private void FindMaxs(out double eval_x, out double eval_y)
        {
            double maxx = points[0].X;
            double maxy = points[0].Y;

            for (int i = 1; i < FilterSize; i++)
            {

                if (points[i].X > maxx)
                    maxx = points[i].X;

                if (points[i].Y > maxy)
                    maxy = points[i].Y;
            }

            eval_x = maxx;
            eval_y = maxy;
        }

        private void FindMaxs(double cval_x, double cval_y, out double eval_x, out double eval_y)
        {
            double maxx = points[0].X;
            double maxy = points[0].Y;

            double cvx = 0;
            double cvy = 0;

            for (int i = 1; i < FilterSize; i++)
            {
                if (i == wincenter_idx)
                {
                    cvx = cval_x;
                    cvy = cval_y;
                }
                else
                {
                    cvx = points[i].X;
                    cvy = points[i].Y;
                }

                if (cvx > maxx)
                    maxx = cvx;

                if (cvy > maxy)
                    maxy = cvy;
            }

            eval_x = maxx;
            eval_y = maxy;
        }



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

            if (points.Count < FilterSize)
            {
                if (points.Count == 0)
                {
                    anchor_lat_rad = inLat_rad;
                    anchor_lon_rad = inLon_rad;
                    points.Add(new MPoint(0, 0));
                } 
                else
                {
                    double deltay = 0, deltax = 0;
                    Algorithms.GetDeltasByGeopoints_WGS84(anchor_lat_rad, anchor_lon_rad,
                        inLat_rad, inLon_rad,
                        out deltay, out deltax);

                    points.Add(new MPoint(deltax, deltay));
                }

                outLat_rad = inLat_rad;
                outLon_rad = inLon_rad;
            }
            else
            {
                points.RemoveAt(0);

                double deltay = 0, deltax = 0;
                Algorithms.GetDeltasByGeopoints_WGS84(anchor_lat_rad, anchor_lon_rad,
                    inLat_rad, inLon_rad,
                    out deltay, out deltax);

                points.Add(new MPoint(deltax, deltay));

                // Тут фильтрация
                FindMins(out double lstage_x, out double lstage_y);
                FindMaxs(lstage_x, lstage_y, out double ustage_x, out double ustage_y);
                //FindMins(ustage_x, ustage_y, out double lustage_x, out double lustage_y);
                //FindMaxs(lustage_x, lustage_y, out double x, out double y);
                double x = ustage_x;
                double y = ustage_y;

                Algorithms.GeopointOffsetByDeltas_WGS84(anchor_lat_rad, anchor_lon_rad, y, x, out outLat_rad, out outLon_rad);

                if (Algorithms.HaversineInverse(anchor_lat_rad, anchor_lon_rad, outLat_rad, outLon_rad, Algorithms.WGS84Ellipsoid.MajorSemiAxis_m) > 1000)
                {
                    Reset();
                }
            }

            return true;
        }


        #endregion
    }
}
