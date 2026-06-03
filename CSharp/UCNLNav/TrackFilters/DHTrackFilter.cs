using System;
using System.Collections.Generic;

namespace UCNLNav.TrackFilters
{

    public class GeoPoint3DTd : GeoPoint3D
    {
        #region Properties

        public DateTime TimeStamp;
        public double Dist2Prev_m;
        public double Time2Prev_s;

        #endregion

        #region Constructor

        public GeoPoint3DTd(double lat, double lon, double dpt, DateTime ts, double dst2prev_m, double time2prev_s)
            : base(lat, lon, dpt)
        {
            TimeStamp = ts;
            Dist2Prev_m = dst2prev_m;
            Time2Prev_s = time2prev_s;
        }

        public GeoPoint3DTd(double lat, double lon, double dpt, DateTime ts, GeoPoint3DTd prevPoint)
            : base(lat, lon, dpt)
        {
            TimeStamp = ts;

            Dist2Prev_m = Algorithms.HaversineInverse(lat, lon,
                prevPoint.Latitude, prevPoint.Longitude,
                Algorithms.WGS84Ellipsoid.MajorSemiAxis_m);

            Time2Prev_s = ts.Subtract(prevPoint.TimeStamp).TotalSeconds;
        }

        public GeoPoint3DTd(double lat, double lon, double dpt, DateTime ts)
            : base(lat, lon, dpt)
        {
            TimeStamp = ts;
            Dist2Prev_m = double.NaN;
            Time2Prev_s = double.NaN;
        }

        #endregion

        #region Methods

        public override string ToString()
        {
            return string.Format("{0}, D2P: {1:F03} m, T2P: {2:F03} s", base.ToString(), Dist2Prev_m, Time2Prev_s);
        }

        #endregion
    }

    public class MPoint3DTd : MPoint3D
    {
        #region Properties

        public DateTime TimeStamp;
        public double Dist2Prev_m;
        public double Time2Prev_s;

        #endregion

        #region Constructor

        public MPoint3DTd(double x, double y, double z, DateTime ts, double dst2prev_m, double time2prev_s)
            : base(x, y, z)
        {
            TimeStamp = ts;
            Dist2Prev_m = dst2prev_m;
            Time2Prev_s = time2prev_s;
        }

        public MPoint3DTd(double x, double y, double z, DateTime ts, MPoint3DTd prevPoint)
            : base(x, y, z)
        {
            TimeStamp = ts;
            Dist2Prev_m = Math.Sqrt((x - prevPoint.X) * (x - prevPoint.X) + (y - prevPoint.Y) * (y - prevPoint.Y));
            Time2Prev_s = ts.Subtract(prevPoint.TimeStamp).TotalSeconds;
        }

        public MPoint3DTd(double x, double y, double z, DateTime ts)
            : base(x, y, z)
        {
            TimeStamp = ts;
            Dist2Prev_m = double.NaN;
            Time2Prev_s = double.NaN;
        }

        #endregion

        #region Methods

        public override string ToString()
        {
            return string.Format("{0}, D2P: {1:F03} m, T2P: {2:F03} s", base.ToString(), Dist2Prev_m, Time2Prev_s);
        }

        #endregion
    }

    public class DHTrackFilterCartesian : ITrackCartesianFilter
    {
        #region Properties

        List<List<MPoint3DTd>> sides;

        int pSideIdx = 0;
        int sSideIdx { get { return pSideIdx == 0 ? 1 : 0; } }

        double maxSpeedMps = 1;
        public double MaxSpeedMps
        {
            get => maxSpeedMps;
            set
            {
                if (value >= 0.0)
                    maxSpeedMps = value;
                else
                    throw new ArgumentOutOfRangeException("value", "Value must be non-negative");
            }
        }

        double dstThreshold_m = 5;
        public double DstThreshold_m
        {
            get => dstThreshold_m;
            set
            {
                if (value >= 0.0)
                    dstThreshold_m = value;
                else
                    throw new ArgumentOutOfRangeException("value", "Value must be non-negative");
            }
        }

        int fifoSize = 0;
        public int FIFOSize { get { return fifoSize; } }

        #endregion

        #region Constructor

        public DHTrackFilterCartesian()
            : this(8, 1, 5)
        {

        }

        public DHTrackFilterCartesian(int fifoSize, double maxSpeedMps, double dstTrhesholdm)
        {
            if (fifoSize < 2)
                throw new ArgumentOutOfRangeException("fifoSize", "should be greater than 1");

            if (maxSpeedMps <= 0)
                throw new ArgumentOutOfRangeException("maxSpeedMps", "should be greater than 0 m/s");

            if (dstThreshold_m <= 0)
                throw new ArgumentOutOfRangeException("dstThresholdm", "should be greater than 0 m");

            this.fifoSize = fifoSize;
            MaxSpeedMps = maxSpeedMps;
            DstThreshold_m = dstTrhesholdm;
            sides = new List<List<MPoint3DTd>>
            {
                new List<MPoint3DTd>(),
                new List<MPoint3DTd>()
            };
        }

        #endregion

        #region Methods

        #region Private

        private void AddPoint(int listIdx, MPoint3DTd newPoint)
        {
            if (sides[listIdx].Count + 1 > fifoSize)
                sides[listIdx].RemoveAt(0);

            sides[listIdx].Add(newPoint);
        }

        private void UpdateStatistics()
        {
            if ((sides[0].Count == fifoSize) && (sides[1].Count == fifoSize))
            {
                double[] meanDst = new double[2];

                for (int i = 0; i < sides.Count; i++)
                {
                    meanDst[i] = 0.0;
                    for (int j = 0; j < sides[i].Count; j++)
                        meanDst[i] += sides[i][j].Dist2Prev_m;

                    meanDst[i] /= sides[i].Count;
                }

                double[] sigmas = new double[2];
                for (int i = 0; i < sides.Count; i++)
                {
                    for (int j = 0; j < sides[i].Count; j++)
                        sigmas[i] += Math.Pow(sides[i][j].Dist2Prev_m - meanDst[i], 2);

                    sigmas[i] = Math.Sqrt(sigmas[i]);
                }

                if (sigmas[pSideIdx] > sigmas[sSideIdx])
                    pSideIdx = sSideIdx;
            }
        }

        #endregion

        #region Public

        public void Reset()
        {
            sides[0].Clear();
            sides[1].Clear();
            pSideIdx = 0;
        }

        public bool Process(double x, double y, double z, DateTime inTS,
            out double xflt, out double yflt, out double zflt, out DateTime outTS)
        {
            bool result = false;
            xflt = double.NaN;
            yflt = double.NaN;
            zflt = z;
            outTS = DateTime.MinValue;

            if (sides[pSideIdx].Count == 0)
            {
                AddPoint(pSideIdx, new MPoint3DTd(x, y, z, inTS));
                xflt = x;
                yflt = y;
                zflt = z;
                outTS = inTS;
                result = true;
            }
            else
            {
                UpdateStatistics();

                MPoint3DTd pLastPoint = sides[pSideIdx][sides[pSideIdx].Count - 1];
                MPoint3DTd pTestPoint = new MPoint3DTd(x, y, z, inTS, pLastPoint);

                if ((pTestPoint.Dist2Prev_m < pTestPoint.Time2Prev_s * maxSpeedMps) ||
                    (pTestPoint.Dist2Prev_m < dstThreshold_m))
                {
                    AddPoint(pSideIdx, pTestPoint);
                    xflt = x;
                    yflt = y;
                    zflt = z;
                    outTS = inTS;
                    result = true;
                }
                else
                {
                    MPoint3DTd sNewPoint;
                    if (sides[sSideIdx].Count > 0)
                    {
                        MPoint3DTd sLastPoint = sides[sSideIdx][sides[sSideIdx].Count - 1];
                        sNewPoint = new MPoint3DTd(x, y, z, inTS, sLastPoint);
                    }
                    else
                    {
                        sNewPoint = new MPoint3DTd(x, y, z, inTS);
                    }

                    AddPoint(sSideIdx, sNewPoint);
                }
            }

            return result;
        }

        #endregion

        #endregion
    }

    public class DHTrackFilter : ITrackFilter
    {
        #region Properties

        List<List<GeoPoint3DTd>> sides;

        int pSideIdx = 0;
        int sSideIdx { get { return pSideIdx == 0 ? 1 : 0; } }

        double maxSpeedMps = 1;
        public double MaxSpeedMps
        {
            get { return maxSpeedMps; }
            set
            {
                if (value >= 0.0)
                    maxSpeedMps = value;
                else
                    throw new ArgumentOutOfRangeException("value", "Value must be non-negative");
            }
        }

        double dstThreshold_m = 5;
        public double DstThreshold_m
        {
            get { return dstThreshold_m; }
            set
            {
                if (value >= 0.0)
                    dstThreshold_m = value;
                else
                    throw new ArgumentOutOfRangeException("value", "Value must be non-negative");
            }
        }

        int fifoSize = 0;
        public int FIFOSize { get { return fifoSize; } }

        #endregion

        #region Constructor

        public DHTrackFilter()
            : this(8, 1, 5)
        {

        }

        public DHTrackFilter(int fifoSize, double maxSpeedMps, double dstTrhesholdm)
        {
            if (fifoSize < 2)
                throw new ArgumentOutOfRangeException("fifoSize", "should be greater than 1");

            if (maxSpeedMps <= 0)
                throw new ArgumentOutOfRangeException("maxSpeedMps", "should be greater than 0 m/s");

            if (dstThreshold_m <= 0)
                throw new ArgumentOutOfRangeException("dstThresholdm", "should be greater than 0 m");

            this.fifoSize = fifoSize;
            MaxSpeedMps = maxSpeedMps;
            DstThreshold_m = dstTrhesholdm;
            sides = new List<List<GeoPoint3DTd>>
            {
                new List<GeoPoint3DTd>(),
                new List<GeoPoint3DTd>()
            };
        }

        #endregion

        #region Methods

        #region Private

        private void AddPoint(int listIdx, GeoPoint3DTd newPoint)
        {
            if (sides[listIdx].Count + 1 > fifoSize)
                sides[listIdx].RemoveAt(0);

            sides[listIdx].Add(newPoint);
        }

        private void UpdateStatistics()
        {
            if ((sides[0].Count == fifoSize) && (sides[1].Count == fifoSize))
            {
                double[] meanDst = new double[2];

                for (int i = 0; i < sides.Count; i++)
                {
                    meanDst[i] = 0.0;
                    for (int j = 0; j < sides[i].Count; j++)
                        meanDst[i] += sides[i][j].Dist2Prev_m;

                    meanDst[i] /= sides[i].Count;
                }

                double[] sigmas = new double[2];
                for (int i = 0; i < sides.Count; i++)
                {
                    for (int j = 0; j < sides[i].Count; j++)
                        sigmas[i] += Math.Pow(sides[i][j].Dist2Prev_m - meanDst[i], 2);

                    sigmas[i] = Math.Sqrt(sigmas[i]);
                }

                if (sigmas[pSideIdx] > sigmas[sSideIdx])
                    pSideIdx = sSideIdx;
            }
        }

        #endregion

        #region Public

        public void Reset()
        {
            sides[0].Clear();
            sides[1].Clear();
            pSideIdx = 0;
        }

        public bool Process(double inLat_rad, double inLon_rad, double inDpt_m, DateTime inTS,
            out double outLat_rad, out double outLon_rad, out double outDpt_m, out DateTime outTS)
        {
            bool result = false;
            outLat_rad = double.NaN;
            outLon_rad = double.NaN;
            outDpt_m = double.NaN;
            outTS = DateTime.MinValue;

            if (sides[pSideIdx].Count == 0)
            {
                AddPoint(pSideIdx, new GeoPoint3DTd(inLat_rad, inLon_rad, inDpt_m, inTS));
                outLat_rad = inLat_rad;
                outLon_rad = inLon_rad;
                outDpt_m = inDpt_m;
                outTS = inTS;
                result = true;
            }
            else
            {
                UpdateStatistics();

                GeoPoint3DTd pLastPoint = sides[pSideIdx][sides[pSideIdx].Count - 1];
                GeoPoint3DTd pTestPoint = new GeoPoint3DTd(inLat_rad, inLon_rad, inDpt_m, inTS, pLastPoint);

                if ((pTestPoint.Dist2Prev_m < pTestPoint.Time2Prev_s * maxSpeedMps) ||
                    (pTestPoint.Dist2Prev_m < dstThreshold_m))
                {
                    AddPoint(pSideIdx, pTestPoint);
                    outLat_rad = inLat_rad;
                    outLon_rad = inLon_rad;
                    outDpt_m = inDpt_m;
                    outTS = inTS;
                    result = true;
                }
                else
                {
                    GeoPoint3DTd sNewPoint;
                    if (sides[sSideIdx].Count > 0)
                    {
                        GeoPoint3DTd sLastPoint = sides[sSideIdx][sides[sSideIdx].Count - 1];
                        sNewPoint = new GeoPoint3DTd(inLat_rad, inLon_rad, inDpt_m, inTS, sLastPoint);
                    }
                    else
                    {
                        sNewPoint = new GeoPoint3DTd(inLat_rad, inLon_rad, inDpt_m, inTS);
                    }

                    AddPoint(sSideIdx, sNewPoint);
                }
            }

            return result;
        }

        #endregion

        #endregion
    }
}
