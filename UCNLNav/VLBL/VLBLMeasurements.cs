using System;
using System.Collections.Generic;

namespace UCNLNav.VLBL
{
    public class VLBLMeasurements<T> where T : GeoPoint3D, IVLBLMeasurement
    {
        #region Properties

        List<T> byAge;
        List<T> byBearingFromRefPoint;
        List<T> byDistanceToRefPoint;

        Comparison<T> byBearingComparison;
        Comparison<T> byDistanceComparison;

        GeoPoint refPoint;
        bool isRefPointSet;

        double aRange;
        double aRangeStart;
        double aRangeEnd;
        double angularGap;
        double angularRange;

        public int Capacity { get; private set; }
        public int BaseSize { get; private set; }

        public bool IsBaseExists
        {
            get
            {
                return byDistanceToRefPoint.Count >= BaseSize + 1;
            }
        }

        public GeoPoint RefPoint
        {
            get { return new GeoPoint(refPoint.Latitude, refPoint.Longitude); }
        }

        #endregion

        #region Constructor

        public VLBLMeasurements(int capacity, int baseSize)
        {
            if (capacity < 3)
                throw new ArgumentOutOfRangeException("Maximal size of measurements list cannot be less than 3");

            if (baseSize < 3)
                throw new ArgumentOutOfRangeException("Base size cannot be less than 3");

            Capacity = capacity;
            BaseSize = baseSize;

            refPoint = new GeoPoint(double.NaN, double.NaN);

            byAge = new List<T>();
            byBearingFromRefPoint = new List<T>();
            byDistanceToRefPoint = new List<T>();

            byBearingComparison = new Comparison<T>((x, y) => x.BearingFromRefPoint.CompareTo(y.BearingFromRefPoint));
            byDistanceComparison = new Comparison<T>((x, y) => x.DistanceToRefPoint.CompareTo(y.DistanceToRefPoint));
        }

        #endregion

        #region Methods

        #region Private

        private void Resort()
        {
            byBearingFromRefPoint.Sort(byBearingComparison);
            byDistanceToRefPoint.Sort(byDistanceComparison);
        }

        private void UpdateRefPoint()
        {
            foreach (var item in byAge)
            {
                item.UpdateRefPoint(refPoint);
            }
        }

        private void UpdateAngularRange()
        {
            if (byAge.Count < 2)
                return;

            int maxGapStartIdx = 0;
            int maxGapEndIdx = 1;
            double maxGap = byBearingFromRefPoint[maxGapEndIdx].BearingFromRefPoint - byBearingFromRefPoint[maxGapStartIdx].BearingFromRefPoint;
            double gap;

            int idx = 2;
            int lIdx, rIdx;
            while (idx <= byBearingFromRefPoint.Count)
            {
                lIdx = idx - 1;
                rIdx = idx;

                if (rIdx >= byBearingFromRefPoint.Count)
                    rIdx = rIdx % byBearingFromRefPoint.Count;

                gap = byBearingFromRefPoint[rIdx].BearingFromRefPoint - byBearingFromRefPoint[lIdx].BearingFromRefPoint;

                if (gap < 0)
                    gap = gap + 360;

                if (gap > maxGap)
                {
                    maxGap = gap;
                    maxGapStartIdx = lIdx - 1;
                    maxGapEndIdx = rIdx;
                }

                idx = idx + 1;
            }

            aRangeStart = byBearingFromRefPoint[maxGapEndIdx].BearingFromRefPoint;
            aRangeEnd = byBearingFromRefPoint[maxGapStartIdx].BearingFromRefPoint;
            aRange = aRangeEnd - aRangeStart;

            if (aRange < 0)
                aRange += 360;

            angularRange = 360 - maxGap;
            angularGap = maxGap;
        }

        #endregion

        #region Public

        public void Clear()
        {
            byAge.Clear();
            byBearingFromRefPoint.Clear();
            byDistanceToRefPoint.Clear();
            refPoint = new GeoPoint(double.NaN, double.NaN);
            isRefPointSet = false;

            aRange = 0;
            aRangeStart = 0;
            aRangeEnd = 0;
            angularGap = 0;
            angularRange = 0;
        }

        public void SetRefPoint(GeoPoint newRefPoint)
        {
            if (double.IsNaN(newRefPoint.Latitude) || (double.IsNaN(newRefPoint.Longitude)))
                throw new ArgumentOutOfRangeException();

            refPoint.Latitude = newRefPoint.Latitude;
            refPoint.Longitude = newRefPoint.Longitude;
            UpdateRefPoint();
            Resort();
            UpdateAngularRange();
            isRefPointSet = true;
        }

        public void Add(T item)
        {
            bool refPointUpdated = false;

            if (!isRefPointSet)
            {
                if (byAge.Count == 0)
                {
                    refPoint = new GeoPoint(item.Latitude, item.Longitude);
                }
                else
                {
                    refPoint = Navigation.GetPointsCentroid2D(byAge);
                }

                refPointUpdated = true;
            }

            if (byAge.Count >= Capacity)
            {
                var oldest = byAge[byAge.Count - 1];
                byAge.Remove(oldest);
                byBearingFromRefPoint.Remove(oldest);
                byDistanceToRefPoint.Remove(oldest);
            }

            byAge.Add(item);
            byDistanceToRefPoint.Add(item);
            byBearingFromRefPoint.Add(item);

            if (refPointUpdated)
            {
                UpdateRefPoint();
            }
            else
            {
                item.UpdateRefPoint(refPoint);
            }

            Resort();
            UpdateAngularRange();
        }

        public List<T> GetBase()
        {
            if (!IsBaseExists)
                throw new InvalidOperationException("Base does not exists");

            List<T> result = new List<T>();

            #region arrange desired angles

            double delta = aRange / (BaseSize - 1);

            if (angularGap < delta)
            {
                double deltaDec = (delta - angularGap) * (BaseSize - 1) / (BaseSize - 2);
                delta = delta - deltaDec / BaseSize;
            }

            double[] desiredAngles = new double[BaseSize];

            for (int i = 0; i < BaseSize; i++)
            {
                desiredAngles[i] = aRangeStart + i * delta;
                while (desiredAngles[i] > 360)
                    desiredAngles[i] -= 360;
            }

            #endregion

            #region select nearest to desired angles

            double minDiff = double.MaxValue;
            int minDiffIdx = 0;
            double diff;
            for (int i = 0; i < BaseSize; i++)
            {
                minDiff = double.MaxValue;
                minDiffIdx = 0;

                for (int j = 0; j < byBearingFromRefPoint.Count; j++)
                {
                    if (!result.Contains(byBearingFromRefPoint[j]))
                    {
                        diff = desiredAngles[i] - byBearingFromRefPoint[j].BearingFromRefPoint;
                        if (diff < 0)
                            diff += 360;

                        if (diff < minDiff)
                        {
                            minDiff = diff;
                            minDiffIdx = j;
                        }
                    }
                }

                result.Add(byBearingFromRefPoint[minDiffIdx]);
            }

            #endregion

            if (!result.Contains(byDistanceToRefPoint[0]))
                result.Add(byDistanceToRefPoint[0]);

            result.Sort(byDistanceComparison);

            return result;
        }

        #endregion

        #endregion
    }
}
