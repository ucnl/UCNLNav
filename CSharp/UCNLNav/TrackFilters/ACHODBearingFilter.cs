using System;
using System.Collections.Generic;

namespace UCNLNav.TrackFilters
{
    /// <summary>
    /// Точка с пеленгом и дальностью для ACHOD
    /// </summary>
    public class BearingPoint3DTd
    {
        #region Properties

        public double Azimuth_deg;        // абсолютный азимут
        public double Range_m;            // дальность
        public DateTime TimeStamp;
        public double AzimuthDelta_deg;   // угловая разница с предыдущей точкой
        public double RangeDelta_m;       // разница дальности с предыдущей точкой

        #endregion

        #region Constructor

        public BearingPoint3DTd(double azimuth_deg, double range_m, DateTime ts,
            double azimuthDelta_deg, double rangeDelta_m)
        {
            Azimuth_deg = azimuth_deg;
            Range_m = range_m;
            TimeStamp = ts;
            AzimuthDelta_deg = azimuthDelta_deg;
            RangeDelta_m = rangeDelta_m;
        }

        public BearingPoint3DTd(double azimuth_deg, double range_m, DateTime ts,
            BearingPoint3DTd prevPoint)
        {
            Azimuth_deg = azimuth_deg;
            Range_m = range_m;
            TimeStamp = ts;

            AzimuthDelta_deg = AngularDifference(azimuth_deg, prevPoint.Azimuth_deg);
            RangeDelta_m = Math.Abs(range_m - prevPoint.Range_m);
        }

        public BearingPoint3DTd(double azimuth_deg, double range_m, DateTime ts)
        {
            Azimuth_deg = azimuth_deg;
            Range_m = range_m;
            TimeStamp = ts;
            AzimuthDelta_deg = double.NaN;
            RangeDelta_m = double.NaN;
        }

        #endregion

        #region Methods

        public static double AngularDifference(double deg1, double deg2)
        {
            double diff = Math.Abs(deg1 - deg2) % 360.0;
            return diff > 180.0 ? 360.0 - diff : diff;
        }

        public override string ToString()
        {
            return string.Format("Az: {0:F03}°, Rng: {1:F03} m, ΔAz: {2:F03}°, ΔR: {3:F03} m",
                Azimuth_deg, Range_m, AzimuthDelta_deg, RangeDelta_m);
        }

        #endregion
    }

    /// <summary>
    /// ACHOD (Adaptive Concurrent Hypotheses Outlier Detector) для фильтрации пеленгов.
    /// Работает в пространстве абсолютных азимутов (heading + relative bearing).
    /// Использует две конкурирующие гипотезы траектории и автоматически выбирает лучшую.
    /// </summary>
    public class ACHODBearingFilter
    {
        #region Constants

        const double DEFAULT_MAX_AZIMUTH_RATE_DPS = 5.0;
        const double DEFAULT_MIN_SECTOR_WIDTH_DEG = 8.0;
        const double DEFAULT_SENSOR_NOISE_FACTOR = 3.0;
        const double DEFAULT_MAX_RANGE_CHANGE_M = 10.0;
        const int DEFAULT_FIFO_SIZE = 8;

        #endregion

        #region Properties

        List<BearingPoint3DTd> primarySide;    // всегда для принятых точек
        List<BearingPoint3DTd> secondarySide;  // всегда для отклонённых точек

        int pSideIdx = 1;  // 1 = primary активна, 2 = secondary активна

        double maxAzimuthRateDps;
        public double MaxAzimuthRateDps
        {
            get => maxAzimuthRateDps;
            set
            {
                if (value > 0.0)
                    maxAzimuthRateDps = value;
                else
                    throw new ArgumentOutOfRangeException("value", "Must be greater than 0 deg/s");
            }
        }

        double minSectorWidthDeg;
        public double MinSectorWidthDeg
        {
            get => minSectorWidthDeg;
            set
            {
                if (value > 0.0)
                    minSectorWidthDeg = value;
                else
                    throw new ArgumentOutOfRangeException("value", "Must be greater than 0 deg");
            }
        }

        double sensorNoiseFactor;
        public double SensorNoiseFactor
        {
            get => sensorNoiseFactor;
            set
            {
                if (value > 0.0)
                    sensorNoiseFactor = value;
                else
                    throw new ArgumentOutOfRangeException("value", "Must be greater than 0");
            }
        }

        double maxRangeChangeM;
        public double MaxRangeChangeM
        {
            get => maxRangeChangeM;
            set
            {
                if (value > 0.0)
                    maxRangeChangeM = value;
                else
                    throw new ArgumentOutOfRangeException("value", "Must be greater than 0 m");
            }
        }

        int fifoSize;
        public int FIFOSize
        {
            get => fifoSize;
            set
            {
                if (value < 2)
                    throw new ArgumentOutOfRangeException("value", "Must be >= 2");

                fifoSize = value;

                while (primarySide.Count > fifoSize)
                    primarySide.RemoveAt(0);
                while (secondarySide.Count > fifoSize)
                    secondarySide.RemoveAt(0);
            }
        }

        public int PrimarySideCount => primarySide.Count;
        public int SecondarySideCount => secondarySide.Count;
        public bool IsPrimaryActive => pSideIdx == 1;

        #endregion

        #region Constructor

        public ACHODBearingFilter()
            : this(DEFAULT_FIFO_SIZE, DEFAULT_MAX_AZIMUTH_RATE_DPS,
                   DEFAULT_MIN_SECTOR_WIDTH_DEG, DEFAULT_SENSOR_NOISE_FACTOR,
                   DEFAULT_MAX_RANGE_CHANGE_M)
        {
        }

        public ACHODBearingFilter(int fifoSize, double maxAzimuthRateDps,
            double minSectorWidthDeg, double sensorNoiseFactor, double maxRangeChangeM)
        {
            if (fifoSize < 2)
                throw new ArgumentOutOfRangeException("fifoSize", "Must be >= 2");

            this.fifoSize = fifoSize;
            MaxAzimuthRateDps = maxAzimuthRateDps;
            MinSectorWidthDeg = minSectorWidthDeg;
            SensorNoiseFactor = sensorNoiseFactor;
            MaxRangeChangeM = maxRangeChangeM;

            primarySide = new List<BearingPoint3DTd>();
            secondarySide = new List<BearingPoint3DTd>();
        }

        #endregion

        #region Private Methods

        private void AddPoint(List<BearingPoint3DTd> side, BearingPoint3DTd point)
        {
            if (side.Count + 1 > fifoSize)
                side.RemoveAt(0);

            side.Add(point);
        }

        private void PredictSector(List<BearingPoint3DTd> side, DateTime currentTime,
            out double center_deg, out double halfWidth_deg)
        {
            if (side.Count < 2)
            {
                center_deg = side[side.Count - 1].Azimuth_deg;
                halfWidth_deg = minSectorWidthDeg;
                return;
            }

            int nPts = Math.Min(4, side.Count);
            var rates = new List<double>();

            for (int k = side.Count - nPts + 1; k < side.Count; k++)
            {
                double dt = (side[k].TimeStamp - side[k - 1].TimeStamp).TotalSeconds;
                if (dt > 0)
                {
                    double rate = BearingPoint3DTd.AngularDifference(
                        side[k].Azimuth_deg, side[k - 1].Azimuth_deg) / dt;
                    rates.Add(rate);
                }
            }

            double estRate = 0.0;
            if (rates.Count > 0)
            {
                rates.Sort();
                estRate = rates[rates.Count / 2];
                estRate = Math.Min(estRate, maxAzimuthRateDps);
            }

            double dt2Pred = (currentTime - side[side.Count - 1].TimeStamp).TotalSeconds;
            center_deg = (side[side.Count - 1].Azimuth_deg + estRate * dt2Pred) % 360.0;
            if (center_deg < 0) center_deg += 360.0;

            double dynamicWidth = estRate * dt2Pred * 1.5;
            double noiseWidth = sensorNoiseFactor;
            halfWidth_deg = Math.Max(minSectorWidthDeg, dynamicWidth + noiseWidth);
        }

        private bool InSector(double azimuth_deg, double center_deg, double halfWidth_deg)
        {
            return BearingPoint3DTd.AngularDifference(azimuth_deg, center_deg) <= halfWidth_deg;
        }

        private double RangeConfidence(double range_m, List<BearingPoint3DTd> side)
        {
            if (side.Count < 3)
                return 0.5;

            int nPts = Math.Min(5, side.Count);
            double sum = 0.0;
            double sumSq = 0.0;

            for (int k = side.Count - nPts; k < side.Count; k++)
            {
                sum += side[k].Range_m;
                sumSq += side[k].Range_m * side[k].Range_m;
            }

            double mean = sum / nPts;
            double std = Math.Sqrt(sumSq / nPts - mean * mean);

            if (std < 0.5)
            {
                double z = Math.Abs(range_m - mean) / Math.Max(std, 0.1);
                return Math.Exp(-z * z / 2.0);
            }
            else
            {
                return 0.3;
            }
        }

        private double CalculateSigma(List<BearingPoint3DTd> side)
        {
            if (side.Count < 2)
                return double.MaxValue;

            double mean = 0.0;
            for (int i = 0; i < side.Count; i++)
                mean += side[i].AzimuthDelta_deg;
            mean /= side.Count;

            double sumSq = 0.0;
            for (int i = 0; i < side.Count; i++)
                sumSq += Math.Pow(side[i].AzimuthDelta_deg - mean, 2);

            return Math.Sqrt(sumSq);
        }

        #endregion

        #region Public Methods

        public void Reset()
        {
            primarySide.Clear();
            secondarySide.Clear();
            pSideIdx = 1;
        }

        /// <summary>
        /// Обработка одного измерения
        /// </summary>
        /// <param name="heading_deg">Курс антенны (магнитометр), градусы 0-360</param>
        /// <param name="relativeBearing_deg">Пеленг на маяк относительно антенны, градусы 0-360</param>
        /// <param name="range_m">Наклонная дальность, метры</param>
        /// <param name="timestamp">Метка времени</param>
        /// <param name="outAzimuth_deg">Отфильтрованный абсолютный азимут (NaN если отклонено)</param>
        /// <param name="outRange_m">Дальность (NaN если отклонено)</param>
        /// <returns>true если точка принята, false если отклонена как выброс</returns>
        public bool Process(double heading_deg, double relativeBearing_deg,
            double range_m, DateTime timestamp,
            out double outAzimuth_deg, out double outRange_m)
        {
            double azimuth_deg = (heading_deg + relativeBearing_deg) % 360.0;
            if (azimuth_deg < 0) azimuth_deg += 360.0;

            List<BearingPoint3DTd> activeSide = (pSideIdx == 1) ? primarySide : secondarySide;

            if (activeSide.Count == 0)
            {
                var newPoint = new BearingPoint3DTd(azimuth_deg, range_m, timestamp);
                AddPoint(activeSide, newPoint);

                outAzimuth_deg = azimuth_deg;
                outRange_m = range_m;
                return true;
            }

            double lastRange = activeSide[activeSide.Count - 1].Range_m;
            bool rangeOk = Math.Abs(range_m - lastRange) < maxRangeChangeM;

            PredictSector(activeSide, timestamp, out double sectorCenter, out double sectorHalfWidth);
            bool azimuthOk = InSector(azimuth_deg, sectorCenter, sectorHalfWidth);

            double rangeConf = RangeConfidence(range_m, activeSide);

            bool shouldAccept;

            if (azimuthOk)
            {
                shouldAccept = true;
            }
            else if (!azimuthOk && rangeOk && rangeConf > 0.7)
            {
                shouldAccept = false;
            }
            else if (!azimuthOk && !rangeOk)
            {
                shouldAccept = false;
            }
            else
            {
                double azimuthDev = BearingPoint3DTd.AngularDifference(azimuth_deg, sectorCenter);
                shouldAccept = azimuthDev < sectorHalfWidth * 1.5;
            }

            if (shouldAccept)
            {
                BearingPoint3DTd newPoint;
                if (primarySide.Count > 0)
                {
                    var lastPrimary = primarySide[primarySide.Count - 1];
                    newPoint = new BearingPoint3DTd(azimuth_deg, range_m, timestamp, lastPrimary);
                }
                else
                {
                    newPoint = new BearingPoint3DTd(azimuth_deg, range_m, timestamp);
                }

                AddPoint(primarySide, newPoint);

                outAzimuth_deg = azimuth_deg;
                outRange_m = range_m;
            }
            else
            {
                BearingPoint3DTd newPoint;
                if (secondarySide.Count > 0)
                {
                    var lastSecondary = secondarySide[secondarySide.Count - 1];
                    newPoint = new BearingPoint3DTd(azimuth_deg, range_m, timestamp, lastSecondary);
                }
                else
                {
                    newPoint = new BearingPoint3DTd(azimuth_deg, range_m, timestamp);
                }

                AddPoint(secondarySide, newPoint);

                outAzimuth_deg = double.NaN;
                outRange_m = double.NaN;
            }

            if (primarySide.Count == fifoSize && secondarySide.Count == fifoSize)
            {
                double sigmaPrimary = CalculateSigma(primarySide);
                double sigmaSecondary = CalculateSigma(secondarySide);

                if (pSideIdx == 1 && sigmaSecondary < sigmaPrimary)
                    pSideIdx = 2;
                else if (pSideIdx == 2 && sigmaPrimary < sigmaSecondary)
                    pSideIdx = 1;
            }

            return shouldAccept;
        }

        #endregion
    }
}