using System;
using System.Collections.Generic;
using System.Linq;

namespace UCNLNav
{
    #region Custom EventArgs

    public class TargetLocationUpdatedEventArgs : EventArgs
    {
        public DateTime TimeStamp { get; private set; }
        public GeoPoint3DE Location { get; private set; }

        public TargetLocationUpdatedEventArgs(GeoPoint3D loc, double rerr, DateTime ts)
            : this(loc.Latitude, loc.Longitude, loc.Depth, rerr, ts)
        {
        }

        public TargetLocationUpdatedEventArgs(double lat, double lon, double dpt, double rerr, DateTime ts)
        {
            Location = new GeoPoint3DE(lat, lon, dpt, rerr);
            TimeStamp = ts;
        }
    }

    public class TargetCourseUpdatedEventArgs : EventArgs
    {
        public double Course { get; private set; }
        public DateTime TimeStamp { get; private set; }

        public TargetCourseUpdatedEventArgs(double crs, DateTime ts)
        {
            Course = crs;
            TimeStamp = ts;
        }
    }

    public class TargetLocationUpdatedExEventArgs : EventArgs
    {
        public DateTime TimeStamp { get; private set; }
        public GeoPoint3DE Location { get; private set; }
        public double Course { get; private set; }

        public TargetLocationUpdatedExEventArgs(GeoPoint3D loc, double rerr, double course, DateTime ts)
            : this(loc.Latitude, loc.Longitude, loc.Depth, rerr, course, ts)
        {
        }

        public TargetLocationUpdatedExEventArgs(double lat, double lon, double dpt, double rerr, double course, DateTime ts)
        {
            Location = new GeoPoint3DE(lat, lon, dpt, rerr);
            Course = course;
            TimeStamp = ts;
        }
    }

    public class BaseQualityUpdatedEventArgs : EventArgs
    {
        #region Properties

        public TBAQuality TBAState { get; private set; }
        public double GDOP { get; private set; }
        public double PDOP { get; private set; }
        public double HDOP { get; private set; }
        public double VDOP { get; private set; }
        public double TDOP { get; private set; }
        public DOPState DopState { get; private set; }

        #endregion

        #region Constructor

        public BaseQualityUpdatedEventArgs(TBAQuality tbaState, double gdop, double pdop, double hdop, double vdop, double tdop, DOPState dState)
        {
            TBAState = tbaState;
            GDOP = gdop;
            PDOP = pdop;
            HDOP = hdop;
            VDOP = vdop;
            TDOP = tdop;
            DopState = dState;
        }

        #endregion
    }

    #endregion

    public class PCore2D<T> where T : GeoPoint3D
    {
        #region Properties

        public static readonly double SimplexSizeDefault = 10.0;
        public static readonly double RadialErrorThresholdDefault = 7.0;
        public static readonly int CourseEstimatorFifoSizeDefault = 8;

        double soundSpeed = 1500.0;
        public double SoundSpeed
        {
            get { return soundSpeed; }
            set
            {
                if (value > 0)
                    soundSpeed = value;
                else
                    throw new ArgumentOutOfRangeException();
            }
        }

        double simplexSize;
        public double SimplexSize
        {
            get { return simplexSize; }
            set
            {
                if (value > 0)
                    simplexSize = value;
                else
                    throw new ArgumentOutOfRangeException();
            }
        }

        double radialErrorThreshold;
        public double RadialErrorThreshold
        {
            get { return radialErrorThreshold; }
            set 
            {
                if (value > 0)
                    radialErrorThreshold = value;
                else
                    throw new ArgumentOutOfRangeException();
            }
        }
    
        public double TargetDepth
        {
            get { return targetLocation.Depth; }
            set { targetLocation.Depth = value; }
        }

        GeoPoint3D targetLocation;
        DateTime targetLocationTS;

        bool IsTargetLocation
        {
            get { return !double.IsNaN(targetLocation.Latitude) && !double.IsNaN(targetLocation.Longitude); }
        }
        
        public delegate void ExternalSolverDelegate(IEnumerable<T> basePoints, GeoPoint3D previousLocation, 
            out double lat_deg, out double lon_deg, out double rErr, out int it_cnt);
        public ExternalSolverDelegate ExternalSolver;

        public Ellipsoid referenceEllipsoid { get; private set; }

        CourseEstimatorLA2D crsEstimator;

        #endregion

        #region Constructor

        public PCore2D()
            : this(RadialErrorThresholdDefault, SimplexSizeDefault, Algorithms.WGS84Ellipsoid, CourseEstimatorFifoSizeDefault)
        {
        }

        public PCore2D(double rErrThreshold, double simplexSize, Ellipsoid refEllipsoid, int courseEstimatorFifoSize)
        {
            targetLocation = new GeoPoint3D(double.NaN, double.NaN, double.NaN);
            targetLocationTS = DateTime.MinValue;

            crsEstimator = new CourseEstimatorLA2D(courseEstimatorFifoSize);

            RadialErrorThreshold = rErrThreshold;
            SimplexSize = simplexSize;            
            referenceEllipsoid = refEllipsoid;
        }

        #endregion

        #region Methods

        #region Public

        public void ProcessBasePoints(IEnumerable<T> basePoints)
        {
            ProcessBasePoints(basePoints, targetLocation.Depth, DateTime.Now);
        }

        public void ProcessBasePoints(IEnumerable<T> basePoints, double depth)
        {
            ProcessBasePoints(basePoints, depth, DateTime.Now);
        }

        public void ProcessBasePoints(IEnumerable<T> basePoints, DateTime timeStamp)
        {
            ProcessBasePoints(basePoints, targetLocation.Depth, timeStamp);
        }

        public void ProcessBasePoints(IEnumerable<T> basePoints, double depth, DateTime timeStamp)
        {
            double lat_deg, lon_deg, rErr;
            int it_Cnt;            

            if (typeof(GeoPoint3DD).IsAssignableFrom(typeof(T)))
            {
                UCNLNav.Navigation.TOA_Locate2D(basePoints.Cast<GeoPoint3DD>().ToArray<GeoPoint3DD>(),
                    targetLocation.Latitude, targetLocation.Longitude, depth,
                    Algorithms.NLM_DEF_IT_LIMIT, Algorithms.NLM_DEF_PREC_THRLD, simplexSize,
                    referenceEllipsoid,
                    out lat_deg, out lon_deg, out rErr, out it_Cnt);
            }
            else if (typeof(GeoPoint3DT).IsAssignableFrom(typeof(T)))
            {
                UCNLNav.Navigation.TDOA_Locate2D(basePoints.Cast<GeoPoint3DT>().ToArray<GeoPoint3DT>(),
                    targetLocation.Latitude, targetLocation.Longitude, depth,
                    Algorithms.NLM_DEF_IT_LIMIT, Algorithms.NLM_DEF_PREC_THRLD, simplexSize,
                    referenceEllipsoid, soundSpeed,
                    out lat_deg, out lon_deg, out rErr, out it_Cnt);
            }
            else
            {
                if (ExternalSolver != null)
                    ExternalSolver(basePoints, targetLocation, out lat_deg, out lon_deg, out rErr, out it_Cnt);
                else
                    throw new NullReferenceException("ExternalSolver not defined");
            }

            if (rErr < radialErrorThreshold)
            {                               
                crsEstimator.AddPoint(new GeoPoint(lat_deg, lon_deg));

                if (crsEstimator.IsCourse)
                {
                    TargetCourseUpdatedHandler.Rise(this,
                        new TargetCourseUpdatedEventArgs(crsEstimator.Course_deg, timeStamp));
                }

                targetLocation.Latitude = lat_deg;
                targetLocation.Longitude = lon_deg;
                targetLocationTS = timeStamp;

                TargetLocationUpdatedHandler.Rise(this, 
                    new TargetLocationUpdatedEventArgs(targetLocation, rErr, timeStamp));

                TargetLocationUpdatedExHandler.Rise(this,
                    new TargetLocationUpdatedExEventArgs(targetLocation, rErr, crsEstimator.Course_deg, timeStamp));


                TBAQuality tbaState = Navigation.GetTBAState(Navigation.GetBasesMaxAngularGapDeg(basePoints, lat_deg, lon_deg));

                DOPState dopState = DOPState.Invalid;
                double gdop = double.NaN, pdop = double.NaN, hdop = double.NaN, vdop = double.NaN, tdop = double.NaN;

                GeoPoint3D tL = new GeoPoint3D(targetLocation.Latitude, targetLocation.Longitude, depth);
                if (Navigation.GetDOPs(basePoints, tL, Algorithms.WGS84Ellipsoid, out gdop, out pdop, out hdop, out vdop, out tdop))
                {
                    dopState = Navigation.GetDOPState(hdop);
                }

                BaseQualityUpdatedHandler.Rise(this,
                    new BaseQualityUpdatedEventArgs(tbaState, gdop, pdop, hdop, vdop, tdop, dopState));
            }
            else
            {
                RadialErrorExeedsThrehsoldEventHandler.Rise(this, new EventArgs());
            }
        }

        #endregion

        #endregion

        #region Events

        public EventHandler<TargetLocationUpdatedEventArgs> TargetLocationUpdatedHandler;
        public EventHandler<TargetCourseUpdatedEventArgs> TargetCourseUpdatedHandler;
        public EventHandler<TargetLocationUpdatedExEventArgs> TargetLocationUpdatedExHandler;
        public EventHandler<BaseQualityUpdatedEventArgs> BaseQualityUpdatedHandler;
        public EventHandler RadialErrorExeedsThrehsoldEventHandler;

        #endregion
    }
}
