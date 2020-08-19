using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace UCNLNav
{
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

    public class TargetCourseAndSpeedUpdatedEventArgs : EventArgs
    {
        public double Course { get; private set; }
        public double Speed { get; private set; }
        public DateTime TimeStamp { get; private set; }

        public bool IsSpeedValid
        {
            get { return !double.IsNaN(Speed); }
        }

        public TargetCourseAndSpeedUpdatedEventArgs(double crs, double spd, DateTime ts)
        {
            Course = crs;
            Speed = spd;
            TimeStamp = ts;
        }
    }

    public class PCore2D<T> where T : GeoPoint3D
    {
        #region Properties

        public static readonly double SimplexSizeDefault = 10.0;
        public static readonly double RadialErrorThresholdDefault = 7.0;
        public static readonly int MaxIntervalToSpeedEstimationSecDefault = 5;

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

        int maxIntervalToSpeedEstimationSec;
        public int MaxIntervalToSpeedEstimationSec
        {
            get { return maxIntervalToSpeedEstimationSec; }
            set
            {
                if (value > 0)
                    maxIntervalToSpeedEstimationSec = value;
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

        #endregion

        #region Constructor

        public PCore2D()
            : this(RadialErrorThresholdDefault, SimplexSizeDefault, Algorithms.WGS84Ellipsoid, MaxIntervalToSpeedEstimationSecDefault)
        {
        }

        public PCore2D(double rErrThreshold, double simplexSize, Ellipsoid refEllipsoid, int maxIntToSpeedEstimationSec)
        {
            targetLocation = new GeoPoint3D(double.NaN, double.NaN, double.NaN);
            targetLocationTS = DateTime.MinValue;

            RadialErrorThreshold = rErrThreshold;
            SimplexSize = simplexSize;            
            referenceEllipsoid = refEllipsoid;
            MaxIntervalToSpeedEstimationSec = maxIntToSpeedEstimationSec;
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
                if (IsTargetLocation)
                {
                    var duration = timeStamp.Subtract(targetLocationTS);
                    if (duration.TotalSeconds <= maxIntervalToSpeedEstimationSec)
                    {
                        var prevLocRad = GeoPoint3D.ToRad(targetLocation);
                        double currLocLatRad = Algorithms.Deg2Rad(lat_deg);
                        double currLocLonRad = Algorithms.Deg2Rad(lon_deg);
                        double d_passed_m = 0, course_fwd_rad = 0, course_rev_rad = 0;
                        int its = 0;

                        Algorithms.VincentyInverse(prevLocRad.Latitude, prevLocRad.Longitude,
                            currLocLatRad, currLocLonRad, referenceEllipsoid,
                            Algorithms.VNC_DEF_EPSILON, Algorithms.VNC_DEF_IT_LIMIT,
                            out d_passed_m, out course_fwd_rad, out course_rev_rad, out its);

                        double estTargetCourse = Algorithms.Rad2Deg(course_fwd_rad);
                        double estTargetSpeed = d_passed_m / duration.TotalSeconds;
                    
                        TargetCourseSpeedAndCourseUpdatedHandler.Rise(this, 
                            new TargetCourseAndSpeedUpdatedEventArgs(estTargetCourse, estTargetSpeed, timeStamp));
                    }
                }

                targetLocation.Latitude = lat_deg;
                targetLocation.Longitude = lon_deg;
                targetLocationTS = timeStamp;

                TargetLocationUpdatedHandler.Rise(this, 
                    new TargetLocationUpdatedEventArgs(targetLocation, rErr, timeStamp));
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
        public EventHandler<TargetCourseAndSpeedUpdatedEventArgs> TargetCourseSpeedAndCourseUpdatedHandler;
        public EventHandler RadialErrorExeedsThrehsoldEventHandler;

        #endregion
    }
}
