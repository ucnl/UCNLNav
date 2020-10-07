using System;
using System.Collections.Generic;
using System.Linq;

namespace UCNLNav
{
    public class LocationResult : GeoPoint3D
    {
        #region Properties

        public double RadialError;
        public int Iterations;

        #endregion

        #region Constructor

        public LocationResult(double lat, double lon, double dpt, double radialError, int iterations)
            : base(lat, lon, dpt)
        {
            RadialError = radialError;
            Iterations = iterations;
        }

        #endregion

        #region Methods

        public override string ToString()
        {
            return string.Format("{0}, RER: {1:F03} m, ITR: {2}", base.ToString(), RadialError, Iterations);
        }

        #endregion
    }

    [Obsolete]
    public class PositioningCore2D<T> where T : GeoPoint3D
    {
        #region Properties

        public double TargetDepth
        {
            get { return PreviousLocation.Depth; }
            set { PreviousLocation.Depth = value; }
        }

        double soundSpeed = 1500;
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

        double simplexSize = 10;
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

        double radialErrorThreshold = 7.0;
        public double RadialErrorThreshold
        {
            get { return radialErrorThreshold; }
            set { radialErrorThreshold = value; }
        }

        public GeoPoint3D PreviousLocation;        

        public delegate void ExternalSolverDelegate(IEnumerable<T> basePoints, GeoPoint3D previousLocation, 
            out double lat_deg, out double lon_deg, out double rErr, out int it_cnt);
        public ExternalSolverDelegate ExternalSolver;

        Ellipsoid referenceEllipsoid;

        #endregion

        #region Constructor

        public PositioningCore2D(double rErrThreshold, double simplexSize, Ellipsoid refEllipsoid)
        {
            PreviousLocation = new GeoPoint3D(double.NaN, double.NaN, 0);
            RadialErrorThreshold = rErrThreshold;
            SimplexSize = simplexSize;
            TargetDepth = 0;
            referenceEllipsoid = refEllipsoid;            
        }

        #endregion

        #region Methods

        #region Public

        public LocationResult Locate(IEnumerable<T> basePoints)
        {
            double lat_deg, lon_deg, radialError;
            int it_Cnt;

            if (typeof(GeoPoint3DD).IsAssignableFrom(typeof(T)))
            {
                UCNLNav.Navigation.TOA_Locate2D(basePoints.Cast<GeoPoint3DD>().ToArray<GeoPoint3DD>(),
                    PreviousLocation.Latitude, PreviousLocation.Longitude, PreviousLocation.Depth,
                    Algorithms.NLM_DEF_IT_LIMIT, Algorithms.NLM_DEF_PREC_THRLD, simplexSize,
                    referenceEllipsoid,
                    out lat_deg, out lon_deg, out radialError, out it_Cnt);
            }
            else if (typeof(GeoPoint3DT).IsAssignableFrom(typeof(T)))
            {
                UCNLNav.Navigation.TDOA_Locate2D(basePoints.Cast<GeoPoint3DT>().ToArray<GeoPoint3DT>(),
                    PreviousLocation.Latitude, PreviousLocation.Longitude, PreviousLocation.Depth,
                    Algorithms.NLM_DEF_IT_LIMIT, Algorithms.NLM_DEF_PREC_THRLD, simplexSize,
                    referenceEllipsoid, soundSpeed,
                    out lat_deg, out lon_deg, out radialError, out it_Cnt);
            }
            else
            {
                if (ExternalSolver != null)                    
                    ExternalSolver(basePoints, PreviousLocation, out lat_deg, out lon_deg, out radialError, out it_Cnt);
                else
                    throw new NullReferenceException("ExternalSolver not defined");
            }

            if (radialError < radialErrorThreshold)
            {
                PreviousLocation.Latitude = lat_deg;
                PreviousLocation.Longitude = lon_deg;
                PreviousLocationUpdatedEventHandler.Rise(this, new EventArgs());
            }

            return new LocationResult(lat_deg, lon_deg, PreviousLocation.Depth, radialError, it_Cnt);
        }

        #endregion

        #endregion

        #region Events

        public EventHandler PreviousLocationUpdatedEventHandler;

        #endregion
    }
}
