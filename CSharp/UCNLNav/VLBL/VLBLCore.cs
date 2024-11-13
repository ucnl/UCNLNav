using System;
using System.Collections.Generic;
using System.Threading;

namespace UCNLNav.VLBL
{
    #region Custom EventArgs

    public class BaseUpdatedEventArgs : EventArgs
    {
        #region Properties

        public GeoPoint ReferencePoint { get; private set; }
        public List<IVLBLMeasurement> BasePoints { get; private set; }

        #endregion

        #region Constructor

        public BaseUpdatedEventArgs(GeoPoint referencePoint, IEnumerable<IVLBLMeasurement> basePoints)
        {
            ReferencePoint = new GeoPoint(referencePoint.Latitude, referencePoint.Longitude);
            BasePoints = new List<IVLBLMeasurement>(basePoints);
        }

        #endregion
    }

    public class TargetPositionUpdatedEventArgs : EventArgs
    {
        #region Properties

        public GeoPoint3D TargetPosition { get; private set; }
        public double RadialError { get; private set; }
        public int Iterations { get; private set; }

        #endregion

        #region Constructor

        public TargetPositionUpdatedEventArgs(double lat, double lon, double dpt, double radialError, int it_cnt)
        {
            TargetPosition = new GeoPoint3D(lat, lon, dpt);
            RadialError = radialError;
            Iterations = it_cnt;
        }

        #endregion
    }
    
    #endregion

    public class VLBLCore<T> where T : GeoPoint3D, IVLBLMeasurement
    {
        #region Properties

        public bool IsBaseExists { get { return measurements.IsBaseExists; } }
        public int Capacity { get { return measurements.Capacity; } }
        public int BaseSize { get { return measurements.BaseSize; } }

        public double TargetDepth
        {
            get { return positioningCore.PreviousLocation.Depth; }
            set { positioningCore.PreviousLocation.Depth = value; }
        }
        public GeoPoint ReferencePoint { get { return measurements.RefPoint; } }

        public double SoundSpeed
        {
            get { return positioningCore.SoundSpeed; }
            set { positioningCore.SoundSpeed = value; }
        }

        public double SimplexSize
        {
            get { return positioningCore.SimplexSize; }
            set { positioningCore.SimplexSize = value; }
        }

        public double RadialErrorThreshold
        {
            get { return positioningCore.RadialErrorThreshold; }
            set { positioningCore.RadialErrorThreshold = value; }
        }

        VLBLMeasurements<T> measurements;

        PositioningCore2D<T> positioningCore;

        int tLock = 0;

        #endregion

        #region Constructor

        public VLBLCore(int capacity, int baseSize, double rErrThreshold, double simplexSize, Ellipsoid refEllipsoid)
        {
            measurements = new VLBLMeasurements<T>(capacity, baseSize);
            
            positioningCore = new PositioningCore2D<T>(rErrThreshold, simplexSize, refEllipsoid);
            positioningCore.PreviousLocationUpdatedEventHandler += (o, e) =>
            {
                measurements.SetRefPoint(positioningCore.PreviousLocation);
                ReferencePointUpdatedEventHandler.Rise(this, new EventArgs());
            };
            TargetDepth = 0;
        }

        #endregion

        #region Methods

        #region Public

        public void DiscardPrevious()
        {
            while (Interlocked.CompareExchange(ref tLock, 1, 0) != 0)
                Thread.SpinWait(1);

            positioningCore.PreviousLocation = new GeoPoint3D(double.NaN, double.NaN, positioningCore.PreviousLocation.Depth);

            Interlocked.Decrement(ref tLock);
        }

        public void ResetMeasurements()
        {
            while (Interlocked.CompareExchange(ref tLock, 1, 0) != 0)
                Thread.SpinWait(1);

            measurements.Clear();

            Interlocked.Decrement(ref tLock);
        }

        public void AddMeasurement(T newMeasurement)
        {
            while (Interlocked.CompareExchange(ref tLock, 1, 0) != 0)
                Thread.SpinWait(1);

            measurements.Add(newMeasurement);
            if (measurements.IsBaseExists)
            {
                var basePoints = measurements.GetBase();
                BaseUpdatedEventHandler.Rise(this, new BaseUpdatedEventArgs(ReferencePoint, basePoints));

                var result = positioningCore.Locate(basePoints);
                TargetPositionUpdatedEventHandler.Rise(this,
                    new TargetPositionUpdatedEventArgs(result.Latitude, result.Longitude, result.Depth,
                    result.RadialError, result.Iterations));
            }

            Interlocked.Decrement(ref tLock);
        }

        #endregion
        
        #endregion

        #region Events

        public EventHandler<BaseUpdatedEventArgs> BaseUpdatedEventHandler;
        public EventHandler ReferencePointUpdatedEventHandler;
        public EventHandler<TargetPositionUpdatedEventArgs> TargetPositionUpdatedEventHandler;

        #endregion
    }
}
