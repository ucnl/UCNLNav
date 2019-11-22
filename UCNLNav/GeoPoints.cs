
using System;
namespace UCNLNav
{
    /// <summary>
    /// Metric point 2D
    /// </summary>
    public class MPoint
    {
        #region Properties

        public double X;
        public double Y;        

        #endregion

        #region Constructor

        public MPoint()
            : this(0, 0)
        {
        }

        public MPoint(double x, double y)
        {
            X = x;
            Y = y;
        }

        #endregion

        #region Methods

        public override string ToString()
        {
            return string.Format("X: {0:F06}, Y: {1:F06}", X, Y);
        }

        #endregion
    }

    /// <summary>
    /// Metric point 3D
    /// </summary>
    public class MPoint3D : MPoint
    {
        #region Properties

        public double Z;

        #endregion

        #region Constructor

        public MPoint3D(double x, double y, double z)
            : base(x, y)
        {
            Z = z;
        }

        #endregion

        #region Methods

        public override string ToString()
        {
            return string.Format("{0}, Z: {1:F06}", base.ToString(), Z);
        }

        #endregion
    }

    public class GeoPoint
    {
        #region Properties

        public double Latitude;
        public double Longitude;        

        #endregion

        #region Constructor

        public GeoPoint()
            : this(0, 0)
        {
        }

        public GeoPoint(double lat, double lon)
        {
            Latitude = lat;
            Longitude = lon;
        }

        #endregion

        #region Methods

        public override string ToString()
        {
            return string.Format("LAT: {0:F06}°, LON: {1:F06}°", Latitude, Longitude);
        }

        #endregion
    }

    public class GeoPoint3D : GeoPoint
    {
        #region Properties

        public double Depth;

        #endregion

        #region Constructor

        public GeoPoint3D(double lat, double lon, double dpt)
            : base(lat, lon)
        {
            Depth = dpt;
        }

        #endregion

        #region Methods

        public override string ToString()
        {
            return string.Format("{0}, DPT: {1:F03} m", base.ToString(), Depth);
        }

        #endregion
    }

    public class GeoPoint3DD : GeoPoint3D
    {
        #region Properties

        public double SlantRange;

        #endregion

        #region Constructor

        public GeoPoint3DD(double lat, double lon, double dpt, double slantRange)
            : base(lat, lon, dpt)
        {
            SlantRange = slantRange;
        }

        #endregion

        #region Methods

        public override string ToString()
        {
            return string.Format("{0}, SRNG: {1:F03} m", base.ToString(), SlantRange);
        }

        #endregion
    }

    public class GeoPoint3DT : GeoPoint3D
    {
        #region Properties

        public double TOASec;

        #endregion

        #region Constructor

        public GeoPoint3DT(double lat, double lon, double dpt, double toaSec)
            : base(lat, lon, dpt)
        {
            TOASec = toaSec;
        }

        #endregion

        #region Methods

        public override string ToString()
        {
            return string.Format("{0}, TOA: {1:F06} sec", base.ToString(), TOASec);
        }

        #endregion
    }

    public class GeoPoint3DTm : GeoPoint3D
    {
        #region Properties

        public DateTime TimeStamp;

        #endregion

        #region Constructor

        public GeoPoint3DTm(double lat, double lon, double dpt, DateTime timeStamp)
            : base(lat, lon, dpt)
        {
            TimeStamp = timeStamp;
        }

        #endregion

        #region Methods

        public override string ToString()
        {
            return string.Format("{0}, TS: {1}", base.ToString(), TimeStamp.ToShortTimeString());
        }

        #endregion
    }

    public struct TOABasePoint
    {
        public double X; 
        public double Y;
        public double Z;
        public double D; // Distance or pseudorange

        public TOABasePoint(double x, double y, double z, double d)
        {
            X = x;
            Y = y;
            Z = z;
            D = d;
        }
    }   

    public struct TDOABaseline
    {
        public double X1;
        public double Y1;
        public double Z1;
        public double X2;
        public double Y2;
        public double Z2;
        public double PRD; // Pseudorange difference

        public TDOABaseline(double x1, double y1, double z1, double x2, double y2, double z2, double prd)
        {
            X1 = x1;
            Y1 = y1;
            Z1 = z1;
            X2 = x2;
            Y2 = y2;
            Z2 = z2;
            PRD = prd;
        }
    }
   
}
