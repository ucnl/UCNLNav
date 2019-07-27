
using System;
namespace UCNLNav
{
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
            return string.Format("{0:F06}, {1:F06}", Latitude, Longitude);
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
            return string.Format("{0}, {1:F03}", base.ToString(), Depth);
        }

        #endregion
    }

    public class GeoPoint3DT : GeoPoint3D
    {
        #region Properties

        public DateTime TimeStamp;

        #endregion

        #region Constructor

        public GeoPoint3DT(double lat, double lon, double dpt, DateTime timeStamp)
            : base(lat, lon, dpt)
        {
            TimeStamp = timeStamp;
        }

        #endregion

        #region Methods

        public override string ToString()
        {
            return string.Format("{0}, {1}", base.ToString(), TimeStamp.ToShortTimeString());
        }

        #endregion
    }
}
