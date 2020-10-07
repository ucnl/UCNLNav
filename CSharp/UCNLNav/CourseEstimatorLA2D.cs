using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace UCNLNav
{
    public class CourseEstimatorLA2D
    {
        #region Properties

        List<MPoint> fifo;

        int fifoSize = 16;
        public int FIFOSIze
        {
            get { return fifoSize; }
            private set 
            {
                if (value < 1)
                    throw new ArgumentOutOfRangeException("value", "Sould be greater than 1");
                else
                    fifoSize = value;                
            }
        }        

        bool initialized = false;

        public bool IsCourse
        {
            get { return !double.IsNaN(Course_deg); }
        }

        public double Course_deg { get; private set; }

        double anchor_lat_rad;
        double anchor_lon_rad;

        #endregion

        #region Constructor

        public CourseEstimatorLA2D(int fifoSize)
        {
            FIFOSIze = fifoSize;
            fifo = new List<MPoint>(fifoSize);
        }

        #endregion

        #region Methods

        private void Clear()
        {
            Course_deg = double.NaN;
            initialized = false;
            fifo.Clear();
        }

        private void AppendPoint(double x, double y)
        {
            if (fifo.Count + 1 >= fifoSize)
                fifo.RemoveAt(0);

            fifo.Add(new MPoint(x, y));

            if (fifo.Count >= 2)
            {
                double k = 0, b = 0, a = 0;
                Algorithms.LinearApproxXY(fifo, out k, out b, out a);

                if (!double.IsNaN(a))
                {
                    a = Algorithms.Rad2Deg(a);
                    a = 90 - a + 180;
                    if (a < 0)
                        a += 360;
                    Course_deg = Algorithms.Wrap360(a);
                }
            }
        }

        public void AddPoint(GeoPoint newPoint)
        {
            double dx = 0, dy = 0;

            if (initialized)
            {                
                Algorithms.GetDeltasByGeopoints_WGS84(anchor_lat_rad, anchor_lon_rad,
                    Algorithms.Deg2Rad(newPoint.Latitude), Algorithms.Deg2Rad(newPoint.Longitude),
                    out dy, out dx);
            }
            else
            {
                anchor_lat_rad = Algorithms.Deg2Rad(newPoint.Latitude);
                anchor_lon_rad = Algorithms.Deg2Rad(newPoint.Longitude);
                initialized = true;                    
            }

            AppendPoint(dx, dy);
        }

        #endregion
    }
}
