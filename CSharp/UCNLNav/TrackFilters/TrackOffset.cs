using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace UCNLNav.TrackFilters
{
    public class TrackOffset : ITrackFilter
    {
        #region Properties

        public double Latitude_offset_m { get; private set; }
        public double Longitude_offset_m { get; private set; }

        #endregion

        #region Constructor

        public TrackOffset(double xoffset_m, double yoffset_m)
        {
            Latitude_offset_m = yoffset_m;
            Longitude_offset_m = xoffset_m;
        }

        #endregion

        #region Methods

        public void Reset()
        {
            // nothing to do
        }

        public bool Process(double inLat_rad, double inLon_rad, double inDpt_m, DateTime inTS,
            out double outLat_rad, out double outLon_rad, out double outDpt_m, out DateTime outTS)
        {
            outDpt_m = inDpt_m;
            outTS = inTS;

            Algorithms.GeopointOffsetByDeltas_WGS84(inLat_rad, inLon_rad, Latitude_offset_m, Longitude_offset_m,
                out outLat_rad, out outLon_rad);

            return true;

        }


        #endregion
    }
}
