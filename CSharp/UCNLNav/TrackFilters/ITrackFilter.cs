using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace UCNLNav.TrackFilters
{
    interface ITrackFilter
    {
        void Reset();
        bool Process(double inLat_rad, double inLon_rad, double inDpt_m, DateTime inTS,
            out double outLat_rad, out double outLon_rad, out double outDpt_m, out DateTime outTS);
    }

    interface ITrackCartesianFilter
    {
        void Reset();
        bool Process(double x, double y, double z, DateTime inTS,
            out double flt_x, out double flt_y, out double flt_z, out DateTime outTS);
    }
}
