
using System;
using UCNLNav;
namespace UCNLNav_Tests
{
    class Program
    {
        static void Main(string[] args)
        {
            int bases_number = 4;  // number of base points
            double start_bases_z = 1.5;  // initial z-coordinate of base points
            double base_z_step = 5; // base z-coordinate step

            // actual target location
            double actual_target_lat_deg = 44.12345678; // singed degrees
            double actual_target_lon_deg = 48.12345678; // signed degrees
            double actual_target_z = 25;                       // meters
            
            double start_dst_projection = 500;          // start distance from target, meters
            double dst_inc_step = 50;                   // distance increment
                        
            // azimuth step
            double az_step = Math.PI * 2 / bases_number;
            double actual_target_lat_rad = Algorithms.Deg2Rad(actual_target_lat_deg);
            double actual_target_lon_rad = Algorithms.Deg2Rad(actual_target_lon_deg);

            // signal propagation speed
            double velocity = 1450.0; // m/s

            double az = 0;
            double base_lat_rad = 0;
            double base_lon_rad = 0;
            double rev_az_rad = 0;
            int it_cnt = 0;
            double dst_projection;
            double slant_range;
            double dZ;// = actual_target_z - start_bases_z;
            bool res = false;
            double dst_proj = 0;
            double fwd_az_rad = 0;
            double target_lat = 0;
            double target_lon = 0;
            double radialError = 0;
            double dst_error = 0;
            double target_z = 0;
            double base_z;

            GeoPoint3DD[] toaBases = new GeoPoint3DD[bases_number];
            GeoPoint3DT[] tdoaBases = new GeoPoint3DT[bases_number];


            Console.WriteLine();
            Console.WriteLine("* Testing Vincenty equations...");
            Console.WriteLine("Generating base points...");
            for (int baseIdx = 0; baseIdx < bases_number; baseIdx++)
            {
                dst_projection = start_dst_projection + dst_inc_step * baseIdx;
                res = Algorithms.VincentyDirect(actual_target_lat_rad, actual_target_lon_rad,
                                          az, dst_projection,
                                          Algorithms.WGS84Ellipsoid,
                                          Algorithms.VNC_DEF_EPSILON, Algorithms.VNC_DEF_IT_LIMIT,
                                          out base_lat_rad, out base_lon_rad, out rev_az_rad, out it_cnt);

                base_z = start_bases_z + base_z_step * baseIdx;
                dZ = actual_target_z - base_z;
                slant_range = Math.Sqrt(dst_projection * dst_projection + dZ * dZ);

                toaBases[baseIdx] = new GeoPoint3DD(Algorithms.Rad2Deg(base_lat_rad), Algorithms.Rad2Deg(base_lon_rad), base_z, slant_range);
                tdoaBases[baseIdx] = new GeoPoint3DT(toaBases[baseIdx].Latitude, toaBases[baseIdx].Longitude, toaBases[baseIdx].Depth, toaBases[baseIdx].SlantRange / velocity);

                Console.WriteLine(string.Format("Base #{0}", baseIdx));
                Console.WriteLine(string.Format("Distance projection: {0:F03} m", dst_projection));
                Console.WriteLine(string.Format("Azimuth: {0:F03}°", Algorithms.Rad2Deg(az)));

                Console.WriteLine(string.Format("Vincenty direct ({0}): LAT: {1:F07}°, LON: {2:F07}°, Iterations: {3}",
                    res,
                    toaBases[baseIdx].Latitude,
                    toaBases[baseIdx].Longitude,
                    it_cnt));
                
                it_cnt = 0;
                res = Algorithms.VincentyInverse(actual_target_lat_rad, actual_target_lon_rad,
                                           base_lat_rad, base_lon_rad,
                                           Algorithms.WGS84Ellipsoid,
                                           Algorithms.VNC_DEF_EPSILON, Algorithms.VNC_DEF_IT_LIMIT,
                                           out dst_proj, out fwd_az_rad, out rev_az_rad, out it_cnt);

                Console.WriteLine("Checking with Vincenty inverse:");
                Console.WriteLine(string.Format("Vincenty inverse ({0}): DST: {1:F03} m ({2:F03}), AZM: {3:F03}° ({4:F03})°, Iterations: {5}",
                    res,
                    dst_proj, dst_projection,
                    Algorithms.Rad2Deg(fwd_az_rad), Algorithms.Rad2Deg(az),                    
                    it_cnt));

                az += az_step;
            }

            Console.WriteLine();
            Console.WriteLine("Press a key to start TOA 2D test...");
            Console.ReadKey();

            Console.WriteLine();
            Console.WriteLine("* Solving a TOA problem with TOA_Locate2D...");
            Navigation.TOA_Locate2D(toaBases,
                                    double.NaN, double.NaN,
                                    actual_target_z,
                                    Algorithms.NLM_DEF_IT_LIMIT, Algorithms.NLM_DEF_PREC_THRLD, 10.0,
                                    Algorithms.WGS84Ellipsoid,
                                    out target_lat, out target_lon, out radialError, out it_cnt);

            Console.WriteLine("LAT: {0:F07}° (Actual {1:F07}°), LON: {2:F07}° (Actual {3:F07}°), Estimated radial error: {4:F03} m, Iterations: {5}",
                target_lat, actual_target_lat_deg, target_lon, actual_target_lon_deg, radialError, it_cnt);

            Algorithms.VincentyInverse(actual_target_lat_rad, actual_target_lon_rad, 
                                       Algorithms.Deg2Rad(target_lat), Algorithms.Deg2Rad(target_lon),
                                       Algorithms.WGS84Ellipsoid, Algorithms.VNC_DEF_EPSILON, Algorithms.VNC_DEF_IT_LIMIT,
                                       out dst_error, out fwd_az_rad, out rev_az_rad, out it_cnt);
            Console.WriteLine(string.Format("Distance on the ellipsoid between actual target location and calculated: {0:F03} m", dst_error));            

            
            Console.WriteLine();
            Console.WriteLine("Press a key to start TDOA test...");
            Console.ReadKey();

            Console.WriteLine();
            Console.WriteLine("* Solving a TDOA problem with TDOA_Locate2D...");
                  
            for (int i = 0; i < tdoaBases.Length; i++)
            {
                tdoaBases[i] = new GeoPoint3DT(toaBases[i].Latitude, toaBases[i].Longitude, toaBases[i].Depth, toaBases[i].SlantRange / velocity);
            }

            Navigation.TDOA_Locate2D(tdoaBases,
                                     double.NaN, double.NaN,
                                     actual_target_z,
                                     Algorithms.NLM_DEF_IT_LIMIT, Algorithms.NLM_DEF_PREC_THRLD, 10,
                                     Algorithms.WGS84Ellipsoid,
                                     velocity,
                                     out target_lat, out target_lon, out radialError, out it_cnt);

            Console.WriteLine("LAT: {0:F07}° (Actual {1:F07}°), LON: {2:F07}° (Actual {3:F07}°), Estimated radial error: {4:F03} m, Iterations: {5}",
               target_lat, actual_target_lat_deg, target_lon, actual_target_lon_deg, radialError, it_cnt);

            Algorithms.VincentyInverse(actual_target_lat_rad, actual_target_lon_rad,
                                       Algorithms.Deg2Rad(target_lat), Algorithms.Deg2Rad(target_lon),
                                       Algorithms.WGS84Ellipsoid, Algorithms.VNC_DEF_EPSILON, Algorithms.VNC_DEF_IT_LIMIT,
                                       out dst_error, out fwd_az_rad, out rev_az_rad, out it_cnt);
            Console.WriteLine(string.Format("Distance on the ellipsoid between actual target location and calculated: {0:F03} m", dst_error));



            Console.WriteLine();
            Console.WriteLine("Press a key to start TOA 3D test...");
            Console.ReadKey();

            Console.WriteLine();
            Console.WriteLine("* Solving a TOA problem with TOA_Locate3D...");
            Navigation.TOA_Locate3D(toaBases,
                                    double.NaN, double.NaN, double.NaN,
                                    Algorithms.NLM_DEF_IT_LIMIT, Algorithms.NLM_DEF_PREC_THRLD, 10.0,
                                    Algorithms.WGS84Ellipsoid,
                                    out target_lat, out target_lon, out target_z, out radialError, out it_cnt);

            Console.WriteLine("LAT: {0:F07}° (Actual {1:F07}°), LON: {2:F07}° (Actual {3:F07}°), DPT: {4:F03} m ({5:F03}), Estimated radial error: {6:F03} m, Iterations: {7}",
                target_lat, actual_target_lat_deg, target_lon, actual_target_lon_deg, target_z, actual_target_z, radialError, it_cnt);

            Algorithms.VincentyInverse(actual_target_lat_rad, actual_target_lon_rad,
                                       Algorithms.Deg2Rad(target_lat), Algorithms.Deg2Rad(target_lon),
                                       Algorithms.WGS84Ellipsoid, Algorithms.VNC_DEF_EPSILON, Algorithms.VNC_DEF_IT_LIMIT,
                                       out dst_error, out fwd_az_rad, out rev_az_rad, out it_cnt);
            Console.WriteLine(string.Format("Distance on the ellipsoid between actual target location and calculated: {0:F03} m", dst_error));




            Console.WriteLine();
            Console.WriteLine("Press a key to exit...");
            Console.ReadKey();            
        }
    }
}
