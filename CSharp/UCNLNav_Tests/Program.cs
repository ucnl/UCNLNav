
using System;
using System.Collections.Generic;
using UCNLNav;
using UCNLNav.VLBL;

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
            double actual_target_z = 25;                // meters
            
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


            #region Vincenty equations testing

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

            #endregion

            Console.WriteLine("\r\nPress a key to start TOA 2D test...");
            Console.ReadKey();

            #region TOA_Locate2D testing

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

            #endregion

            Console.WriteLine("\r\nPress a key to start TDOA test...");
            Console.ReadKey();

            #region TDOA_Locate2D testing

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

            #endregion

            Console.WriteLine("\r\nPress a key to start TOA 3D test...");
            Console.ReadKey();

            #region TOA_Locate3D testing

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

            #endregion

            Console.WriteLine("\r\nPress a key to start TDOA 3D test...");
            Console.ReadKey();

            #region TDOA_Locate3D testing

            Console.WriteLine();
            Console.WriteLine("* Solving a TDOA problem with TOA_Locate3D...");

            Navigation.TDOA_Locate3D(tdoaBases,
                                     double.NaN, double.NaN,
                                     double.NaN,
                                     Algorithms.NLM_DEF_IT_LIMIT, Algorithms.NLM_DEF_PREC_THRLD, 10,
                                     Algorithms.WGS84Ellipsoid,
                                     velocity,
                                     out target_lat, out target_lon, out target_z, out radialError, out it_cnt);

            Console.WriteLine("LAT: {0:F07}° (Actual {1:F07}°), LON: {2:F07}° (Actual {3:F07}°), DPT: {4:F02} m, Estimated radial error: {4:F03} m, Iterations: {5}",
               target_lat, actual_target_lat_deg, target_lon, actual_target_lon_deg, radialError, it_cnt);

            Algorithms.VincentyInverse(actual_target_lat_rad, actual_target_lon_rad,
                                       Algorithms.Deg2Rad(target_lat), Algorithms.Deg2Rad(target_lon),
                                       Algorithms.WGS84Ellipsoid, Algorithms.VNC_DEF_EPSILON, Algorithms.VNC_DEF_IT_LIMIT,
                                       out dst_error, out fwd_az_rad, out rev_az_rad, out it_cnt);
            Console.WriteLine(string.Format("Distance on the ellipsoid between actual target location and calculated: {0:F03} m", dst_error));

            #endregion


            Console.WriteLine("\r\nPress a key to start VLBL tests...");
            Console.ReadKey();

            #region Testing VLBL
                        
            int vlbl_fifoSize = 256;
            int vlbl_baseSize = 3;
            double vlbl_radialErrorThreshold = 7;
            double vlbl_simplexSize = 10;
            GeoPoint3D actual_target_location = new GeoPoint3D(actual_target_lat_deg, actual_target_lon_deg, actual_target_z);

            Console.WriteLine(string.Format("Number of base points: {0}", vlbl_baseSize));
            Console.WriteLine(string.Format("VLBL FIFO size: {0}", vlbl_fifoSize));

            VLBLCore<VLBLTOAMeasurement> vlblCore = new VLBLCore<VLBLTOAMeasurement>(vlbl_fifoSize, vlbl_baseSize, vlbl_radialErrorThreshold,
                vlbl_simplexSize, Algorithms.WGS84Ellipsoid);
            vlblCore.BaseUpdatedEventHandler += (o, e) =>
                {
                    Console.WriteLine(string.Format("Base updated, refPoint: {0}", e.ReferencePoint));
                    Console.WriteLine("Base points:");
                    foreach (var basePoint in e.BasePoints)
                    {
                        Console.WriteLine(basePoint.ToString());
                    }
                };

            vlblCore.TargetDepth = actual_target_z;

            vlblCore.ReferencePointUpdatedEventHandler += (o, e) =>
                {
                    Console.WriteLine(string.Format("Reference point updated: {0}", vlblCore.ReferencePoint));
                };

            vlblCore.TargetPositionUpdatedEventHandler += (o, e) =>
                {
                    Console.WriteLine(string.Format("Target position updated: {0}, RERR: {1:F03} m, Iterations: {2}", e.TargetPosition, e.RadialError, e.Iterations));
                    Console.WriteLine(string.Format("Actual target position: {0}", actual_target_location));

                    Algorithms.VincentyInverse(actual_target_lat_rad, actual_target_lon_rad,
                                       Algorithms.Deg2Rad(e.TargetPosition.Latitude), Algorithms.Deg2Rad(e.TargetPosition.Longitude),
                                       Algorithms.WGS84Ellipsoid, Algorithms.VNC_DEF_EPSILON, Algorithms.VNC_DEF_IT_LIMIT,
                                       out dst_error, out fwd_az_rad, out rev_az_rad, out it_cnt);
                    Console.WriteLine(string.Format("Distance on the ellipsoid between actual target location and calculated: {0:F03} m", dst_error));
                };

            for (int i = 0; i < bases_number; i++)
            {
                Console.WriteLine(string.Format("Adding new base point: {0}", toaBases[i]));
                vlblCore.AddMeasurement(new VLBLTOAMeasurement(toaBases[i].Latitude, toaBases[i].Longitude, toaBases[i].Depth, toaBases[i].SlantRange));
            }

            #endregion

            Console.WriteLine();
            


            #region Misc

            Console.WriteLine("\r\nPress a key to start TOA_NLM3D_Solve() tests...");
            Console.ReadKey();

            List<TOABasePoint> base_points = new List<TOABasePoint>();
    
            // random relevant values for actual point position
            double actual_x = 20.0;
            double actual_y = -20.0;
            double actual_z = 20.0;

            double x_prev = 0.0;
            double y_prev = 0.0;
            double z_prev = 0.0;
            
            double x, y, z, d;

            x = -60.0;
            y = 20.0;
            z = 10.0;
            d =  Algorithms.Dist3D(x, y, z, actual_x, actual_y, actual_z);
            base_points.Add(new TOABasePoint(x, y, z, d));

            x = -20.0;
            y = 60.0;
            z = 15.0;
            d = Algorithms.Dist3D(x, y, z, actual_x, actual_y, actual_z);
            base_points.Add(new TOABasePoint(x, y, z, d));
        
            x = 20.0;
            y = 60.0;
            z = 20.0;
            d = Algorithms.Dist3D(x, y, z, actual_x, actual_y, actual_z);
            base_points.Add(new TOABasePoint(x, y, z, d));
        
            x = 40.0;
            y = 10.0;
            z = 25.0;
            d = Algorithms.Dist3D(x, y, z, actual_x, actual_y, actual_z);
            base_points.Add(new TOABasePoint(x, y, z, d));

            x = 50.0;
            y = -20.0;
            z = 30.0;
            d = Algorithms.Dist3D(x, y, z, actual_x, actual_y, actual_z);
            base_points.Add(new TOABasePoint(x, y, z, d));

            x = 30.0;
            y = -50.0;
            z = 35.0;
            d = Algorithms.Dist3D(x, y, z, actual_x, actual_y, actual_z);
            base_points.Add(new TOABasePoint(x, y, z, d));

            x = -20.0;
            y = -40.0;
            z = 40.0;
            d = Algorithms.Dist3D(x, y, z, actual_x, actual_y, actual_z);
            base_points.Add(new TOABasePoint(x, y, z, d));

            x = -50.0;
            y = -20.0;
            z = 45.0;
            d = Algorithms.Dist3D(x, y, z, actual_x, actual_y, actual_z);
            base_points.Add(new TOABasePoint(x, y, z, d));


            foreach (var base_point in base_points)
            {
                x_prev += base_point.X;
                y_prev += base_point.Y;
                z_prev += base_point.Z;
            }

            x_prev /= base_points.Count;
            y_prev /= base_points.Count;
            z_prev /= base_points.Count;

            double simplexSize = 10.0;
            double x_best = x_prev;
            double y_best = y_prev;
            double z_best = z_prev;
            

            Algorithms.TOA_NLM3D_Solve(base_points.ToArray(), x_prev, y_prev, z_prev, Algorithms.NLM_DEF_IT_LIMIT, Algorithms.NLM_DEF_PREC_THRLD, simplexSize,
                out x_best, out y_best, out z_best, out radialError, out it_cnt);

            Console.WriteLine();
            Console.WriteLine(string.Format("Actual location (x, y, z): ({0:F03}, {1:F03}, {2:F03})", actual_x, actual_y, actual_z));
            Console.WriteLine(string.Format("First approximation (x, y, z): ({0:F03}, {1:F03}, {2:F03})", x_prev, y_prev, z_prev));
            Console.WriteLine(string.Format("Estimated location (x, y, z): ({0:F03}, {1:F03}, {2:F03})", x_best, y_best, z_best));
            Console.WriteLine(string.Format("Radial error: {0:F03} m, iterations: {1}", radialError, it_cnt));

            Console.WriteLine();
            Console.WriteLine("\r\nPress a key to start TDOA_NLM3D_Solve() tests...");            
            Console.ReadKey();
            
            //                        
            List<TDOABaseline> base_lines = new List<TDOABaseline>();                                     

            for (int i = 0; i < base_points.Count - 1; i++)
            {
                for (int j = i + 1; j < base_points.Count; j++)
                {
                    base_lines.Add(new TDOABaseline(base_points[i].X, base_points[i].Y, base_points[i].Z, 
                                                    base_points[j].X, base_points[j].Y, base_points[j].Z, base_points[i].D - base_points[j].D));
                }
            }

            radialError = 0.0;
            it_cnt = 0;
            Algorithms.TDOA_NLM3D_Solve(base_lines.ToArray(), x_prev, y_prev, z_prev, Algorithms.NLM_DEF_IT_LIMIT, 
                Algorithms.NLM_DEF_PREC_THRLD, simplexSize, out x_best, out y_best, out z_best, out radialError, out it_cnt);

            Console.WriteLine();
            Console.WriteLine(string.Format("Actual location (x, y, z): ({0:F03}, {1:F03}, {2:F03})", actual_x, actual_y, actual_z));
            Console.WriteLine(string.Format("First approximation (x, y, z): ({0:F03}, {1:F03}, {2:F03})", x_prev, y_prev, z_prev));
            Console.WriteLine(string.Format("Estimated location (x, y, z): ({0:F03}, {1:F03}, {2:F03})", x_best, y_best, z_best));
            Console.WriteLine(string.Format("Radial error: {0:F03} m, iterations: {1}", radialError, it_cnt));

            Console.WriteLine();
            Console.WriteLine("Press a key to exit...");
            Console.ReadKey();

            #endregion

        }
    }
}
