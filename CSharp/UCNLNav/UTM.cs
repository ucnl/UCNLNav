using System;

namespace UCNLNav
{
    public class UTM
    {
        public double Easting { get; set; }
        public double Northing { get; set; }
        public int Zone { get; set; }
        public char Hemisphere { get; set; }

        private static readonly Ellipsoid _wgs84 = Algorithms.WGS84Ellipsoid;
        private const double K0 = 0.9996; // Масштабный коэффициент UTM

        /// <summary>
        /// Конвертирует географические координаты (широта/долгота) в UTM (WGS84)
        /// </summary>
        /// <param name="latitude">Широта в градусах (-90..90)</param>
        /// <param name="longitude">Долгота в градусах (-180..180)</param>
        /// <returns>UTM координаты</returns>
        public static UTM FromLatLon(double latitude, double longitude)
        {
            // Проверка диапазонов
            if (latitude < -90 || latitude > 90)
                throw new ArgumentOutOfRangeException(nameof(latitude), "Широта должна быть в диапазоне -90..90");

            if (longitude < -180 || longitude > 180)
                throw new ArgumentOutOfRangeException(nameof(longitude), "Долгота должна быть в диапазоне -180..180");

            var utm = new UTM();

            // Вычисляем зону UTM (1-60)
            utm.Zone = CalculateZone(latitude, longitude);
            utm.Hemisphere = latitude >= 0 ? 'N' : 'S';

            // Конвертируем в радианы
            double latRad = Algorithms.Deg2Rad(latitude);
            double lonRad = Algorithms.Deg2Rad(longitude);

            // Центральный меридиан зоны
            double lon0 = (utm.Zone - 1) * 6 - 180 + 3;
            double lon0Rad = Algorithms.Deg2Rad(lon0);

            // Параметры эллипсоида
            double a = _wgs84.MajorSemiAxis_m;
            double e = _wgs84.Eccentricity;
            double eSq = _wgs84.EccentricitySq;

            // Вычисляем параметры для преобразования
            double nu = a / Math.Sqrt(1 - eSq * Math.Pow(Math.Sin(latRad), 2));
            double t = Math.Pow(Math.Tan(latRad), 2);
            double c = eSq * Math.Pow(Math.Cos(latRad), 2) / (1 - eSq);
            double l = (lonRad - lon0Rad) * Math.Cos(latRad);

            // Меридианная дуга
            double m = CalculateMeridionalArc(latRad, a, eSq);

            // Вычисляем Easting
            double l2 = l * l;
            double l4 = l2 * l2;
            double l6 = l4 * l2;

            double easting = K0 * nu * (l + (1 - t + c) * l2 * l / 6 +
                                       (5 - 18 * t + t * t + 72 * c - 58 * eSq) * l4 * l / 120) + 500000;

            // Вычисляем Northing
            double northing = K0 * (m + nu * Math.Tan(latRad) * (l2 / 2 +
                                   (5 - t + 9 * c + 4 * c * c) * l4 / 24 +
                                   (61 - 58 * t + t * t + 600 * c - 330 * eSq) * l6 / 720));

            // Для южного полушария добавляем 10,000,000 метров
            if (utm.Hemisphere == 'S')
                northing += 10000000;

            utm.Easting = Math.Round(easting, 2);
            utm.Northing = Math.Round(northing, 2);

            return utm;
        }

        /// <summary>
        /// Конвертирует UTM координаты обратно в географические (широта/долгота)
        /// </summary>
        /// <param name="easting">Координата X в метрах (0..1000000)</param>
        /// <param name="northing">Координата Y в метрах (0..10000000 для южного полушария)</param>
        /// <param name="zone">Номер зоны UTM (1..60)</param>
        /// <param name="hemisphere">Полушарие: 'N' или 'S'</param>
        /// <returns>Широта и долгота в градусах</returns>
        public static (double Latitude, double Longitude) ToLatLon(double easting, double northing, int zone, char hemisphere)
        {
            // Проверка параметров
            if (zone < 1 || zone > 60)
                throw new ArgumentOutOfRangeException(nameof(zone), "Зона UTM должна быть в диапазоне 1..60");

            if (hemisphere != 'N' && hemisphere != 'S')
                throw new ArgumentException("Полушарие должно быть 'N' или 'S'", nameof(hemisphere));

            // Корректировка для южного полушария
            if (hemisphere == 'S')
                northing -= 10000000;

            // Параметры эллипсоида WGS84
            double a = _wgs84.MajorSemiAxis_m;
            double e = _wgs84.Eccentricity;
            double eSq = _wgs84.EccentricitySq;

            // Центральный меридиан зоны
            double lon0 = (zone - 1) * 6 - 180 + 3;
            double lon0Rad = Algorithms.Deg2Rad(lon0);

            // Убираем масштабный коэффициент и смещение по X
            double x = easting - 500000;
            double y = northing;

            // Вычисляем меридианную дугу для начального приближения
            double m = y / K0;

            // Находим широту методом итераций (обратная меридианная дуга)
            double mu = m / (a * (1 - eSq / 4 - 3 * eSq * eSq / 64 - 5 * Math.Pow(eSq, 3) / 256));

            double e1 = (1 - Math.Sqrt(1 - eSq)) / (1 + Math.Sqrt(1 - eSq));

            double latRad = mu +
                           (3 * e1 / 2 - 27 * Math.Pow(e1, 3) / 32) * Math.Sin(2 * mu) +
                           (21 * e1 * e1 / 16 - 55 * Math.Pow(e1, 4) / 32) * Math.Sin(4 * mu) +
                           (151 * Math.Pow(e1, 3) / 96) * Math.Sin(6 * mu) +
                           (1097 * Math.Pow(e1, 4) / 512) * Math.Sin(8 * mu);

            // Итерационное уточнение
            for (int i = 0; i < 10; i++)
            {
                double sinLat_i = Math.Sin(latRad);
                double cosLat_i = Math.Cos(latRad);
                double tanLat_i = Math.Tan(latRad);

                double nu_i = a / Math.Sqrt(1 - eSq * sinLat_i * sinLat_i);
                double t_i = tanLat_i * tanLat_i;
                double c_i = eSq * cosLat_i * cosLat_i / (1 - eSq);

                double x2_i = x * x;
                double x3_i = x2_i * x;
                double x4_i = x3_i * x;
                double x5_i = x4_i * x;
                double x6_i = x5_i * x;

                double term1_i = 1;
                double term2_i = (1 + 2 * t_i + c_i) * x2_i / (2 * nu_i * nu_i);
                double term3_i = (5 + 28 * t_i + 24 * t_i * t_i + 6 * c_i + 8 * eSq) * x4_i / (24 * Math.Pow(nu_i, 4));
                double term4_i = (61 + 662 * t_i + 1320 * t_i * t_i + 720 * Math.Pow(t_i, 3)) * x6_i / (720 * Math.Pow(nu_i, 6));

                double latNew = latRad - (y / K0 - CalculateMeridionalArc(latRad, a, eSq)) / (a * (1 - eSq) / Math.Pow(1 - eSq * sinLat_i * sinLat_i, 1.5));

                if (Math.Abs(latNew - latRad) < 1e-12)
                {
                    latRad = latNew;
                    break;
                }
                latRad = latNew;
            }

            // Вычисляем долготу
            double tanLat = Math.Tan(latRad);
            double cosLat = Math.Cos(latRad);
            double nu = a / Math.Sqrt(1 - eSq * Math.Sin(latRad) * Math.Sin(latRad));
            double t = tanLat * tanLat;
            double c = eSq * cosLat * cosLat / (1 - eSq);

            double x2 = x * x;
            double x3 = x2 * x;
            double x4 = x3 * x;
            double x5 = x4 * x;

            double term1 = x / (nu * K0);
            double term2 = -(1 + 2 * t + c) * x3 / (6 * Math.Pow(nu, 3) * K0);
            double term3 = (5 + 28 * t + 24 * t * t + 6 * c + 8 * eSq) * x5 / (120 * Math.Pow(nu, 5) * K0);

            double lonDiff = term1 + term2 + term3;
            double lonRad = lon0Rad + lonDiff / cosLat;

            // Конвертируем в градусы
            double latitude = Algorithms.Rad2Deg(latRad);
            double longitude = Algorithms.Rad2Deg(lonRad);

            return (latitude, longitude);
        }

        /// <summary>
        /// Вычисляет зону UTM с учётом особых случаев
        /// </summary>
        private static int CalculateZone(double latitude, double longitude)
        {
            int zone = (int)Math.Floor((longitude + 180) / 6) + 1;

            // Особые случаи для Норвегии
            if (latitude >= 56 && latitude < 64 && longitude >= 3 && longitude < 12)
                zone = 32;

            // Особые случаи для Шпицбергена
            if (latitude >= 72 && latitude < 84 && longitude >= 0 && longitude < 42)
            {
                if (longitude < 9) zone = 31;
                else if (longitude < 21) zone = 33;
                else if (longitude < 33) zone = 35;
                else zone = 37;
            }

            return zone;
        }

        /// <summary>
        /// Вычисляет меридианную дугу от экватора до заданной широты
        /// </summary>
        private static double CalculateMeridionalArc(double latRad, double a, double eSq)
        {
            double e4 = eSq * eSq;
            double e6 = e4 * eSq;
            double e8 = e6 * eSq;

            double a0 = 1 - eSq / 4 - 3 * e4 / 64 - 5 * e6 / 256 - 175 * e8 / 16384;
            double a2 = 3 * eSq / 8 + 3 * e4 / 32 + 45 * e6 / 1024 + 105 * e8 / 4096;
            double a4 = 15 * e4 / 256 + 45 * e6 / 1024 + 525 * e8 / 16384;
            double a6 = 35 * e6 / 3072 + 175 * e8 / 12288;
            double a8 = 315 * e8 / 131072;

            // Вычисление меридианной дуги
            double m = a * (a0 * latRad -
                           a2 * Math.Sin(2 * latRad) / 2 +
                           a4 * Math.Sin(4 * latRad) / 4 -
                           a6 * Math.Sin(6 * latRad) / 6 +
                           a8 * Math.Sin(8 * latRad) / 8);

            return m;
        }

        public override string ToString()
        {
            return $"{Zone}{Hemisphere} E:{Easting:F2} N:{Northing:F2}";
        }
    }
}