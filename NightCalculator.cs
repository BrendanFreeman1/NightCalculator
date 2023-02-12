public class NightGameProcess

    {
        /// <summary>
        /// Checks and sets the IsNight flag on each event passed in
        /// </summary>
        /// <param name="events">The events you want to check</param>
        /// <returns></returns>
        public void Execute(IEnumerable<Event> events)
        {
            foreach(var @event in events) 
            {
                //if (IsEventNightGame(@event))
                //    //@event.IsNight = true;
                //else
                //    //@event.IsNight = false;
            }
        }

        /// <summary>
        /// Calculates if the passed event took place at night or not
        /// </summary>
        /// <param name="event">The Event you wish to check</param>
        /// <returns>A boolean representing if the event took place at night or not</returns>
        public bool IsEventNightGame(Event @event)
         {
            var eventCoordinates = @event.Location.Coordinates;
            var eventDate = @event.Date.Date;

            var nightTime = GetNightTimeInUTC(eventCoordinates, eventDate);

            if(@event.Date.TimeOfDay > nightTime)
                return true;

            return false;
        }

        /// <summary>
        /// Calculates the time night starts at the passed location on the passed date. I'm defining night as the sun being 12 or more degrees below the horizon
        /// </summary>
        /// <param name="coordinates">Latitude: North is a Positive number, South is a Negative number. Longitude: East is a Positive number, West is a Negative number</param>
        /// <param name="eventDate">The date the event took place on. Should be UTC</param>
        /// <returns>A TimeSpan representing the time night starts on that day</returns>
        public TimeSpan GetNightTimeInUTC(CoOrdinate coordinates, DateTime eventDate)
        {
            // Constants represent the suns angle below the horizon.
            const int civilTwilight = 0;
            const int nauticalTwilight = -6;
            const int astronomicalTwilight = -12;
            const int night = -18;

            var sunSetTime = GetSunSetTimeInUTC(coordinates, eventDate);
            var sunAngle = GetSunPosition(coordinates, sunSetTime);

            var timeAdjust = new TimeSpan(0, 1, 0);
            var testTime = sunSetTime;

            // Loop until the suns angle is less than 12degrees
            while(sunAngle >= astronomicalTwilight)
            {
                // Increase time by 1 minute
                testTime += timeAdjust;

                // Get new sun angle
                sunAngle = GetSunPosition(coordinates, testTime);
            }

            var nightTime = testTime.TimeOfDay;

            return nightTime;            
        }

        /* Code Converted to C# from JavaScript from https://gist.github.com/ruiokada/b28076d4911820ddcbbc */
        /// <summary>
        /// Calculates the Sunset time at the passed location on the passed date(UTC)
        /// </summary>
        /// <param name="coordinates">Latitude: North is a Positive number, South is a Negative number. Longitude: East is a Positive number, West is a Negative number</param>
        /// <param name="date">Date should be converted to UTC before passing</param>
        /// <returns>A DateTime approximating the sunset on the passed date</returns>
        public DateTime GetSunSetTimeInUTC(CoOrdinate coordinates, DateTime date)
        {
            var latitude = coordinates.Latitude;
            var longitude = coordinates.Longitude;

            double UTC_Rise;
            double UTC_Set;

            var radians = Math.PI / 180.0;
            var degrees = 180.0 / Math.PI;

            var a = Math.Floor((14 - (date.Month + 1.0)) / 12);
            var y = date.Year + 4800 - a;
            var m = (date.Month) + 12 * a - 3;

            var j_day = date.Date.Day + Math.Floor((153 * m + 2) / 5) + 365 * y + Math.Floor(y / 4) - Math.Floor(y / 100) + Math.Floor(y / 400) - 32045;
            var n_star = j_day - 2451545.0009 - longitude / 360.0;
            var n = Math.Floor(n_star + 0.5);
            var solar_noon = 2451545.0009 - longitude / 360.0 + n;
            var M = 356.0470 + 0.9856002585 * n;
            var C = 1.9148 * Math.Sin(M * radians) + 0.02 * Math.Sin(2 * M * radians) + 0.0003 * Math.Sin(3 * M * radians);
            var L = (M + 102.9372 + C + 180) % 360;
            var j_transit = solar_noon + 0.0053 * Math.Sin(M * radians) - 0.0069 * Math.Sin(2 * L * radians);
            var D = Math.Asin(Math.Sin(L * radians) * Math.Sin(23.45 * radians)) * degrees;
            var cos_omega = (Math.Sin(-0.83 * radians) - Math.Sin(latitude * radians) * Math.Sin(D * radians)) / (Math.Cos(latitude * radians) * Math.Cos(D * radians));

            // Get Julian dates of sunrise/sunset
            var omega = Math.Acos(cos_omega) * degrees;
            var j_set = j_transit + omega / 360.0;
            var j_rise = j_transit - omega / 360.0;

            // Get sunrise and sunset times in UTC
            var utcTimeRise = 24 * (j_rise - j_day) + 12;
            var utcTimeSet = 24 * (j_set - j_day) + 12;

            UTC_Rise = utcTimeRise % 24;
            UTC_Set = utcTimeSet % 24;

            //Turn result into a DateTime of the sunset
            var hour = Math.Truncate(UTC_Set);
            var minute = Math.Truncate(60 * (UTC_Set - Math.Truncate(UTC_Set)));

            var UTC_SunSet = new DateTime(date.Year, date.Month, date.Day, (int)hour, (int)minute, 00);

            return UTC_SunSet;
        }

        #region GetSunAltitude
        /* Code from https://physics.stackexchange.com/a/696035 */
        /// <summary>
        /// Get The Altitude of the sun at the passed location at the passed Date/Time
        /// </summary>
        /// <param name="coordinates">Latitude: North is a Positive number, South is a Negative number. Longitude: East is a Positive number, West is a Negative number</param>        /// <param name="date">Date should be converted to UTC before passing</param>
        /// <returns>A double representing the angle of the sun. 0 is Sunset, Negative is below the horizon, Positive is above the horizon</returns>
        private double GetSunPosition(CoOrdinate coordinates, DateTime date)
        {
            var latitude = coordinates.Latitude;
            var longitude = coordinates.Longitude;

            var dayNumber = GetDayNumber(date);
            var argumentOfPerihelion = Sun_ArgumentOfPerihelion(dayNumber);
            var eclipticObliquity = EclipticObliquity(dayNumber);

            // First, compute the eccentric anomaly E from the mean anomaly M and from the eccentricity e(degrees):
            //
            // E = M + e * (180 / pi) * sin(M) * (1.0 + e * cos(M))
            //
            // or(if E and M are expressed in radians):
            //
            // E = M + e * sin(M) * (1.0 + e * cos(M))

            var meanAnomaly = Sun_MeanAnomaly(dayNumber);
            var eccentricity = Sun_Eccentricity(dayNumber);
            var eccentricAnomaly = meanAnomaly + (180 / Math.PI) * eccentricity * Math.Sin(Deg2Rad(meanAnomaly)) * (1.0 + eccentricity * Math.Cos(Deg2Rad(meanAnomaly)));

            // Then compute the Sun's distance r and its true anomaly v from:

            // xv = r * cos(v) = cos(E) - e
            // yv = r * sin(v) = sqrt(1.0 - e * e) * sin(E)
            // v = atan2(yv, xv)
            // r = sqrt(xv * xv + yv * yv)

            // (note that the r computed here is later used as rs)

            var xv = Math.Cos(Deg2Rad(eccentricAnomaly)) - eccentricity;
            var yv = Math.Sqrt(1.0 - eccentricity * eccentricity) * Math.Sin(Deg2Rad(eccentricAnomaly));
            var v = Rad2Deg(Math.Atan2(yv, xv));
            var r = Math.Sqrt(xv * xv + yv * yv);

            // Now, compute the Sun's true longitude:

            // lonsun = v + w

            var sunTrueLongitude = Rev(v + argumentOfPerihelion);

            // Convert lonsun, r to ecliptic rectangular geocentric coordinates xs,ys:

            // xs = r * cos(lonsun)
            // ys = r * sin(lonsun)

            var xs = r * Math.Cos(Deg2Rad(sunTrueLongitude));
            var ys = r * Math.Sin(Deg2Rad(sunTrueLongitude));

            // Since the Sun always is in the ecliptic plane, zs is of course zero.
            // xs,ys is the Sun's position in a coordinate system in the plane of the ecliptic.
            // To convert this to equatorial, rectangular, geocentric coordinates, compute:

            // xe = xs
            // ye = ys * cos(ecl)
            // ze = ys * sin(ecl)

            var xe = xs;
            var ye = ys * Math.Cos(Deg2Rad(eclipticObliquity));
            var ze = ys * Math.Sin(Deg2Rad(eclipticObliquity));

            // Finally, compute the Sun's Right Ascension (RA) and Declination (Dec):
            // RA = atan2(ye, xe)
            // Dec = atan2(ze, sqrt(xe * xe + ye * ye))

            var rightAscension = Rad2Deg(Math.Atan2(ye, xe));
            var declination = Rad2Deg(Math.Atan2(ze, Math.Sqrt(xe * xe + ye * ye)));

            // Calculate Greenwich Sidereal Time, Sidereal Time and the Sun's Hour Angle

            var sunMeanLongitude = Sun_MeanLongitude(dayNumber);
            var gmst0 = sunMeanLongitude / 15 + 12;

            var siderealTime = RevTime(gmst0 + date.Hour + (date.Minute / 60F) + longitude / 15);
            var hourAngle = RevTime(siderealTime - rightAscension / 15);

            // Convert the Sun's Hour Angle and Declination to a rectangular coordinate system where the X
            // axis points to the celestial equator in the south, the Y axis to the horizon in the west,
            // and the Z axis to the north celestial pole.

            var x = Math.Cos(Deg2Rad(hourAngle * 15)) * Math.Cos(Deg2Rad(declination));
            var z = Math.Sin(Deg2Rad(declination));

            // Rotate this x,y,z axis system along an axis going east-west(Y axis) in such a way that the
            // Z axis will point to the zenith. At the North Pole, the angle of rotation will be zero since
            // there the north celestial pole already is in the zenith. At other latitudes the angle of
            // rotation becomes 90 - latitude.

            var zhor = x * Math.Cos(Deg2Rad(latitude)) + z * Math.Sin(Deg2Rad(latitude));

            // Compute altitude.
            var altitude = Rad2Deg(Math.Asin(zhor));

            return altitude;
        }
        private static double GetDayNumber(DateTime dt)
        {
            //http://www.stjarnhimlen.se/comp/ppcomp.html#5
            // The time scale in these formulae are counted in days.Hours, minutes, seconds are expressed as fractions of a day.

            // Day 0.0 occurs at 2000 Jan 0.0 UT(or 1999 Dec 31, 0:00 UT).

            // This "day number" d is computed as follows (y = year, m = month, D = date, UT = UT in hours + decimals):

            // d = 367 * y - 7 * (y + (m + 9) / 12) / 4 + 275 * m / 9 + D - 730530

            // Note that ALL divisions here should be INTEGER divisions.
            // Finally, include the time of the day, by adding:

            // d = d + UT / 24.0(this is a floating - point division)

            var d = 367 * dt.Year - 7 * (dt.Year + (dt.Month + 9) / 12) / 4 + 275 * dt.Month / 9 + dt.Day - 730530;
            double hm = dt.Hour + (dt.Minute / 60F);
            return d + hm / 24;
        }

        /// <summary>
        /// Longitude of Perihelion (w1)
        /// </summary>
        private static double Sun_LongitudeOfPerihelion(double dayNumber)
        {
            var n = Sun_LongitudeOfAscendingNode();
            var w = Sun_ArgumentOfPerihelion(dayNumber);
            return Rev(n + w);
        }

        /// <summary>
        /// Longitude of the Ascending Node (N)
        /// </summary>
        private static double Sun_LongitudeOfAscendingNode()
        {
            return 0.0D;
        }

        /// <summary>
        /// Argument of Perihelion (w)
        /// </summary>
        private static double Sun_ArgumentOfPerihelion(double dayNumber)
        {
            return 282.9404 + 4.70935e-5 * dayNumber;
        }

        /// <summary>
        /// Eccentricity (e) where 0 = circle, 0-1 = ellipse, and 1 = parabola
        /// </summary>
        private static double Sun_Eccentricity(double dayNumber)
        {
            return 0.016709 - 1.151e-9 * dayNumber;
        }

        /// <summary>
        /// Mean anomaly (M) (0 at perihelion; increases uniformly with time)
        /// </summary>
        private static double Sun_MeanAnomaly(double dayNumber)
        {
            return Rev(356.0470 + 0.9856002585 * dayNumber);
        }

        /// <summary>
        /// Ecliptic Obliquity (ecl)
        /// </summary>
        private static double EclipticObliquity(double dayNumber)
        {
            return 23.4393 - 3.563e-7 * dayNumber;
        }

        /// <summary>
        /// Mean Longitude (L)
        /// </summary>
        private static double Sun_MeanLongitude(double dayNumber)
        {
            var m = Sun_MeanAnomaly(dayNumber);
            var w1 = Sun_LongitudeOfPerihelion(dayNumber);
            return Rev(m + w1);
        }

        /// <summary>
        /// Convert degrees to radians
        /// </summary>
        private static double Deg2Rad(double angleDegrees)
        {
            return angleDegrees * (Math.PI / 180.0);
        }

        /// <summary>
        /// Convert radians to degrees
        /// </summary>
        private static double Rad2Deg(double angleRadians)
        {
            return angleRadians * (180.0 / Math.PI);
        }

        /// <summary>
        /// Revolution function, normalizes an angle to between 0 and 360 degrees by adding or subtracting even multiples of 360.
        /// </summary>
        private static double Rev(double x)
        {
            return x - Math.Floor(x / 360.0) * 360.0;
        }

        /// <summary>
        /// Revolution function, normalizes a time value (in hours) to between 0 and 24 by adding or subtracting even multiples of 24.
        /// </summary>
        private static double RevTime(double x)
        {
            return x - Math.Floor(x / 24.0) * 24.0;
        }
        #endregion
    }
