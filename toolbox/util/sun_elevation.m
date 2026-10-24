function elevation = sun_elevation(lat, lon, time_utc)
%SUN_ELEVATION Solar elevation angle (deg) at LAT/LON (arrays) and one UTC time.
%   NOAA approximation (declination, equation of time, hour angle), about
%   0.5 deg accuracy, without refraction. Enough to cut nights at twilight.

time_utc.TimeZone = "";
day_fraction = (day(time_utc, "dayofyear") - 1) + (hour(time_utc) - 12) / 24 ...
    + minute(time_utc) / 1440 + second(time_utc) / 86400;
gamma = 2 * pi / 365 * day_fraction;
equation_of_time = 229.18 * (0.000075 + 0.001868 * cos(gamma) - 0.032077 * sin(gamma) ...
    - 0.014615 * cos(2 * gamma) - 0.040849 * sin(2 * gamma));
declination = 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) ...
    - 0.006758 * cos(2 * gamma) + 0.000907 * sin(2 * gamma) ...
    - 0.002697 * cos(3 * gamma) + 0.00148 * sin(3 * gamma);

% True solar time (minutes) and hour angle (deg)
minutes_utc = hour(time_utc) * 60 + minute(time_utc) + second(time_utc) / 60;
solar_time = minutes_utc + equation_of_time + 4 * lon;
hour_angle = solar_time / 4 - 180;

elevation = asind(sind(lat) .* sin(declination) + cosd(lat) .* cos(declination) .* cosd(hour_angle));
end
