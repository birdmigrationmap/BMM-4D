function [sunset, sunrise] = getSunriseSunset(lat, lon, target_time)
% Return the sunset/sunrise bounds for the night nearest the target time.

day_start = dateshift(target_time, "start", "day");
sunrise_today = solar_event_utc(lat, lon, day_start, true);
sunset_today = solar_event_utc(lat, lon, day_start, false);

if target_time <= sunrise_today
    sunset = solar_event_utc(lat, lon, day_start - days(1), false);
    sunrise = sunrise_today;
elseif target_time >= sunset_today
    sunset = sunset_today;
    sunrise = solar_event_utc(lat, lon, day_start + days(1), true);
elseif target_time - sunrise_today <= sunset_today - target_time
    sunset = solar_event_utc(lat, lon, day_start - days(1), false);
    sunrise = sunrise_today;
else
    sunset = sunset_today;
    sunrise = solar_event_utc(lat, lon, day_start + days(1), true);
end

end

function event_time = solar_event_utc(lat, lon, day_start, is_sunrise)
% NOAA solar calculation in UTC for the requested date and location.

n = day(day_start, "dayofyear");
lng_hour = lon / 15;

if is_sunrise
    t = n + (6 - lng_hour) / 24;
else
    t = n + (18 - lng_hour) / 24;
end

M = 0.9856 * t - 3.289;
L = mod(M + 1.916 * sind(M) + 0.020 * sind(2 * M) + 282.634, 360);

RA = mod(rad2deg(atan(0.91764 * tand(L))), 360);
L_quadrant = floor(L / 90) * 90;
RA_quadrant = floor(RA / 90) * 90;
RA = (RA + L_quadrant - RA_quadrant) / 15;

sin_dec = 0.39782 * sind(L);
cos_dec = cosd(asind(sin_dec));
cosH = (cosd(90.833) - sin_dec * sind(lat)) / (cos_dec * cosd(lat));

if is_sunrise
    H = 360 - acosd(cosH);
else
    H = acosd(cosH);
end

T = H / 15 + RA - 0.06571 * t - 6.622;
UT = mod(T - lng_hour, 24);

event_time = day_start + hours(UT);
event_time.TimeZone = target_timezone(day_start);
end

function tz = target_timezone(t)
if isdatetime(t)
    tz = t.TimeZone;
else
    tz = "";
end
end
