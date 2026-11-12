function names = ar_holdout_radars()
%AR_HOLDOUT_RADARS Radars left out together in the night-set prior tests (dev).
%   NAMES = AR_HOLDOUT_RADARS() returns 8 us_east radars chosen by a fixed
%   rule: inland (onwater50 < 1%) radars with a neighbour in each of the
%   four quadrants within twice the median nearest-radar distance
%   (171 km), then picked greedily to be as far as possible from each
%   other, starting from the one nearest the centroid. They are 727-875 km
%   apart and each keeps a radar neighbour 117-221 km away, so leaving them
%   out together tests interpolation between radars, not extrapolation.

names = ["KPAH" "KBIS" "KENX" "KLBB" "KVAX" "KARX" "KPOE" "KFCX"];
end
