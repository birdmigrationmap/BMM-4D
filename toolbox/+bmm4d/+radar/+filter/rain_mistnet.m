function scan = rain_mistnet(scan, parms)
%RAIN_MISTNET Remove rain classified by MistNet.
%   With mistnetScreen, MistNet runs only on volumes with rain-like echo
%   (has_rain_like_echo); others are kept whole. With removeRain
%   "mistnet_rho", gates whose median-smoothed rhoHV exceeds
%   rainRhoThreshold are removed too (rain edges, drizzle, unclassified
%   cells that MistNet misses). Sets scan.mistnet_used.

if parms.mistnetScreen && ~has_rain_like_echo(scan)
    % Screened as rain-free: MistNet would remove almost nothing
    scan.mistnet_used = false;
    return
end
scan.mistnet_used = true;
[preds, classes] = bmm4d.radar.mistnet_polar_radar(scan);
for i_elev = 1:numel(scan.gates)
    mask = preds{i_elev} == classes.rain;
    if parms.removeRain == "mistnet_rho"
        % The median over ~1 km x 3 radials keeps isolated high rhoHV gates
        % of biology.
        [~, i_elev_rh] = min(abs(scan.gates(i_elev).elev - [scan.rh.sweeps.elev]));
        rh_data = sweep2mat(scan.rh.sweeps(i_elev_rh));
        if size(rh_data, 2) == size(mask, 2)
            rh_data = smoothdata(smoothdata(rh_data, 1, "movmedian", 5, "omitnan"), 2, "movmedian", 3, "omitnan");
            min_nbins = min(size(rh_data, 1), size(mask, 1));
            mask(1:min_nbins, :) = mask(1:min_nbins, :) | rh_data(1:min_nbins, :) > parms.rainRhoThreshold;
        end
    end
    scan.gates(i_elev) = bmm4d.radar.filter.remove_rain(scan.gates(i_elev), mask);
end
end

function rainy = has_rain_like_echo(radar_scan)
%HAS_RAIN_LIKE_ECHO Cheap screen for precipitation before running MistNet.
%   Mean area (km2) per sweep below 2 deg, within 3-150 km, of gates with
%   dBZ >= 10 and rhoHV >= 0.97 (widespread weak rain) or dBZ >= 20 and
%   rhoHV >= 0.95 (rain cells). Thresholds from 287 random migration
%   volumes (studies/cleaning_study): volumes below both limits (42%) lose at
%   most 2.6 birds km^-2 of rain that MistNet would have removed.
area_weak = 0;
area_cell = 0;
low_sweeps = find([radar_scan.dz.sweeps.elev] < 2);
for i_elev = low_sweeps
    [dz, sweep_range, sweep_az] = sweep2mat(radar_scan.dz.sweeps(i_elev));
    [~, i_elev_rh] = min(abs(radar_scan.dz.sweeps(i_elev).elev - [radar_scan.rh.sweeps.elev]));
    rh = sweep2mat(radar_scan.rh.sweeps(i_elev_rh));
    if size(rh, 2) ~= size(dz, 2)
        rainy = true;   % cannot screen: run MistNet
        return
    end
    n = min(size(dz, 1), size(rh, 1));
    dz = dz(1:n, :);
    rh = rh(1:n, :);
    sweep_range = sweep_range(1:n);
    gate_area = sweep_range * deg2rad(median(diff(sweep_az))) * median(diff(sweep_range)) / 1e6;
    near = sweep_range >= 3000 & sweep_range <= 150000;
    area_weak = area_weak + sum(gate_area .* near .* sum(dz >= 10 & rh >= 0.97, 2));
    area_cell = area_cell + sum(gate_area .* near .* sum(dz >= 20 & rh >= 0.95, 2));
end
n_sweep = max(numel(low_sweeps), 1);
rainy = area_weak / n_sweep >= 200 || area_cell / n_sweep >= 30;
end
