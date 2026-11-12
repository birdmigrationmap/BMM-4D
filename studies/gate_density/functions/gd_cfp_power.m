function power = gd_cfp_power(raw_file, bin_edges_km, max_elevation_deg)
%GD_CFP_POWER Reported and removed echo power per cut and range bin.
%   POWER = GD_CFP_POWER(RAW_FILE, BIN_EDGES_KM, MAX_ELEVATION_DEG) reads REF
%   and CFP of every cut below MAX_ELEVATION_DEG (one per elevation: the
%   surveillance pass of a split cut, as GD_I and GD_L) and sums, per range
%   bin, the linear reflectivity the radar reported (P_OUT, no echo = 0) and
%   the same plus the power the clutter filter removed (P_ALL; CFP codes
%   >= 8 remove code - 8 dB). P_OUT / P_ALL is the share of the echo the
%   radar kept: 1 if the filter removed only clutter, the lowest possible
%   share of the birds kept if it removed only birds.
%   POWER has elev [n_cut 1], is_split [n_cut 1] (the elevation is scanned
%   twice), n_az [n_cut 1], and p_out, p_all, n_gate [n_bin n_cut].

cuts = bmm4d.radar.read_level2_moments(raw_file, ["REF" "CFP"], max_elevation_deg);
cuts = cuts(arrayfun(@(cut) ~isempty(cut.REF), cuts));
elevation_key = round([cuts.elev] * 2.5) / 2.5;
is_split = arrayfun(@(k) sum(elevation_key == elevation_key(k)) > 1, 1:numel(cuts));
[~, keep] = unique(elevation_key, "stable");
cuts = cuts(keep);
is_split = is_split(keep);
n_bin = numel(bin_edges_km) - 1;
n_cut = numel(cuts);
power = struct("elev", [cuts.elev]', "is_split", is_split(:), "n_az", arrayfun(@(cut) numel(cut.az), cuts)', ...
    "p_out", zeros(n_bin, n_cut), "p_all", zeros(n_bin, n_cut), "n_gate", zeros(n_bin, n_cut));
for k = 1:n_cut
    cut = cuts(k);
    if isempty(cut.CFP_code)
        continue
    end
    n_gate = min(size(cut.REF, 1), size(cut.CFP_code, 1));
    z = 10 .^ (cut.REF(1:n_gate, :) / 10);
    z(cut.REF_code(1:n_gate, :) == 0) = 0;
    z(cut.REF_code(1:n_gate, :) == 1) = nan;          % range folded: unknown
    code = double(cut.CFP_code(1:n_gate, :));
    removed_db = max(code - 8, 0);
    bin = discretize(double(cut.REF_range(1:n_gate)) / 1000, bin_edges_km);
    rows = ~isnan(bin);
    valid = ~isnan(z(rows, :));
    z_out = z(rows, :);
    z_out(~valid) = 0;
    z_all = z_out .* 10 .^ (removed_db(rows, :) / 10);
    power.p_out(:, k) = accumarray(bin(rows), sum(z_out, 2), [n_bin 1]);
    power.p_all(:, k) = accumarray(bin(rows), sum(z_all, 2), [n_bin 1]);
    power.n_gate(:, k) = accumarray(bin(rows), sum(valid, 2), [n_bin 1]);
end
end
