function cut_code = nr_cut_types(observations, radar_selected, manifest)
%NR_CUT_TYPES Waveform of every observation entry from the Level II headers (dev).
%   CUT_CODE = NR_CUT_TYPES(OBSERVATIONS, RADAR_SELECTED, MANIFEST) is
%   uint8 [n_row n_time]: the waveform of the cut each observation row
%   comes from at each time, read from the raw file's scan strategy
%   (nr_read_vcp): 1 split (contiguous surveillance reflectivity), 2-3
%   Doppler, 4 batch, 5 staggered pulse pair, 0 unknown (no raw file, or
%   no cut within 0.1 deg of the row's elevation). The VCP can change
%   during a night, so the code is per time. MANIFEST is the case's
%   manifest.mat (radars, raw_files [n_radar n_time]); RADAR_SELECTED the
%   case's radar table (row_radar_id indexes it).

n_time = size(observations.values, 2);
cut_code = zeros(size(observations.values), "uint8");
for i_r = 1:height(radar_selected)
    i_manifest = find(manifest.radars == radar_selected.name(i_r), 1);
    rows = find(observations.row_radar_id == i_r);
    if isempty(i_manifest) || isempty(rows), continue, end
    [row_elev, ~, elev_id] = unique(round(double(observations.row_elev_deg(rows)), 2));
    for i_t = 1:n_time
        raw_file = manifest.raw_files(i_manifest, i_t);
        if raw_file == "" || ~isfile(raw_file), continue, end
        cuts = nr_read_vcp(raw_file).cuts;
        if isempty(cuts), continue, end
        code = zeros(numel(row_elev), 1, "uint8");
        for k = 1:numel(row_elev)
            near = abs(cuts.elev - row_elev(k)) < 0.1;
            if ~any(near), continue, end
            if any(cuts.split(near))
                code(k) = 1;
            else
                code(k) = cuts.waveform(find(near, 1));
            end
        end
        cut_code(rows, i_t) = code(elev_id);
    end
end
end
