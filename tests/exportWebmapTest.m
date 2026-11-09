classdef exportWebmapTest < matlab.unittest.TestCase
    %EXPORTWEBMAPTEST Export of a field for the web map viewer (bmm4d.export_webmap).
    % Synthetic case: 11 x 11 columns of 3 layers (0.1 deg, 200 m) over flat
    % ground at sea level, one radar with a uniform wind (5, -3) m/s, two
    % frames of uniform density 10 and 20 birds km^-3 (columns of 6 and
    % 12 birds km^-2). The export is decoded as the viewer does.

    properties
        ProjectRoot
        OutDir
        Meta
    end

    methods (TestClassSetup)
        function exportSyntheticCase(testCase)
            testCase.ProjectRoot = fileparts(fileparts(mfilename("fullpath")));
            if isempty(mpmlist("bmm4d"))
                mpminstall(fullfile(testCase.ProjectRoot, "toolbox"), InPlace=true, Temporary=true, Prompt=false, Verbosity="quiet");
            end
            % Grid, one radar and the case settings the export reads
            g.lon = -80:0.1:-79;
            g.lat = 35:0.1:36;
            g.x = (g.lon + 80) * 91000;
            g.y = (g.lat - 35) * 111000;
            g.z = reshape([100 300 500], 1, 1, 3);
            g.dz = 200;
            g.sz = [numel(g.lon) numel(g.lat) numel(g.z)];
            g.dem = zeros(g.sz(1:2));
            g.f_X = @(lon) (lon + 80) * 91000;
            g.f_Y = @(lat) (lat - 35) * 111000;
            scan = struct("lon", -79.5, "lat", 35.5, "height", 0, "uv", [5 -3]);
            radar_scans = {scan scan};
            radar_selected = table("KTST", -79.5, 35.5, VariableNames=["name" "lon" "lat"]);
            c = struct("name", "webmap_test", "times", datetime(2026, 10, 22, 2, [0 30], 0), ...
                "sun_limit_deg", -6, "inversion", struct("velocity_interpolation", "gaussian", "velocity_length_scale_m", 150000, ...
                "velocity_idw_neighbors", 4, "velocity_idw_power", 2));
            field = [10 20] .* ones(prod(g.sz), 1);
            testCase.OutDir = fullfile(tempname, "webmap", "test_run");
            bmm4d.export_webmap(c, field, active_id=(1:prod(g.sz))', run="test_run", note="unit test", ...
                case_data=struct("g", g, "radar_scans", {radar_scans}, "radar_selected", radar_selected), out_dir=testCase.OutDir);
            testCase.Meta = jsondecode(fileread(fullfile(testCase.OutDir, "meta.json")));
        end
    end

    methods (TestClassTeardown)
        function removeExport(testCase)
            if isfolder(fileparts(fileparts(testCase.OutDir)))
                rmdir(fileparts(fileparts(testCase.OutDir)), "s");
            end
        end
    end

    methods (Test)
        function testFilesAndRun(testCase)
            % One column image per frame, one wind image per interval, the run described
            testCase.verifyEqual(numel(dir(fullfile(testCase.OutDir, "col_*.png"))), 2)
            testCase.verifyEqual(numel(dir(fullfile(testCase.OutDir, "wind_*.png"))), 1)
            meta = testCase.Meta;
            testCase.verifyEqual(string(meta.format), "bmm4d-webmap-1")
            testCase.verifyEqual([string(meta.case) string(meta.run) string(meta.note)], ["webmap_test" "test_run" "unit test"])
            testCase.verifyEqual(string(meta.times), ["2026-10-22T02:00:00Z"; "2026-10-22T02:30:00Z"])
            % One radar is still a JSON list (the viewer maps over it)
            raw = fileread(fullfile(testCase.OutDir, "meta.json"));
            testCase.verifyTrue(contains(raw, regexpPattern('"radars":\s*\[')))
        end

        function testColumnDecodes(testCase)
            % The viewer's decoding gives the column density within the 2-3% encoding steps
            meta = testCase.Meta;
            decode = @(code) meta.column_floor * (exp((double(code) - 1) / 254 * log1p(meta.column_max / meta.column_floor)) - 1);
            for i_t = 1:2
                image = imread(fullfile(testCase.OutDir, sprintf("col_%02d.png", i_t - 1)));
                inside = image(:, :, 2) == 255;
                testCase.verifyGreaterThan(nnz(inside), 0.8 * numel(inside))
                red = image(:, :, 1);
                testCase.verifyEqual(decode(red(inside)), repmat(6 * i_t, nnz(inside), 1), RelTol=0.03)
            end
        end

        function testWindDecodes(testCase)
            % The birds' wind is the radar's (5, -3) m/s everywhere, within the 0.4 m/s encoding step
            meta = testCase.Meta;
            image = double(imread(fullfile(testCase.OutDir, "wind_00.png")));
            decode = @(code) code / 255 * 2 * meta.wind_limit - meta.wind_limit;
            u = decode(image(:, :, 1));
            v = decode(image(:, :, 2));
            testCase.verifyEqual(u(:), repmat(5, numel(u), 1), AbsTol=0.4)
            testCase.verifyEqual(v(:), repmat(-3, numel(v), 1), AbsTol=0.4)
        end

        function testMapPlacement(testCase)
            % Longitude: the grid cells' edges. Latitude: rows of the longitude
            % step in Mercator y from the northernmost cell centre, so each
            % edge is within one row (< 0.1 deg) outside the outer centres
            meta = testCase.Meta;
            testCase.verifyEqual([meta.lon_west meta.lon_east], [-80.05 -78.95], AbsTol=1e-9)
            testCase.verifyEqual(meta.width, 11)
            testCase.verifyTrue(meta.lat_north > 36 && meta.lat_north < 36.1)
            testCase.verifyTrue(meta.lat_south > 34.9 && meta.lat_south < 35)
        end
    end
end
