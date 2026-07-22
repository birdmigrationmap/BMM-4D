function compute_F_elev(x, y, z, nxx, nyy, nzz, elevation_kernel, thr_cumsum, mode)
%COMPUTE_F_ELEV Precompute one elevation-specific forward kernel.
%   COMPUTE_F_ELEV(X, Y, Z, NXX, NYY, NZZ, ELEVATION_KERNEL, THR_CUMSUM, MODE)
%   approximates the radar beam footprint on a fine g and aggregates it
%   back to the coarse reconstruction g. MODE can be "fast" (default)
%   or "memory_efficient".
if nargin < 9 || isempty(mode)
    mode = "fast";
else
    mode = string(mode);
end
assert(any(mode == ["fast" "memory_efficient"]), "mode must be 'fast' or 'memory_efficient'.")

disp(elevation_kernel.name)
disp("mode: " + mode)

%% Define the g
assert(isscalar(unique(diff(x))))
dx = diff(x(1:2));
assert(isscalar(unique(diff(y))))
dy = diff(y(1:2));
assert(isscalar(unique(diff(z))))
dz = diff(z(1:2));

% Size of the coarse g
sz_coarse_grid = [numel(x) numel(y) numel(z)];

% Fine g
dxx=dx/nxx; dyy=dy/nyy; dzz=dz/nzz;
xx=(x(1)-dx/2+dxx/2):dxx:(x(end)+dx/2-dxx/2);
xx=single(reshape(xx,1,[]));
yy=(y(1)-dy/2+dyy/2):dyy:(y(end)+dy/2-dyy/2);
yy=single(reshape(yy,[],1));
zz=(z(1)-dz/2+dzz/2):dzz:(z(end)+dz/2-dzz/2);
zz=single(reshape(zz,1,1,[]));

sz_fine_grid = [numel(xx) numel(yy) numel(zz)];
sz_fine_xy = sz_fine_grid(1:2);
fprintf("g: %d x %d x %d coarse | %d x %d x %d fine\n", numel(x), numel(y), numel(z), numel(xx), numel(yy), numel(zz))


%% Convert the Cartesian g to polar coordinates accounting for Earth curvature
% xx,yy -> sigma, phi

% range along the ground (great circle distance)
SS = sqrt(xx.^2 + yy.^2);

% set cst
earth_radius = single(6371200); % from NARR GRIB file
multiplier = single(4/3); % correct for atmospheric refraction
r_e = earth_radius * multiplier; % earth effective radius

sigma = SS / r_e; clear SS
PP = single(pol2cmp(angle(xx + 1i * yy)));


%% Compute weighting matrix

% compute the std dev of the beam along phi and theta
hpbw = single(0.925); % half power beamwiths
simga_tp = hpbw/2/sqrt(2*log(2));

theta0 = elevation_kernel.elev;
simga_r = 0.35*median(diff(elevation_kernel.range));

if mode == "fast"
    [i_range_max, RR, TT] = build_full_polar_grid(sigma, zz, r_e, theta0, simga_tp, elevation_kernel.range);
    Felevf = cell(i_range_max, numel(elevation_kernel.az));
    Felevid = cell(i_range_max, numel(elevation_kernel.az));
    fprintf("fast mode: building full polar g, %d ranges x %d azimuths\n", i_range_max, numel(elevation_kernel.az))

    id_theta = find(TT>(theta0-3*simga_tp) & TT<(theta0+3*simga_tp));
    PP = repmat(PP, 1, 1, numel(zz));

    % Loop through all the ranges
    for i_range = 1:i_range_max
        r0 = elevation_kernel.range(i_range);
        id_range = id_theta(abs(RR(id_theta)-r0)<=3*simga_r);

        for i_azim = 1:numel(elevation_kernel.az)
            phi0 = elevation_kernel.az(i_azim);
            diff_phi = mod(PP(id_range)-phi0,360);
            diff_phi = min(360-diff_phi, diff_phi);
            id_range_phi = diff_phi<3*simga_tp;
            id_phi = id_range(id_range_phi);

            if numel(id_phi)==0
                continue
            end

            f = exp(-(RR(id_phi)-r0).^2/2/simga_r^2 -(TT(id_phi)-theta0).^2/simga_tp^2 - (diff_phi(id_range_phi)).^2/simga_tp^2 );
            [id_xx,id_yy,id_zz] = ind2sub(sz_fine_grid,id_phi);
            id_x = floor((id_xx-1)/(dx/dxx))+1;
            id_y = floor((id_yy-1)/(dy/dyy))+1;
            id_z = floor((id_zz-1)/(dz/dzz))+1;
            id_coarse = sub2ind(sz_coarse_grid,id_x,id_y,id_z);

            [f_keep, id_keep] = collapse_beam_weights(id_coarse, f, thr_cumsum);
            if isempty(f_keep)
                continue
            end

            Felevf{i_range,i_azim} = f_keep;
            Felevid{i_range,i_azim} = uint32([ ...
                i_range*ones(numel(id_keep),1), ...
                i_azim*ones(numel(id_keep),1), ...
                id_keep(:)]);

            if should_report_azimuth_progress(i_azim, numel(elevation_kernel.az))
                fprintf("  range %d/%d | azimuth %d/%d\n", i_range, i_range_max, i_azim, numel(elevation_kernel.az))
            end
        end

        fprintf("completed range %d/%d\n", i_range, i_range_max)
    end
else
    i_range_max = find_range_limit_by_layer(sigma, zz, r_e, theta0, simga_tp, elevation_kernel.range);
    Felevf = cell(i_range_max, numel(elevation_kernel.az));
    Felevid = cell(i_range_max, numel(elevation_kernel.az));
    fprintf("memory_efficient mode: processing %d fine-z layers, %d ranges x %d azimuths\n", numel(zz), i_range_max, numel(elevation_kernel.az))

    for i_z = 1:numel(zz)
        [RR_layer, TT_layer] = build_polar_layer(sigma, zz(i_z), r_e);
        id_theta = find(TT_layer>(theta0-3*simga_tp) & TT_layer<(theta0+3*simga_tp));
        if isempty(id_theta)
            continue
        end

        i_range_max_layer = sum(elevation_kernel.range<=max(RR_layer(id_theta)));
        id_z = floor((i_z-1)/nzz)+1;
        fprintf("layer %d/%d | coarse z %d/%d | active ranges %d/%d\n", i_z, numel(zz), id_z, numel(z), i_range_max_layer, i_range_max)

        for i_range = 1:i_range_max_layer
            r0 = elevation_kernel.range(i_range);
            id_range = id_theta(abs(RR_layer(id_theta)-r0)<=3*simga_r);
            if isempty(id_range)
                continue
            end

            for i_azim = 1:numel(elevation_kernel.az)
                phi0 = elevation_kernel.az(i_azim);
                diff_phi = mod(PP(id_range)-phi0,360);
                diff_phi = min(360-diff_phi, diff_phi);
                id_range_phi = diff_phi<3*simga_tp;
                id_phi = id_range(id_range_phi);

                if isempty(id_phi)
                    continue
                end

                f = exp(-(RR_layer(id_phi)-r0).^2/2/simga_r^2 -(TT_layer(id_phi)-theta0).^2/simga_tp^2 - (diff_phi(id_range_phi)).^2/simga_tp^2 );
                [id_yy,id_xx] = ind2sub(sz_fine_xy,id_phi);
                id_x = floor((id_xx-1)/nxx)+1;
                id_y = floor((id_yy-1)/nyy)+1;
                id_coarse = sub2ind(sz_coarse_grid, id_x, id_y, id_z * ones(size(id_x)));

                Felevf{i_range,i_azim} = [Felevf{i_range,i_azim}; f(:)];
                Felevid{i_range,i_azim} = [Felevid{i_range,i_azim}; uint32(id_coarse(:))];

                if should_report_azimuth_progress(i_azim, numel(elevation_kernel.az))
                    fprintf("  layer %d/%d | range %d/%d | azimuth %d/%d\n", i_z, numel(zz), i_range, i_range_max_layer, i_azim, numel(elevation_kernel.az))
                end
            end
        end
    end

    for i_range = 1:i_range_max
        for i_azim = 1:numel(elevation_kernel.az)
            if isempty(Felevid{i_range,i_azim})
                continue
            end

            [f_keep, id_keep] = collapse_beam_weights(Felevid{i_range,i_azim}, Felevf{i_range,i_azim}, thr_cumsum);
            Felevf{i_range,i_azim} = f_keep;
            Felevid{i_range,i_azim} = uint32([ ...
                i_range*ones(numel(id_keep),1), ...
                i_azim*ones(numel(id_keep),1), ...
                id_keep(:)]);
        end

        fprintf("collapsed range %d/%d\n", i_range, i_range_max)
    end
end

% Indices of F
tmp = reshape(Felevid,[],1);
id=cellfun(@(x) size(x,2),tmp)==3;
Felevidc=cell2mat(tmp(id));

% Value of F
tmp = reshape(Felevf,[],1);
Felevfc=cell2mat(tmp(id));

% Save
save("F_" + elevation_kernel.name + "_" + num2str(max(x) / 1000) + "km.mat", ...
    'Felevidc', 'Felevfc', 'x', 'y', 'z', 'elevation_kernel')
fprintf("saved F_%s_%dkm.mat\n", elevation_kernel.name, round(max(x) / 1000))

end


function [i_range_max, RR, TT] = build_full_polar_grid(sigma, zz, r_e, theta0, simga_tp, ranges)
h_p = (r_e + zz) .* cos(sigma) - r_e;
A = (r_e + zz) .* sin(sigma);
RR = sqrt(h_p.^2 + A.^2);
TT = rad2deg(asin(h_p ./ RR));
id_theta = TT>(theta0-3*simga_tp) & TT<(theta0+3*simga_tp);
i_range_max = sum(ranges<=max(RR(id_theta)));
end


function [RR_layer, TT_layer] = build_polar_layer(sigma, z_layer, r_e)
h_p = (r_e + z_layer) .* cos(sigma) - r_e;
A = (r_e + z_layer) .* sin(sigma);
RR_layer = sqrt(h_p.^2 + A.^2);
TT_layer = rad2deg(asin(h_p ./ RR_layer));
end


function i_range_max = find_range_limit_by_layer(sigma, zz, r_e, theta0, simga_tp, ranges)
rr_max = [];
for i_z = 1:numel(zz)
    [RR_layer, TT_layer] = build_polar_layer(sigma, zz(i_z), r_e);
    id_theta = TT_layer>(theta0-3*simga_tp) & TT_layer<(theta0+3*simga_tp);
    if any(id_theta(:))
        rr_max(end+1) = max(RR_layer(id_theta)); %#ok<AGROW>
    end
end
i_range_max = sum(ranges<=max(rr_max));
end


function [f_keep, id_keep] = collapse_beam_weights(id_coarse, f, thr_cumsum)
[G, id_g] = findgroups(id_coarse);
f_g = splitapply(@sum, f, G);
[f_gs, id_s] = sort(f_g);
id_ss = cumsum(f_gs)./sum(f_gs) > thr_cumsum;
f_keep = f_gs(id_ss) ./ sum(f_gs(id_ss));
id_keep = id_g(id_s(id_ss));
end


function tf = should_report_azimuth_progress(i_azim, n_azim)
step = max(1, ceil(n_azim / 4));
tf = i_azim == 1 || i_azim == n_azim || mod(i_azim, step) == 0;
end
