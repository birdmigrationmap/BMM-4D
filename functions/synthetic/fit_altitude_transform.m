function [obs_values_gaussian, trans_finv, trans_ffor] = fit_altitude_transform(obs_values, beam_height_m, grid, vp, varargin)
%FIT_ALTITUDE_TRANSFORM Transform density residuals to Gaussian space by altitude.
% A separate transform is fitted for each altitude bin because the possible
% density range changes with altitude. Neighboring bins are included to
% stabilize the kernel distribution and reduce discontinuities.

p = inputParser;
addParameter(p, "min_points", 100);
addParameter(p, "n_axis", 1001);
parse(p, varargin{:});

gz = squeeze(grid.z);
obs_values_detrended = obs_values(:) - interp1(gz, vp(:), beam_height_m(:), "linear");
altitude_bin = discretize(beam_height_m(:), [-grid.dz/2; gz(:)] + grid.dz/2);

dens_axis = nan(p.Results.n_axis, grid.sz(3));
denst_axis = nan(p.Results.n_axis, grid.sz(3));
for i_z = 1:grid.sz(3)
    id_z = altitude_bin >= i_z - 1 & altitude_bin <= i_z + 1 & isfinite(obs_values_detrended);
    if sum(altitude_bin == i_z) > p.Results.min_points
        pd = fitdist(obs_values_detrended(id_z), "Kernel");
        [~, dens_axis(:, i_z)] = discretize(obs_values_detrended(id_z), p.Results.n_axis - 1);
        denst_axis(:, i_z) = norminv(cdf(pd, dens_axis(:, i_z)));
    end
end

z_axis = reshape(repmat(gz(:)', p.Results.n_axis, 1), [], 1);
id = isfinite(dens_axis(:)) & isfinite(denst_axis(:));

trans_ffor = scatteredInterpolant(dens_axis(id), z_axis(id), denst_axis(id), "linear", "none");
obs_values_gaussian = trans_ffor(obs_values_detrended, beam_height_m(:));

trans_finv = scatteredInterpolant(denst_axis(id), z_axis(id), dens_axis(id), "linear", "nearest");

if false
    % Check the detrending and Gaussian transformation.
    figure;
    tiledlayout("flow");
    nexttile; scatter(obs_values_detrended, beam_height_m, 2, "."); xlabel("Detrended density"); ylabel("Altitude (m)")
    nexttile; scatter(obs_values_gaussian, beam_height_m, 2, "."); xlabel("Gaussian density"); ylabel("Altitude (m)")
    nexttile; histogram(obs_values_gaussian); xlabel("Gaussian density")

    % Inspect the kernel distribution fitted within one altitude bin.
    i_z = round(grid.sz(3)/2);
    id_z = altitude_bin >= i_z - 1 & altitude_bin <= i_z + 1 & isfinite(obs_values_detrended);
    pd = fitdist(obs_values_detrended(id_z), "Kernel");
    figure; histogram(obs_values_detrended(id_z), Normalization="pdf"); hold on; plot(pd)
end
end
