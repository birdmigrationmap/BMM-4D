function beam = beam(varargin)
%BEAM Gate and bird model of the gate density estimator.
%   Each bird adds sigma * w to a gate, with w the two-way beam weight
%   f^4 W^2: f^4 = exp(-8 ln2 a^2 / width^2) at angular offset a from the
%   beam axis and W^2 = exp(-x^2 / (2 sigma_r^2)) at range offset x from the
%   gate centre, sigma_r = 0.35 GATE_LENGTH. The reflectivity is
%   eta = sum(sigma * w) / V_eff with V_eff = r^2 pi width width_az / (8 ln2)
%   * sqrt(2 pi) sigma_r, the integral of w, so that E[eta] = density *
%   mean(sigma). Birds are counted within SAMPLE_RADIUS beam widths of the
%   axis and RANGE_REACH sigma_r of the gate centre (f^4 < 4e-4 and
%   W^2 < 0.012 beyond), a counting volume V_sample.
%   Bird cross-sections are lognormal with mean SIGMA_MEAN (cm^2, 11 as in
%   the density conversion) and log standard deviation SIGMA_LOG_SD. The
%   default 1.25 is an assumption, not a calibration: the best-fitting
%   spread where gates hold 10-30 birds (studies/gate_density, GD_G); it
%   sets the absolute level (x0.92 for 1, x1.27 for 1.75), not the bias
%   with distance.
%   The antenna rotates while a radial is collected, which widens the
%   effective beam in azimuth only: WIDTH_AZ = AZIMUTH_BROADENING * WIDTH,
%   1.088 for WSR-88D super-resolution radials (Torres and Curtis 2007).
%   Fields: width_rad, width_az_rad, gate_length_m, range_sd_m,
%   sample_radius, range_reach, sigma_mean_cm2, sigma_log_sd, and handles
%   volume_eff_km3(r_m) and volume_sample_km3(r_m).

p = inputParser;
addParameter(p, "width_deg", 0.95);
addParameter(p, "gate_length_m", 250);
addParameter(p, "sample_radius", 1.2);         % in beam widths; f^4 < 4e-4 beyond
addParameter(p, "range_reach", 3);             % in range standard deviations; W^2 < 0.012 beyond
addParameter(p, "azimuth_broadening", 1.088);  % effective / antenna beamwidth in azimuth
addParameter(p, "sigma_mean_cm2", 11);
addParameter(p, "sigma_log_sd", 1.25);
parse(p, varargin{:});
beam = p.Results;
beam.width_rad = deg2rad(beam.width_deg);
beam.width_az_rad = beam.azimuth_broadening * beam.width_rad;
beam.range_sd_m = 0.35 * beam.gate_length_m;
beam.volume_eff_km3 = @(r_m) r_m .^ 2 .* pi .* beam.width_rad .* beam.width_az_rad ./ (8 * log(2)) .* sqrt(2 * pi) .* beam.range_sd_m / 1e9;
beam.volume_sample_km3 = @(r_m) r_m .^ 2 .* pi .* beam.sample_radius ^ 2 * beam.width_rad .* beam.width_az_rad .* 2 .* beam.range_reach .* beam.range_sd_m / 1e9;
end
