function beam = gd_beam(varargin)
%GD_BEAM Gate and bird model shared by the simulator and the estimator.
%   A gate at range r collects the birds within a disc of angular radius
%   SAMPLE_RADIUS * WIDTH around the beam axis and a range interval of
%   GATE_LENGTH. Each bird adds sigma * w to the gate, with w the two-way
%   Gaussian beam weight exp(-8 ln2 a^2 / width^2) at angular offset a, so
%   the reflectivity is eta = sum(sigma * w) / V_eff with
%   V_eff = GATE_LENGTH r^2 pi width^2 / (8 ln2), the integral of w; then
%   E[eta] = density * mean(sigma). Bird cross-sections are lognormal with
%   mean SIGMA_MEAN (cm^2, 11 as in the density conversion) and log
%   standard deviation SIGMA_LOG_SD.
%   Fields: width_rad, gate_length_m, sample_radius, sigma_mean_cm2,
%   sigma_log_sd, and handles volume_eff_km3(r_m) and volume_sample_km3(r_m).

p = inputParser;
addParameter(p, "width_deg", 0.95);
addParameter(p, "gate_length_m", 250);
addParameter(p, "sample_radius", 1.2);       % in beam widths; w < 4e-4 beyond
addParameter(p, "sigma_mean_cm2", 11);
addParameter(p, "sigma_log_sd", 1);
parse(p, varargin{:});
beam = p.Results;
beam.width_rad = deg2rad(beam.width_deg);
beam.volume_eff_km3 = @(r_m) beam.gate_length_m .* r_m .^ 2 .* pi .* beam.width_rad .^ 2 ./ (8 * log(2)) / 1e9;
beam.volume_sample_km3 = @(r_m) beam.gate_length_m .* r_m .^ 2 .* pi .* (beam.sample_radius * beam.width_rad) .^ 2 / 1e9;
end
