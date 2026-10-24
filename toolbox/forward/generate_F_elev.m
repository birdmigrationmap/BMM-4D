
% Script to generate missing F_elev kernels from saved scan patterns.
% Run this from the project root after updating data/reference/scan_pattern_new.mat.

% Run from server with
% module load matlab/R2021a
% matlab -nodisplay -nosplash -nodesktop -nojvm -r "elev_iter = $1 ;run('generate_F_elev.m');exit;"

%% Coarse grid of the kernel set ("150km" or "50km")
kernel_set = "50km";
dz = 100;
z = reshape(dz/2:dz:3000-dz/2, 1, 1, []);
nzz = 5;
if kernel_set == "150km"
    dx = 2000;
    xy_max = 150000;
    nxx = 40;
else
    dx = 500;
    xy_max = 50000;
    nxx = 10;
end
dy = dx;
nyy = nxx;
x = reshape(-xy_max:dx:xy_max, 1, []);
y = reshape(-xy_max:dy:xy_max, [], 1);

%% Compute the missing kernels
load(data_path("reference", "scan_pattern_new.mat"), "dfnew")

% On the server, elev_iter selects one kernel per job; otherwise compute all.
if exist("elev_iter", "var")
    kernel_ids = elev_iter;
else
    kernel_ids = 1:numel(dfnew);
end
for i_kernel = kernel_ids
    compute_F_elev(x, y, z, nxx, nyy, nzz, dfnew(i_kernel), 0.1)
end
