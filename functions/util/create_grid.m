function grid = create_grid(elevation_kernels, radar_scans, varargin)
%CREATE_GRID Build the Cartesian reconstruction grid from radar geometry.
%   GRID = CREATE_GRID(ELEVATION_KERNELS, RADAR_SCANS, ...) returns the
%   shared grid used for reconstruction. The grid uses x-fastest ordering
%   to match SUB2IND with GRID.SZ = [numel(GRID.X) numel(GRID.Y) numel(GRID.Z)].

empty_scan_id = cellfun(@isempty, radar_scans);
if any(empty_scan_id)
    warning("Empty radar_scans at index: " + join(string(find(empty_scan_id)), ", "))
end
radar_scans = radar_scans(~empty_scan_id);

radar_lon = cellfun(@(x) x.lon, radar_scans(:,1));
radar_lat = cellfun(@(x) x.lat, radar_scans(:,1));

tmp_dlat = deg2km(distance(mean(radar_lat),mean(radar_lon),mean(radar_lat)+1,mean(radar_lon)))*1000;
tmp_dlon = deg2km(distance(mean(radar_lat),mean(radar_lon),mean(radar_lat),mean(radar_lon)+1))*1000;

radar_dlat = max(elevation_kernels{1}.y)/tmp_dlat; 
radar_dlon = max(elevation_kernels{1}.x)/tmp_dlon; 

lonmin=min(radar_lon)-radar_dlon;
lonmax=max(radar_lon)+radar_dlon;
latmin=min(radar_lat)-radar_dlat;
latmax=max(radar_lat)+radar_dlat;

% Overwrite the grid bounds if provided.
if nargin>2
    lonmin=varargin{1};
end
if nargin>3
    lonmax=varargin{2};
end
if nargin>4
    latmin=varargin{3};
end
if nargin>5
    latmax=varargin{4};
end
if nargin>6
    heightmax=varargin{5};
else
    heightmax=max(elevation_kernels{1}.z);
end

% Find how much 1km=degree
lat_mean=(latmax+latmin)/2;
lon_mean=(lonmax+lonmin)/2;
grid.dlat = deg2km(distance(lat_mean,lon_mean,lat_mean+1,lon_mean))*1000;
grid.dlon = deg2km(distance(lat_mean,lon_mean,lat_mean,lon_mean+1))*1000;

% 
xmax = (lonmax-lonmin)*grid.dlon;
xmax = ceil(xmax/10000)*10000;
ymax = (latmax-latmin)*grid.dlat;
ymax = ceil(ymax/10000)*10000;

grid.dx = median(diff(elevation_kernels{1}.x));
grid.x=(0:ceil(xmax/grid.dx)+3)*grid.dx;
grid.dy = median(diff(elevation_kernels{1}.y));
grid.y=(0:ceil(ymax/grid.dy)+3)'*grid.dy;
grid.z=elevation_kernels{1}.z(elevation_kernels{1}.z <= heightmax);
grid.dz = median(diff(grid.z));

grid.sz = [numel(grid.x) numel(grid.y) numel(grid.z)];

% get lat lon coordinates
grid.f_lon = @(X) X./grid.dlon+lonmin;
grid.f_lat = @(Y) Y./grid.dlat+latmin;
grid.f_X = @(lon) (lon-lonmin)*grid.dlon;
grid.f_Y = @(lat) (lat-latmin)*grid.dlat;

grid.lat=grid.f_lat(grid.y);
grid.lon=grid.f_lon(grid.x);

% for i_r=1:numel(radar_scans)
%     [~, radar_scans{i_r}.radar_lon_id] = min(abs(grid.lon(:)' - radar_scans{i_r}.lon(:)),[],2);
%     [~, radar_scans{i_r}.radar_lat_id] = min(abs(grid.lat(:)' - radar_scans{i_r}.lat(:)),[],2);
% end

grid.time = cellfun(@(x) x.datetime, radar_scans(1,:));

% Load shared water mask and transpose it so the first dimension matches x
% and the second dimension matches y.
water_mask_data = load("data/reference/onwater/Fonwater.mat", "Fonwater");
grid.mask_water = logical(water_mask_data.Fonwater({grid.f_lat(grid.y),grid.f_lon(grid.x)}))';
% figure; hold on; imagesc(grid.f_lon(grid.x),grid.f_lat(grid.y),grid.mask_water'); tmp=axis;borders('states','w'); plot(radar_lon, radar_lat,".r"); axis equal tight; axis(tmp);

% Load DEM on the same x-by-y orientation as the reconstruction grid.
dem_data = load('data/reference/DEM/dem', 'dem', 'dem_lat', 'dem_lon');
dem = double(dem_data.dem); dem(dem==-9999)=-2;
% demsm = imgaussfilt(dem,[100 100]);
demsmf = griddedInterpolant({dem_data.dem_lat,dem_data.dem_lon},dem,'nearest');
grid.dem = demsmf({grid.f_lat(grid.y),grid.f_lon(grid.x)})';
% figure; hold on; imagesc(grid.f_lon(grid.x),grid.f_lat(grid.y),grid.dem); tmp=axis; borders('states','w'); plot(radar_lon, radar_lat,".r"); axis equal tight; axis(tmp); demcmap([0 max(grid.dem)]); 
