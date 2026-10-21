function g = create_grid(elevation_kernels, radar_scans, varargin)
%CREATE_GRID Build the Cartesian reconstruction g from radar geometry.
%   GRID = CREATE_GRID(ELEVATION_KERNELS, RADAR_SCANS, ...) returns the
%   shared g used for reconstruction. The g uses x-fastest ordering
%   to match SUB2IND with GRID.SZ = [numel(GRID.X) numel(GRID.Y) numel(GRID.Z)].

% Preserve the time axis before flattening scans for spatial geometry.
g.time = cellfun(@(x) x.datetime, radar_scans(1,:));
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

% Overwrite the g bounds if provided.
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
g.dlat = deg2km(distance(lat_mean,lon_mean,lat_mean+1,lon_mean))*1000;
g.dlon = deg2km(distance(lat_mean,lon_mean,lat_mean,lon_mean+1))*1000;

% 
xmax = (lonmax-lonmin)*g.dlon;
xmax = ceil(xmax/10000)*10000;
ymax = (latmax-latmin)*g.dlat;
ymax = ceil(ymax/10000)*10000;

g.dx = median(diff(elevation_kernels{1}.x));
g.x=(0:ceil(xmax/g.dx)+3)*g.dx;
g.dy = median(diff(elevation_kernels{1}.y));
g.y=(0:ceil(ymax/g.dy)+3)'*g.dy;
g.z=elevation_kernels{1}.z(elevation_kernels{1}.z <= heightmax);
g.dz = median(diff(g.z));

g.sz = [numel(g.x) numel(g.y) numel(g.z)];

% get lat lon coordinates
g.f_lon = @(X) X./g.dlon+lonmin;
g.f_lat = @(Y) Y./g.dlat+latmin;
g.f_X = @(lon) (lon-lonmin)*g.dlon;
g.f_Y = @(lat) (lat-latmin)*g.dlat;

g.lat=g.f_lat(g.y);
g.lon=g.f_lon(g.x);

% for i_r=1:numel(radar_scans)
%     [~, radar_scans{i_r}.radar_lon_id] = min(abs(g.lon(:)' - radar_scans{i_r}.lon(:)),[],2);
%     [~, radar_scans{i_r}.radar_lat_id] = min(abs(g.lat(:)' - radar_scans{i_r}.lat(:)),[],2);
% end

% Load shared water mask and transpose it so the first dimension matches x
% and the second dimension matches y.
water_mask_data = load("data/reference/onwater/Fonwater.mat", "Fonwater");
g.mask_water = logical(water_mask_data.Fonwater({g.f_lat(g.y),g.f_lon(g.x)}))';
% figure; hold on; imagesc(g.f_lon(g.x),g.f_lat(g.y),g.mask_water'); tmp=axis;borders('states','w'); plot(radar_lon, radar_lat,".r"); axis equal tight; axis(tmp);

% Load DEM on the same x-by-y orientation as the reconstruction g.
dem_data = load('data/reference/DEM/dem', 'dem', 'dem_lat', 'dem_lon');
dem = double(dem_data.dem); dem(dem==-9999)=-2;
% demsm = imgaussfilt(dem,[100 100]);
demsmf = griddedInterpolant({dem_data.dem_lat,dem_data.dem_lon},dem,'nearest');
g.dem = demsmf({g.f_lat(g.y),g.f_lon(g.x)})';
% figure; hold on; imagesc(g.f_lon(g.x),g.f_lat(g.y),g.dem); tmp=axis; borders('states','w'); plot(radar_lon, radar_lat,".r"); axis equal tight; axis(tmp); demcmap([0 max(g.dem)]); 
