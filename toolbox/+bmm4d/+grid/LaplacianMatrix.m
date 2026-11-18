function L3D = LaplacianMatrix(g, varargin)
%LAPLACIANMATRIX Build the 3D smoothness operator used in inversion.
%   L3D = LAPLACIANMATRIX(G, ...) returns a sparse 3D Laplacian on the
%   reconstruction g G, optionally disconnecting water and land cells.
%   The second differences are in physical units: the operator
%   approximates l_h^2 (d2/dx2 + d2/dy2) + l_z^2 d2/dz2, so it does not
%   depend on the grid spacing, and the horizontal-to-vertical balance is
%   set by l_h / l_z ('length_scale_m', default [2000 100] m, equal to the
%   2 km x 100 m grid where the weights are 1).

p = inputParser;
addParameter(p,'water',true); 
addParameter(p,'coeff_horver',0); % (+) more horizontal smoother, (-) more vertical smoothing, 0 equal smoothing. -1 only vertical, 0.5 only horizontal
addParameter(p,'length_scale_m',[2000 100]); % [l_h l_z] of the physical-unit operator
% addParameter(p,'order_derivative_smoothness', 2);
parse(p,varargin{:});


if p.Results.coeff_horver>0.5 || p.Results.coeff_horver<-1
    warning("coeff_horver should be between -1 (only vertical) and 0.5 (only horizontal)")
end

T1 = bmm4d.grid.LaplacianMatrix1D(g.sz(1));
% T1 = T1 + .7*LaplacianMatrix1D(g.sz(1),3);
% T1 = T1 + .3*LaplacianMatrix1D(g.sz(1),20);
T1 = T1 - diag(sum(T1,2));

T2 = bmm4d.grid.LaplacianMatrix1D(g.sz(2));
% T2 = T2 + .7*LaplacianMatrix1D(g.sz(2),3);
% T2 = T2 + .3*LaplacianMatrix1D(g.sz(2),20);
T2 = T2 - diag(sum(T2,2));


T3 = bmm4d.grid.LaplacianMatrix1D(g.sz(3));
% T3 = T3 + 1*LaplacianMatrix1D(g.sz(3),3);
% T3 = T3 + .5*LaplacianMatrix1D(g.sz(3),20);
T3 = T3 - diag(sum(T3,2));

% Physical-unit second differences: (l / spacing)^2 along each axis.
length_scale_m = p.Results.length_scale_m;
T1 = (length_scale_m(1) / g.dx)^2 .* T1;
T2 = (length_scale_m(1) / g.dy)^2 .* T2;
T3 = (length_scale_m(2) / g.dz)^2 .* T3;

% d = p.Results.order_derivative_smoothness;
% T1 = [difference_operator(g.sz(1), d); zeros(d,g.sz(1))];
% T2 = [difference_operator(g.sz(2), d); zeros(d,g.sz(2))];
% T3 = [difference_operator(g.sz(3), d); zeros(d,g.sz(3))];


L2D = kron(T2, speye(g.sz(1))) + kron(speye(g.sz(2)), T1);

if p.Results.water
    onwater = bmm4d.grid.get_water_mask(g);
    % figure; hold on; imagesc(g.f_lon(g.x),g.f_lat(g.y),onwater'); tmp=axis;borders('states','w'); axis equal tight; axis(tmp);
    
    L2Dwater=L2D;
    [i,j] = find(L2Dwater);
    id = onwater(i) ~= onwater(j);
    % Cut only the land-water pairs; indexing L(i,j) with vectors would
    % zero the whole submatrix and also cut along-shore links.
    L2Dwater(sub2ind(size(L2Dwater), i(id), j(id))) = 0;
    L2Dwater = L2Dwater - diag(sum(L2Dwater,2));
    % figure; hold on; imagesc(g.lon,g.lat,reshape(diag(L2Dwater),g.sz(1),g.sz(2))); tmp=axis;borders('states','w'); axis equal tight; axis(tmp);
    
    L2D = L2Dwater;
end
L3D = (kron(T3, speye(prod(g.sz(1:2))))*(1-p.Results.coeff_horver*2) + kron(speye(g.sz(3)), L2D)*(1+p.Results.coeff_horver))/6;
end
