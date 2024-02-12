% -----------------------------------------------------------------------------------------
% Radar to Geographic Coordinate Transformation
% -----------------------------------------------------------------------------------------
%
% Author: Chisheng Wang
% Affiliation: Shenzhen University
% Creation Date: 2024-01-30
%
% Description:
%   This script transforms radar pixel coordinates to geographic coordinates (latitude,
%   longitude, and altitude). It reads parameters from specific data files, processes
%   radar data, and calculates the corresponding geographic coordinates. The script
%   also includes functionality to output the results in UTM format and visualize them.
%
% Dependencies:
%   - MATLAB 2022
%   - Mapping Toolbox for geospatial operations

%
% Inputs:
%   - Radar data parameters file (e.g., 'data/20191103.slc.par')
%   - PS location data (e.g., loaded from 'data/ps_loc.mat')
%
% Output:
%   - UTM coordinates output in 'output/movingutm.las'
%   - 3D scatter plot of the geographic coordinates
%
% Citations:
%   - Wang, C., L. Chang, X. S. Wang, B. Zhang, and A. Stein (2024), Interferometric Synthetic 
%     Aperture Radar Statistical Inference in Deformation Measurement and Geophysical Inversion: 
%     A Review, IEEE Geoscience and Remote Sensing Magazine, 2-29.
%   - Wang, C., M. Wei, X. Qin, T. Li, S. Chen, C. Zhu, P. Liu, and L. Chang (2024), Three-dimensional
%     lookup table for more precise SAR scatterers positioning in urban scenarios, ISPRS Journal of 
%     Photogrammetry and Remote Sensing, 209, 133-149.
% -----------------------------------------------------------------------------------------

%% Read Parameters
% Define file paths
rslcpar = 'data/20191103.slc.par';
load data/ps_loc.mat;

% Read ellipsoid parameters
ellipsoid = zeros(1, 2);
ellipsoid(1) = readparm(rslcpar, 'earth_semi_major_axis'); % Semi-major axis
ellipsoid(2) = readparm(rslcpar, 'earth_semi_minor_axis'); % Semi-minor axis

% Define constants and read image dimensions
SOL = 299792458;
img_width = readparm(rslcpar, 'range_samples');
img_height = readparm(rslcpar, 'azimuth_lines');

%% Prepare Sensor Parameters
sensor_para.ApproxCenterOriginal.x = readparm(rslcpar, 'center_longitude');
sensor_para.ApproxCenterOriginal.y = readparm(rslcpar, 'center_latitude');
sensor_para.RSR2x = 1 / ((readparm(rslcpar, 'far_range_slc') - readparm(rslcpar, 'near_range_slc')) / (readparm(rslcpar, 'range_samples') * SOL));
sensor_para.t_Range1 = readparm(rslcpar, 'near_range_slc') / SOL;

ell.a = ellipsoid(1);
ell.b = ellipsoid(2);

%% Calculate and Transform
MAXITER = 100;
CRITERPOS = 1e-6;
tic
[x1, y1] = ind2sub([img_height, img_width], ps_loc(:, 1));

xyz = zeros(length(x1), 3);
parfor i = 1:length(x1)
    rn = y1(i);
    ln = x1(i);
    hn = ps_loc(i, 2);
    [orbit_position, orbit_velocity] = intpertsateorbit(rslcpar, ln);
    xyz(i, :) = lph2xyz(orbit_position, orbit_velocity, sensor_para, ln, rn, hn, ell, CRITERPOS, MAXITER);

    if mod(i, 10000) == 0
        fprintf('Progress: %d from %d\n', i, length(x1));
    end
end
toc

%% Translate to Longitude and Latitude
wgs84 = wgs84Ellipsoid;
[lat_s, lon_s, alt] = ecef2geodetic(wgs84, xyz(:, 1), xyz(:, 2), xyz(:, 3));

%% Write to UTM and Visualization
[utm_x, utm_y] = ll2utm(lat_s, lon_s, 'wgs84', 50);
moving1 = pointCloud([utm_x, utm_y, alt]);

% Ensure output directory exists
if ~exist('output', 'dir')
    mkdir('output');
end

lasWriter = lasFileWriter('output/movingutm.las', 'PointDataFormat', 1);
writePointCloud(lasWriter, moving1);

% Plotting
figure;
scatter3(lat_s, lon_s, alt, 3, alt, 'filled');
colormap(jet);
xlabel('latitude');
ylabel('longitude');
c = colorbar;
ylabel(c, 'elevation (m)');
