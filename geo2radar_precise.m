% -----------------------------------------------------------------------------------------
% Geographic to Radar Pixel Coordinate Transformation
% -----------------------------------------------------------------------------------------
%
% Author: Chisheng Wang
% Affiliation: Shenzhen University
% Creation Date: 2024-01-30
%
% Description:
%   This script transforms geographic coordinates to radar pixel coordinates. 
%   It involves reading parameters from satellite orbit data, processing point cloud data from a LAS file,
%   applying a transformation matrix, and calculating corresponding radar pixel coordinates.
%   The script outputs a 3D lookup table and visualizes this transformation.
%
% Dependencies:
%   - MATLAB 2022
%   - Mapping Toolbox (for geospatial operations)
%
% Inputs:
%   - Satellite orbit parameters (e.g., from 'data/20191103.slc.par')
%   - Point cloud data (from 'data/Merged_Data_sub2m.las')
%   - Transformation matrix (loaded from 'data/Tmatrix.mat')
%
% Outputs:
%   - 3D lookup table for the radar pixel coordinates ('output/3Dlookuptable_precise.mat')
%   - Visualizations of the 3D lookup table and radar-coded height
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
rslcpar = 'data/20191103.slc.par';
load data/Tmatrix.mat;
ellipsoid(1) = readparm(rslcpar, 'earth_semi_major_axis');
ellipsoid(2) = readparm(rslcpar, 'earth_semi_minor_axis');
SOL = 299792458;
orbit = allsateorbit(rslcpar);

%% Generate 3D Lookup Table
lasfile = 'data/Merged_Data_sub2m.las';
lasReader = lasFileReader(lasfile);
[ptCloud, atri] = readPointCloud(lasReader, 'Attributes', 'Classification');
ind = 1:atri.Count;
pt1 = [ptCloud.Location(ind, 1:3), 1+0*ptCloud.Location(ind, 1)];
pt1 = pt1';
pt2 = T\pt1;
sp = 1;

% Geographic conversion
[lat_l, lon_l] = utm2ll(pt2(1, 1:sp:end), pt2(2, 1:sp:end), 50);
wgs84 = wgs84Ellipsoid;
[xyz_t(:, 1), xyz_t(:, 2), xyz_t(:, 3)] = geodetic2ecef(wgs84, lat_l', lon_l', pt2(3, 1:sp:end)');

%% Transform to Radar Pixel Coordinates
MAXITER = 100;
CRITERPOS = 1e-6;
coeff = orbit_info(rslcpar);
ta1 = readparm(rslcpar, 'start_time');
PRF = readparm(rslcpar, 'prf');
midLine = round(readparm(rslcpar, 'azimuth_lines') / 2);
coeff.timeAzimuth = (midLine - 1) / PRF + ta1;
tr1 = readparm(rslcpar, 'near_range_slc') / SOL;
RSR = 1 / ((readparm(rslcpar, 'far_range_slc') - readparm(rslcpar, 'near_range_slc')) / (readparm(rslcpar, 'range_samples') * SOL));

parfor i = 1:size(xyz_t, 1)
    xyz_vec = [xyz_t(i, 1), xyz_t(i, 2), xyz_t(i, 3)];
    [t, solution] = xyz2t(xyz_vec, coeff, MAXITER, CRITERPOS);
    t_azi = t.timeAzimuth;
    t_ran = t.timeRange;
    line1(i) = 1 + PRF * (t_azi - ta1);
    pixel1(i) = 1 + RSR * (t_ran - tr1);

    if mod(i, 10000) == 0
        fprintf('Progress: %d from %d \n', i, size(xyz_t, 1))
    end
end

save output/3Dlookuptable_precise.mat line1 pixel1

%% Visualization
% Azimuth
figure;
scatter3(pt1(1, :), pt1(2, :), pt1(3, :), 1, line1);
colormap(jet); colorbar; title('3DLUT-azimuth');

% Range
figure;
scatter3(pt1(1, :), pt1(2, :), pt1(3, :), 1, pixel1);
colormap(jet); colorbar; title('3DLUT-range');

%% Process for Maximum Height Visualization
sar_xy_c = pt2(:, 1:sp:end)';
img_width = readparm(rslcpar, 'range_samples');
img_height = readparm(rslcpar, 'azimuth_lines');

% Initialize height matrices
hgt_img_max = zeros(img_height, img_width) - 1000;
hgt_img_min = zeros(img_height, img_width) + 10000;
num_img = zeros(img_height, img_width);

% Populate height matrices
for ii = 1:length(line1)
    ix = round(line1(ii));
    iy = round(pixel1(ii));
    if ix > 0 && iy > 0 && ix <= img_height && iy <= img_width
        num_img(ix, iy) = num_img(ix, iy) + 1;
        hgt_img_max(ix, iy) = max(hgt_img_max(ix, iy), sar_xy_c(ii, 3));
        hgt_img_min(ix, iy) = min(hgt_img_min(ix, iy), sar_xy_c(ii, 3));
    end
end

%% Display Maximum Height Visualization
figure;
subplot(1, 2, 1);
imagesc(hgt_img_max);
title('Radar-coded Height');
colormap(jet); clim([-10 200]); axis equal; axis off;

% Display corresponding MLI image
mli = imread('data/20191103.slc.mli.bmp');
subplot(1, 2, 2);
imshow(mli);
title('MLI Image');
