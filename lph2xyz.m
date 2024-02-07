function XYZ = lph2xyz(orbit_position, orbit_velocity, sensor_para, line, pixel, h, ell, CRITERPOS, MAXITER)
    % This script has been modified from the original version available at:
    % https://github.com/maerabubu/InSAR-Baseline-Estimation/blob/master/lph2xyz.m

    % Function to calculate XYZ coordinates from radar line and pixel values
    %
    % Inputs:
    %   orbit_position - Position of the satellite
    %   orbit_velocity - Velocity of the satellite
    %   sensor_para - Sensor parameters
    %   line, pixel - Radar line and pixel values
    %   h - Altitude
    %   ell - Ellipsoid parameters
    %   CRITERPOS - Convergence criterion for position calculation
    %   MAXITER - Maximum number of iterations for convergence
    %
    % Output:
    %   XYZ - Calculated XYZ coordinates

    % Default values and checks
    if nargin < 5, h = zeros(size(line)); end
    if nargin < 6
        ell.a = 6378137.0 + h; % semimajor axis WGS84
        ell.b = 6356752.3142451794975639665996337 + h; % semiminor axis WGS84
    else
        ell.a = ell.a + h; % semimajor axis WGS84
        ell.b = ell.b + h; % semiminor axis WGS84
    end
    if nargin < 7, CRITERPOS = 1e-6; end
    if nargin < 8, MAXITER = 20; end
    if size(line) ~= size(pixel), error('line and pixel should have same dimensions!'); end

    % Constants
    VLIGHT = 299792458; % Speed of light in m/s

    % Extracting orbit position and velocity
    pos.X = orbit_position(1); pos.Y = orbit_position(2); pos.Z = orbit_position(3);
    posv.X = orbit_velocity(1); posv.Y = orbit_velocity(2); posv.Z = orbit_velocity(3);

    % Sensor parameters
    RSR2x = sensor_para.RSR2x;
    tr1 = sensor_para.t_Range1;  
    r_time = tr1 + (pixel - 1) / RSR2x;

    % Initialize geodetic to ECEF conversion
    wgs84 = wgs84Ellipsoid;
    [X, Y, Z] = geodetic2ecef(wgs84, sensor_para.ApproxCenterOriginal.y, sensor_para.ApproxCenterOriginal.x, h);

    % Prepare matrices for iterations
    [row, col] = size(line);
    equationset = zeros(3, 1, row, col);
    solxyz = zeros(3, row, col);
    partialsxyz = zeros(3, 3, row, col);

    % Iterative solution
    for iter = 1:MAXITER
        dsat_P.x = X - pos.X;
        dsat_P.y = Y - pos.Y;
        dsat_P.z = Z - pos.Z;

        % Set up equations for the iterative process
        equationset(1, 1, :, :) = (posv.X * dsat_P.x + posv.Y * dsat_P.y + posv.Z * dsat_P.z) * (-1);
        equationset(2, 1, :, :) = (dsat_P.x * dsat_P.x + dsat_P.y * dsat_P.y + dsat_P.z * dsat_P.z - (VLIGHT * r_time) * (VLIGHT * r_time)) * (-1);
        equationset(3, 1, :, :) = ((X * X + Y * Y) / (ell.a * ell.a) + (Z / ell.b) * (Z / ell.b) - 1) * (-1);

        % Calculate partial derivatives
        partialsxyz(1, 1, :, :) = posv.X;
        partialsxyz(1, 2, :, :) = posv.Y;
        partialsxyz(1, 3, :, :) = posv.Z;
        partialsxyz(2, 1, :, :) = 2 * dsat_P.x;
        partialsxyz(2, 2, :, :) = 2 * dsat_P.y;
        partialsxyz(2, 3, :, :) = 2 * dsat_P.z;
        partialsxyz(3, 1, :, :) = (2 * X) / (ell.a * ell.a);
        partialsxyz(3, 2, :, :) = (2 * Y) / (ell.a * ell.a);
        partialsxyz(3, 3, :, :) = (2 * Z) / (ell.b * ell.b);

        % Solve the equations
        for i = 1:row
            for j = 1:col
                solxyz(:, i, j) = (partialsxyz(:, :, i, j)) \ equationset(:, :, i, j);
            end
        end

        % Update XYZ coordinates
        X = shiftdim(solxyz(1, :, :)) + X;   
        Y = shiftdim(solxyz(2, :, :)) + Y;   
        Z = shiftdim(solxyz(3, :, :)) + Z;   

        % Check for convergence
        if abs(solxyz) < CRITERPOS
            break;
        end
        if iter > MAXITER
            disp('Reached the max iteration number in lph2xyz()');
            break;
        end
    end

    % Final XYZ coordinates
    XYZ = [X Y Z];
end
