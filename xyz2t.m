function [time, solution] = xyz2t(xyz_vec, coeff, MAXITER, CRITERTIM)
    % This script has been modified from the original version available at:
    % https://github.com/maerabubu/InSAR-Baseline-Estimation/blob/master/xyz2t.m
    % Function to calculate radar signal travel time
    %
    % Inputs:
    %   xyz_vec - Vector containing X, Y, Z coordinates on Earth's surface
    %   coeff - Coefficients containing satellite information
    %   MAXITER - Maximum number of iterations for convergence
    %   CRITERTIM - Convergence criterion for time calculation
    %
    % Outputs:
    %   time - Structure containing timeAzimuth and timeRange
    %   solution - Converged solution after iterations

    % Constants
    SOL = 299792458; % Speed of light in m/s

    % Default values for MAXITER and CRITERTIM if not provided
    if nargin < 4
        CRITERTIM = 10e-10;
    end
    if nargin < 3
        MAXITER = 20;
    end

    % Extract XYZ coordinates
    pointOnEllips.X = xyz_vec(:, 1);
    pointOnEllips.Y = xyz_vec(:, 2);
    pointOnEllips.Z = xyz_vec(:, 3);

    % Initialize timeAzimuth from coefficients
    timeAzimuth = coeff.timeAzimuth;

    % Iterative solution
    for iter = 1:MAXITER
        % Interpolate satellite position, velocity, and acceleration
        satellitePosition.X = interp1(coeff.time, coeff.state_vector_positions(:, 1), timeAzimuth);
        satellitePosition.Y = interp1(coeff.time, coeff.state_vector_positions(:, 2), timeAzimuth);
        satellitePosition.Z = interp1(coeff.time, coeff.state_vector_positions(:, 3), timeAzimuth);

        satelliteAcceleration.X = interp1(coeff.time, coeff.state_vector_acceleration(:, 1), timeAzimuth);
        satelliteAcceleration.Y = interp1(coeff.time, coeff.state_vector_acceleration(:, 2), timeAzimuth);
        satelliteAcceleration.Z = interp1(coeff.time, coeff.state_vector_acceleration(:, 3), timeAzimuth);

        satelliteVelocity.X = interp1(coeff.time, coeff.state_vector_velocity(:, 1), timeAzimuth);
        satelliteVelocity.Y = interp1(coeff.time, coeff.state_vector_velocity(:, 2), timeAzimuth);
        satelliteVelocity.Z = interp1(coeff.time, coeff.state_vector_velocity(:, 3), timeAzimuth);

        % Calculate delta position
        delta.X = pointOnEllips.X - satellitePosition.X;
        delta.Y = pointOnEllips.Y - satellitePosition.Y;
        delta.Z = pointOnEllips.Z - satellitePosition.Z;

        % Update solution
        solution{iter} = -(satelliteVelocity.X * delta.X + satelliteVelocity.Y * delta.Y + satelliteVelocity.Z * delta.Z) / ...
                         (satelliteAcceleration.X * delta.X + satelliteAcceleration.Y * delta.Y + satelliteAcceleration.Z * delta.Z - ...
                          satelliteVelocity.X^2 - satelliteVelocity.Y^2 - satelliteVelocity.Z^2);

        % Update timeAzimuth
        timeAzimuth = timeAzimuth + solution{iter};

        % Check for convergence
        if max(abs(solution{iter})) < CRITERTIM
            break;
        end
        if iter >= MAXITER
            disp('ERROR: Maximum iterations reached without convergence');
            break;
        end
    end

    % Interpolate final satellite position
    satellitePosition.X = interp1(coeff.time, coeff.state_vector_positions(:, 1), timeAzimuth, 'spline');
    satellitePosition.Y = interp1(coeff.time, coeff.state_vector_positions(:, 2), timeAzimuth, 'spline');
    satellitePosition.Z = interp1(coeff.time, coeff.state_vector_positions(:, 3), timeAzimuth, 'spline');

    % Calculate delta position for final time
    delta.X = pointOnEllips.X - satellitePosition.X;
    delta.Y = pointOnEllips.Y - satellitePosition.Y;
    delta.Z = pointOnEllips.Z - satellitePosition.Z;

    % Calculate timeRange
    timeRange = sqrt(delta.X.^2 + delta.Y.^2 + delta.Z.^2) / SOL;

    % Output time structure
    time.timeAzimuth = timeAzimuth;
    time.timeRange = timeRange;
end
