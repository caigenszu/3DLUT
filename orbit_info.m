
% Open the text file

function [coeff] = orbit_info(slcpar)
fid = fopen(slcpar, 'r');
% Initialize an empty cell array to store the lines
lines = {};

% Read the file line by line
while ~feof(fid)
    lines{end+1, 1} = fgetl(fid);
end

% Close the file
fclose(fid);

% Extract start time
center_time_line = lines{contains(lines, 'center_time')};
center_time = sscanf(center_time_line, 'center_time: %f');

prf_line = lines{contains(lines, 'prf')};
prf= sscanf(prf_line, 'prf: %f');



% Extract number of state vectors
num_state_vectors_line = lines{contains(lines, 'number_of_state_vectors')};
num_state_vectors = sscanf(num_state_vectors_line, 'number_of_state_vectors: %d');

% Extract time interval
time_interval_line = lines{contains(lines, 'state_vector_interval')};
time_interval = sscanf(time_interval_line, 'state_vector_interval: %f');

% Extract time of first state vector
first_state_vector_time_line = lines{contains(lines, 'time_of_first_state_vector')};
first_state_vector_time = sscanf(first_state_vector_time_line, 'time_of_first_state_vector: %f');

% Initialize arrays to store state vector positions and corresponding times
state_vector_positions = zeros(num_state_vectors, 3);
state_vector_times = zeros(num_state_vectors, 1);

% Extract state vector positions and corresponding times
for i = 1:num_state_vectors
    state_vector_position_line = lines{contains(lines, ['state_vector_position_' num2str(i)])};
    state_vector_positions(i, :) = sscanf(state_vector_position_line, ['state_vector_position_' num2str(i) ': %f %f %f']);
    state_vector_velocity_line = lines{contains(lines, ['state_vector_velocity_' num2str(i)])};
    state_vector_velocity(i, :) = sscanf(state_vector_velocity_line, ['state_vector_velocity_' num2str(i) ': %f %f %f']);    
    % Calculate time for each state vector
    state_vector_times(i) = first_state_vector_time + (i-1)*time_interval; % Note: 10s is the state_vector_interval
end
coeff.X=state_vector_positions(:, 1);
coeff.Y=state_vector_positions(:, 2);
coeff.Z=state_vector_positions(:, 3);
coeff.time=state_vector_times;
coeff.timeAzimuth=center_time;
coeff.state_vector_positions=state_vector_positions;
%% estimate velocity
precorb(:,1)=state_vector_times;
precorb(:,2:4)=state_vector_positions(:, 1:3);

% compute normalization factors
px = precorb(:,1); % time
f  = min(px);
g  = (max(px)-min(px));
px = (px-f)/g;
polyDegree = 2;

coef_x1 = (polyfit(px,precorb(:,2),polyDegree));
a = coef_x1(3);
b = coef_x1(2);
c = coef_x1(1);
coef_x = [c/(g^2) b/g-(2*c*f)/(g^2) a-b*f/g+c*(f/g)^2];

coef_y1 = (polyfit(px,precorb(:,3),polyDegree));
a = coef_y1(3);
b = coef_y1(2);
c = coef_y1(1);
coef_y = [c/(g^2) b/g-(2*c*f)/(g^2) a-b*f/g+c*(f/g)^2];

coef_z1 = (polyfit(px,precorb(:,4),polyDegree));
a = coef_z1(3);
b = coef_z1(2);
c = coef_z1(1);
coef_z = [c/(g^2) b/g-(2*c*f)/(g^2) a-b*f/g+c*(f/g)^2];

vel_x = polyval(polyder(coef_x),precorb(:,1));
vel_y = polyval(polyder(coef_y),precorb(:,1));
vel_z = polyval(polyder(coef_z),precorb(:,1));
state_vector_velocity = [vel_x vel_y vel_z];
coeff.state_vector_velocity =state_vector_velocity;
%% estimate acceleration
precorb(:,1)=state_vector_times;
precorb(:,2:4)=state_vector_velocity(:, 1:3);

% compute normalization factors
px = precorb(:,1); % time
f  = min(px);
g  = (max(px)-min(px));
px = (px-f)/g;
polyDegree = 2;

coef_x1 = (polyfit(px,precorb(:,2),polyDegree));
a = coef_x1(3);
b = coef_x1(2);
c = coef_x1(1);
coef_x = [c/(g^2) b/g-(2*c*f)/(g^2) a-b*f/g+c*(f/g)^2];

coef_y1 = (polyfit(px,precorb(:,3),polyDegree));
a = coef_y1(3);
b = coef_y1(2);
c = coef_y1(1);
coef_y = [c/(g^2) b/g-(2*c*f)/(g^2) a-b*f/g+c*(f/g)^2];

coef_z1 = (polyfit(px,precorb(:,4),polyDegree));
a = coef_z1(3);
b = coef_z1(2);
c = coef_z1(1);
coef_z = [c/(g^2) b/g-(2*c*f)/(g^2) a-b*f/g+c*(f/g)^2];

vel_x = polyval(polyder(coef_x),precorb(:,1));
vel_y = polyval(polyder(coef_y),precorb(:,1));
vel_z = polyval(polyder(coef_z),precorb(:,1));
state_vector_acceleration = [vel_x vel_y vel_z];
coeff.state_vector_acceleration =state_vector_acceleration;

