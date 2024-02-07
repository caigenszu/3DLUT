
% Open the text file

function [orbit_position] = allsateorbit(slcpar)
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
start_time_line = lines{contains(lines, 'start_time')};
start_time = sscanf(start_time_line, 'start_time: %f');



% Extract number of state vectors
num_state_vectors_line = lines{contains(lines, 'number_of_state_vectors')};
num_state_vectors = sscanf(num_state_vectors_line, 'number_of_state_vectors: %d');

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
    state_vector_times(i) = first_state_vector_time + (i-1)*10; % Note: 10s is the state_vector_interval
end

orbit_position=[state_vector_times state_vector_positions state_vector_velocity];
