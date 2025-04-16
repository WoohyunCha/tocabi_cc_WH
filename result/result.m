close all; clc;
d1 = load('data.csv');

plot_start_time = 2;
plot_end_time = 12;

torque_data_yongarry = d2(:, 13:20);
plot_indices_yongarry = 250*plot_start_time:250*plot_end_time;

torque_data_yongarry_standstill = [d3(:, 10:21); d3(:, 10:21)];
grf_data_yongarry_standstill = [d3(:, 8:9);d3(:, 8:9)];
control_hz = 125;

plot_indices = control_hz*plot_start_time:control_hz*plot_end_time;


ylimit = [0, 0.7];


figure(1);
inference_dt_len = 1;
start_index = 1;
end_index = start_index + inference_dt_len - 1;
inference_dt = d1(:,start_index:end_index);
start_index = end_index + 1;

foot_force_len = 6;
end_index = start_index + foot_force_len - 1;
Lfoot_force_global = d1(:, start_index:end_index);
start_index = end_index + 1;
end_index = start_index + foot_force_len - 1;
Rfoot_force_global = d1(:, start_index:end_index);
start_index = end_index + 1;

command_torque_len = 33;
end_index = start_index + command_torque_len - 1;
command_torque = d1(:, start_index:end_index);
start_index = end_index + 1;
q_noise_len = 33;
end_index = start_index + q_noise_len - 1;
q_noise_ = d1(:, start_index:end_index);
start_index = end_index + 1;
q_dot_len = 33;
end_index = start_index + q_dot_len - 1;
q_dot_lpf_ = d1(:, start_index:end_index);
start_index = end_index + 1;
q_dot_virtual_len = 39;
end_index = start_index + q_dot_virtual_len - 1;
q_dot_virtual = d1(:, start_index:end_index);

start_index = end_index + 1;
q_virtual_len = 40;
end_index = start_index + q_virtual_len - 1;
q_virtual = d1(:, start_index:end_index);

start_index = end_index+1;
heading_len = 1;
end_index = start_index + heading_len - 1;
heading = d1(:, start_index:end_index);

start_index = end_index+1;
value_len = 1;
end_index = start_index + value_len - 1;
value_data = d1(:, start_index:end_index);

start_index = end_index+1;
flag_len = 1;
end_index = start_index + flag_len - 1;
flag_data = d1(:, start_index:end_index);

start_index = end_index+1;
command_len = 3;
end_index = start_index + command_len - 1;
command_data = d1(:, start_index:end_index);

data_len = size(d1, 1); % number of rows in d1
time = zeros(data_len);
for i=2:size(d1,1)-1
    time(i) = time(i-1) + inference_dt(i);
end
time_yongarry = zeros(size(d2, 1));
for i=2:size(d2,1)-1
    time_yongarry(i) = time_yongarry(i-1) + 1/250;
end

figure(1)
sgtitle("Foot Force Measurements in Global Frame")
for i = 1:foot_force_len
    subplot(2,3,i);
    plot(time(plot_indices), Lfoot_force_global(plot_indices, i));
    hold on;
    plot(time(plot_indices), Rfoot_force_global(plot_indices, i));
    legend("left foot", "right foot");
    if i == 1
        title("X direction")
    end
    if i == 2
        title("Y direction")
    end
    if i == 3
        title("Z direction")
    end
    if i == 4
        title("X moment")
    end
    if i == 5
        title("Y moment")
    end
    if i == 6
        title("Z moment")
    end
end


figure(7)
sgtitle('joint pose')
for i = 1:12
    subplot(2,6,i);
    plot(time(plot_indices), q_noise_(plot_indices, i));
end

figure(3)
sgtitle('Velocity and Heading');
subplot(2,2,1);
lin_vel_data = q_dot_virtual(:, 1);
plot(time(plot_indices), lin_vel_data(plot_indices));
hold on;
plot(time(plot_indices), command_data(plot_indices, 1));
hold on;
plot(time(plot_indices), mean(lin_vel_data(plot_indices)) * ones(size(plot_indices)), '--', 'LineWidth', 1.5) % Plot the mean of lin_vel_data
legend('Measured', 'Target', 'Mean');
title('Linear velocity')
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
grid on;

subplot(2,2,2);

lin_vel_data = q_dot_virtual(:, 2);  % Your input signal
alpha = 0.01;  % Adjust alpha (smoothing factor) based on desired filtering strength
lin_vel_data_lpf = causal_lpf(lin_vel_data, alpha);

plot(time(plot_indices),lin_vel_data(plot_indices));

title('Y drift velocity');
subplot(2,2,3);
ang_vel_data = q_dot_virtual(:, 6);
plot(time(plot_indices), ang_vel_data(plot_indices));
hold on;
plot(time(plot_indices), command_data(plot_indices, 2));
title('Angular velocity')
hold on;
plot(time(plot_indices), mean(ang_vel_data(plot_indices)) * ones(size(plot_indices)), '--', 'LineWidth', 1.5) % Plot the mean of lin_vel_data
legend('Measured', 'Target', 'Mean');
subplot(2,2,4);
plot(time(plot_indices), heading(plot_indices));
hold on;
plot(time(plot_indices), command_data(plot_indices, 3));
legend('Measured', 'Target');
title('Heading')

figure(4)
plot(time(plot_indices), value_data(plot_indices));
title('Value function');

figure(5)
sgtitle("Command joint torques");
for i=1:12
    subplot(2,6,i);
    plot(time(plot_indices), command_torque(plot_indices, i));
    hold on;
    plot(time_yongarry(plot_indices_yongarry), torque_data_yongarry_standstill(plot_indices_yongarry, i));
    legend('WH', 'yongarry');
end

figure(6)
sgtitle ("Joint velocities");
for i=1:12
    subplot(2,6,i);
    plot(time(plot_indices), q_dot_lpf_(plot_indices, i));
end

% figure();
% for i=1:12
%     subplot(2,6,i);
%     plot(d1(:,i));
%     hold on
%     plot(d2(:,i));
% end
% 
% figure();
% for i=13:33
%     subplot(4,6,i-12);
%     plot(d1(:,i));
%     hold on
%     plot(d2(:,i));
% end
% 
% %%
% clear all
% d3 = load('data.csv');
% 
% figure();
% for i=1:33
%     subplot(6,6,i);
%     plot(d3(:,1),d3(:,7+i))
% end
% 
% 
% figure();
% for i=1:33
%     subplot(6,6,i);
%     plot(d3(:,1),d3(:,73+i))
%     hold on
% plot(d3(:,1),d3(:,106+i))
% end

function lin_vel_data_lpf = causal_lpf(lin_vel_data, alpha)
    % Initialize the filtered data array
    lin_vel_data_lpf = zeros(size(lin_vel_data));
    
    % Set the first element equal to the input signal
    lin_vel_data_lpf(1) = lin_vel_data(1);
    
    % Apply the causal low-pass filter
    for t = 2:length(lin_vel_data)
        lin_vel_data_lpf(t) = alpha * lin_vel_data(t) + (1 - alpha) * lin_vel_data_lpf(t-1);
    end
end


figure(2)
sgtitle("Compare GRF with yongarry");
subplot(1,2,1);
plot(time(plot_indices), -Lfoot_force_global(plot_indices, 3));
hold on;
plot(time_yongarry(plot_indices_yongarry), grf_data_yongarry_standstill(plot_indices_yongarry, 1));
legend("WH", "yongarry");
title("Left foot Z direction GRF");
subplot(1,2,2);
plot(time(plot_indices), -Rfoot_force_global(plot_indices, 3));
hold on;
plot(time_yongarry(plot_indices_yongarry), grf_data_yongarry_standstill(plot_indices_yongarry, 2));
legend("WH", "yongarry");
title("Right foot Z direction GRF");