clear all; close all; clc;
d1 = load('data.csv');
% d2 = load('data_low_frictionloss.csv');

control_hz = 125;
plot_start_time = 2;
plot_end_time = 13;
plot_indices = control_hz*plot_start_time:control_hz*plot_end_time;

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
heading = d1(:, start_index);

data_len = size(d1, 1); % number of rows in d1
time = zeros(data_len);
for i=2:size(d1,1)-1
    time(i) = time(i-1) + inference_dt(i);

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
% 
% figure(2)
% sgtitle('Morph params')
% start_index = end_index+1;
% end_index = start_index+4;
% morph_params = d1(:, start_index:end_index);
% for i = 1:4
%     subplot(1,4,i);
%     smoothed_data = smooth(morph_params(:, i), 0.05, 'moving');
%     plot(time(plot_indices), smoothed_data(plot_indices, :));
% end

figure(3)
sgtitle('Velocity and Heading');
subplot(1,3,1);
lin_vel_data = q_dot_virtual(:, 1);
plot(time(plot_indices), lin_vel_data(plot_indices));
title('Linear velocity')
subplot(1,3,2);
ang_vel_data = q_dot_virtual(:, 2);
plot(time(plot_indices), ang_vel_data(plot_indices));
title('Angular velocity');
subplot(1,3,3);
plot(time(plot_indices), heading(plot_indices));
title('Heading')

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
