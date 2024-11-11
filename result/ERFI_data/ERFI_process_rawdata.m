close all; clc;

%% Load Data. Do only once
% clear all;
% dr_act_4_data = load('DR/actuator_gap/0.4ms/data.csv');
% dr_act_6_data = load('DR/actuator_gap/0.6ms/data.csv');
% dr_env_4_data = load('DR/environment_gap/0.4ms/data.csv');
% dr_env_6_data = load('DR/environment_gap/0.6ms/data.csv');
% dr_mod_4_data = load('DR/model_gap/0.4ms/data.csv');
% dr_mod_6_data = load('DR/model_gap/0.6ms/data.csv');
% dr_nom_4_data = load('DR/nominal/0.4ms/data.csv');
% dr_nom_6_data = load('DR/nominal/0.6ms/data.csv');
% 
% %% Load noise lim 50 data
% erfi_50_act_4_data = load('ERFI_50_uniform/actuator_gap/0.4/data.csv');
% erfi_50_act_6_data = load('ERFI_50_uniform/actuator_gap/0.6/data.csv');
% erfi_50_env_4_data = load('ERFI_50_uniform/environment_gap/0.4/data.csv');
% erfi_50_env_6_data = load('ERFI_50_uniform/environment_gap/0.6/data.csv');
% erfi_50_mod_4_data = load('ERFI_50_uniform/model_gap/0.4/data.csv');
% erfi_50_mod_6_data = load('ERFI_50_uniform/model_gap/0.6/data.csv');
% erfi_50_nom_4_data = load('ERFI_50_uniform/nominal/0.4/data.csv');
% erfi_50_nom_6_data = load('ERFI_50_uniform/nominal/0.6/data.csv');
% 
% nn_50_act_4_data = load('NN_exp_50_gaussian_1.5/actuator_gap/0.4/data.csv');
% nn_50_act_6_data = load('NN_exp_50_gaussian_1.5/actuator_gap/0.6/data.csv');
% nn_50_env_4_data = load('NN_exp_50_gaussian_1.5/environment_gap/0.4/data.csv');
% nn_50_env_6_data = load('NN_exp_50_gaussian_1.5/environment_gap/0.6/data.csv');
% nn_50_mod_4_data = load('NN_exp_50_gaussian_1.5/model_gap/0.4/data.csv');
% nn_50_mod_6_data = load('NN_exp_50_gaussian_1.5/model_gap/0.6/data.csv');
% nn_50_nom_4_data = load('NN_exp_50_gaussian_1.5/nominal/0.4/data.csv');
% nn_50_nom_6_data = load('NN_exp_50_gaussian_1.5/nominal/0.6/data.csv');
% 
% nq_50_act_4_data = load('NN_nq_50_gaussian_1.5/actuator_gap/0.4/data.csv');
% nq_50_act_6_data = load('NN_nq_50_gaussian_1.5/actuator_gap/0.6/data.csv');
% nq_50_env_4_data = load('NN_nq_50_gaussian_1.5/environment_gap/0.4/data.csv');
% nq_50_env_6_data = load('NN_nq_50_gaussian_1.5/environment_gap/0.6/data.csv');
% nq_50_mod_4_data = load('NN_nq_50_gaussian_1.5/model_gap/0.4/data.csv');
% nq_50_mod_6_data = load('NN_nq_50_gaussian_1.5/model_gap/0.6/data.csv');
% nq_50_nom_4_data = load('NN_nq_50_gaussian_1.5/nominal/0.4/data.csv');
% nq_50_nom_6_data = load('NN_nq_50_gaussian_1.5/nominal/0.6/data.csv');
% 
% q_50_act_4_data = load('NN_50_gaussian_1.5/actuator_gap/0.4/data.csv');
% q_50_act_6_data = load('NN_50_gaussian_1.5/actuator_gap/0.6/data.csv');
% q_50_env_4_data = load('NN_50_gaussian_1.5/environment_gap/0.4/data.csv');
% q_50_env_6_data = load('NN_50_gaussian_1.5/environment_gap/0.6/data.csv');
% q_50_mod_4_data = load('NN_50_gaussian_1.5/model_gap/0.4/data.csv');
% q_50_mod_6_data = load('NN_50_gaussian_1.5/model_gap/0.6/data.csv');
% q_50_nom_4_data = load('NN_50_gaussian_1.5/nominal/0.4/data.csv');
% q_50_nom_6_data = load('NN_50_gaussian_1.5/nominal/0.6/data.csv');
% 
% %% Load noise lim 30 data
% % erfi_30_act_4_data = load('ERFI_30_uniform/actuator_gap/0.4/data.csv');
% % erfi_30_act_6_data = load('ERFI_30_uniform/actuator_gap/0.6/data.csv');
% % erfi_30_env_4_data = load('ERFI_30_uniform/environment_gap/0.4/data.csv');
% % erfi_30_env_6_data = load('ERFI_30_uniform/environment_gap/0.6/data.csv');
% % erfi_30_mod_4_data = load('ERFI_30_uniform/model_gap/0.4/data.csv');
% % erfi_30_mod_6_data = load('ERFI_30_uniform/model_gap/0.6/data.csv');
% % erfi_30_nom_4_data = load('ERFI_30_uniform/nominal/0.4/data.csv');
% % erfi_30_nom_6_data = load('ERFI_30_uniform/nominal/0.6/data.csv');
% % 
% % nn_30_act_4_data = load('NN_exp_30_gaussian_1.5/actuator_gap/0.4/data.csv');
% % nn_30_act_6_data = load('NN_exp_30_gaussian_1.5/actuator_gap/0.6/data.csv');
% % nn_30_env_4_data = load('NN_exp_30_gaussian_1.5/environment_gap/0.4/data.csv');
% % nn_30_env_6_data = load('NN_exp_30_gaussian_1.5/environment_gap/0.6/data.csv');
% % nn_30_mod_4_data = load('NN_exp_30_gaussian_1.5/model_gap/0.4/data.csv');
% % nn_30_mod_6_data = load('NN_exp_30_gaussian_1.5/model_gap/0.6/data.csv');
% % nn_30_nom_4_data = load('NN_exp_30_gaussian_1.5/nominal/0.4/data.csv');
% % nn_30_nom_6_data = load('NN_exp_30_gaussian_1.5/nominal/0.6/data.csv');
% % 
% % nq_30_act_4_data = load('NN_nq_30_gaussian_1.5/actuator_gap/0.4/data.csv');
% % nq_30_act_6_data = load('NN_nq_30_gaussian_1.5/actuator_gap/0.6/data.csv');
% % nq_30_env_4_data = load('NN_nq_30_gaussian_1.5/environment_gap/0.4/data.csv');
% % nq_30_env_6_data = load('NN_nq_30_gaussian_1.5/environment_gap/0.6/data.csv');
% % nq_30_mod_4_data = load('NN_nq_30_gaussian_1.5/model_gap/0.4/data.csv');
% % nq_30_mod_6_data = load('NN_nq_30_gaussian_1.5/model_gap/0.6/data.csv');
% % nq_30_nom_4_data = load('NN_nq_30_gaussian_1.5/nominal/0.4/data.csv');
% % nq_30_nom_6_data = load('NN_nq_30_gaussian_1.5/nominal/0.6/data.csv');
% % 
% % q_30_act_4_data = load('NN_30_gaussian_1.5/actuator_gap/0.4/data.csv');
% % q_30_act_6_data = load('NN_30_gaussian_1.5/actuator_gap/0.6/data.csv');
% % q_30_env_4_data = load('NN_30_gaussian_1.5/environment_gap/0.4/data.csv');
% % q_30_env_6_data = load('NN_30_gaussian_1.5/environment_gap/0.6/data.csv');
% % q_30_mod_4_data = load('NN_30_gaussian_1.5/model_gap/0.4/data.csv');
% % q_30_mod_6_data = load('NN_30_gaussian_1.5/model_gap/0.6/data.csv');
% % q_30_nom_4_data = load('NN_30_gaussian_1.5/nominal/0.4/data.csv');
% % q_30_nom_6_data = load('NN_30_gaussian_1.5/nominal/0.6/data.csv');
% 
% 
% %% Load noise lim 20 data
% erfi_20_act_4_data = load('ERFI_20_uniform/actuator_gap/0.4/data.csv');
% erfi_20_act_6_data = load('ERFI_20_uniform/actuator_gap/0.6/data.csv');
% erfi_20_env_4_data = load('ERFI_20_uniform/environment_gap/0.4/data.csv');
% erfi_20_env_6_data = load('ERFI_20_uniform/environment_gap/0.6/data.csv');
% erfi_20_mod_4_data = load('ERFI_20_uniform/model_gap/0.4/data.csv');
% erfi_20_mod_6_data = load('ERFI_20_uniform/model_gap/0.6/data.csv');
% erfi_20_nom_4_data = load('ERFI_20_uniform/nominal/0.4/data.csv');
% erfi_20_nom_6_data = load('ERFI_20_uniform/nominal/0.6/data.csv');
% 
% nn_20_act_4_data = load('NN_exp_20_gaussian_1.5/actuator_gap/0.4/data.csv');
% nn_20_act_6_data = load('NN_exp_20_gaussian_1.5/actuator_gap/0.6/data.csv');
% nn_20_env_4_data = load('NN_exp_20_gaussian_1.5/environment_gap/0.4/data.csv');
% nn_20_env_6_data = load('NN_exp_20_gaussian_1.5/environment_gap/0.6/data.csv');
% nn_20_mod_4_data = load('NN_exp_20_gaussian_1.5/model_gap/0.4/data.csv');
% nn_20_mod_6_data = load('NN_exp_20_gaussian_1.5/model_gap/0.6/data.csv');
% nn_20_nom_4_data = load('NN_exp_20_gaussian_1.5/nominal/0.4/data.csv');
% nn_20_nom_6_data = load('NN_exp_20_gaussian_1.5/nominal/0.6/data.csv');
% 
% nq_20_act_4_data = load('NN_nq_20_gaussian_1.5/actuator_gap/0.4/data.csv');
% nq_20_act_6_data = load('NN_nq_20_gaussian_1.5/actuator_gap/0.6/data.csv');
% nq_20_env_4_data = load('NN_nq_20_gaussian_1.5/environment_gap/0.4/data.csv');
% nq_20_env_6_data = load('NN_nq_20_gaussian_1.5/environment_gap/0.6/data.csv');
% nq_20_mod_4_data = load('NN_nq_20_gaussian_1.5/model_gap/0.4/data.csv');
% nq_20_mod_6_data = load('NN_nq_20_gaussian_1.5/model_gap/0.6/data.csv');
% nq_20_nom_4_data = load('NN_nq_20_gaussian_1.5/nominal/0.4/data.csv');
% nq_20_nom_6_data = load('NN_nq_20_gaussian_1.5/nominal/0.6/data.csv');
% 
% q_20_act_4_data = load('NN_20_gaussian_1.5/actuator_gap/0.4/data.csv');
% q_20_act_6_data = load('NN_20_gaussian_1.5/actuator_gap/0.6/data.csv');
% q_20_env_4_data = load('NN_20_gaussian_1.5/environment_gap/0.4/data.csv');
% q_20_env_6_data = load('NN_20_gaussian_1.5/environment_gap/0.6/data.csv');
% q_20_mod_4_data = load('NN_20_gaussian_1.5/model_gap/0.4/data.csv');
% q_20_mod_6_data = load('NN_20_gaussian_1.5/model_gap/0.6/data.csv');
% q_20_nom_4_data = load('NN_20_gaussian_1.5/nominal/0.4/data.csv');
% q_20_nom_6_data = load('NN_20_gaussian_1.5/nominal/0.6/data.csv');

% Load and process each data file


% Helper function to parse data with the same structure
% Load and process each data file



%% Helper function to parse data with the same structure and save to workspace
function process_data(file_data, prefix)
    % Define the lengths of each data segment
    inference_dt_len = 1;
    foot_force_len = 6;
    command_torque_len = 33;
    q_noise_len = 33;
    q_dot_len = 33;
    q_dot_virtual_len = 39;
    q_virtual_len = 40;
    heading_len = 1;
    value_len = 1;
    flag_len = 1;
    command_len = 3;
    latent_len = 512;

    start_index = 1;
    end_index = start_index + inference_dt_len - 1;
    assignin('base', [prefix, '_inference_dt'], file_data(:, start_index:end_index));
    assignin('base', [prefix, '_time'], cumsum(file_data(:, start_index:end_index)));
    start_index = end_index + 1;
    
    end_index = start_index + foot_force_len - 1;
    assignin('base', [prefix, '_Lfoot_force_global'], file_data(:, start_index:end_index));
    start_index = end_index + 1;
    end_index = start_index + foot_force_len - 1;
    assignin('base', [prefix, '_Rfoot_force_global'], file_data(:, start_index:end_index));
    start_index = end_index + 1;
    
    end_index = start_index + command_torque_len - 1;
    assignin('base', [prefix, '_command_torque'], file_data(:, start_index:end_index));
    start_index = end_index + 1;
    
    end_index = start_index + q_noise_len - 1;
    assignin('base', [prefix, '_q_noise'], file_data(:, start_index:end_index));
    start_index = end_index + 1;
    
    end_index = start_index + q_dot_len - 1;
    assignin('base', [prefix, '_q_dot_lpf'], file_data(:, start_index:end_index));
    start_index = end_index + 1;
    
    end_index = start_index + q_dot_virtual_len - 1;
    assignin('base', [prefix, '_q_dot_virtual'], file_data(:, start_index:end_index));
    start_index = end_index + 1;
    
    end_index = start_index + q_virtual_len - 1;
    assignin('base', [prefix, '_q_virtual'], file_data(:, start_index:end_index));
    start_index = end_index + 1;
    
    end_index = start_index + heading_len - 1;
    assignin('base', [prefix, '_heading'], file_data(:, start_index:end_index));
    start_index = end_index + 1;
    
    end_index = start_index + value_len - 1;
    assignin('base', [prefix, '_value_data'], file_data(:, start_index:end_index));
    start_index = end_index + 1;
    
    end_index = start_index + flag_len - 1;
    assignin('base', [prefix, '_flag_data'], file_data(:, start_index:end_index));
    start_index = end_index + 1;
    
    end_index = start_index + command_len - 1;
    assignin('base', [prefix, '_command_data'], file_data(:, start_index:end_index));
    start_index = end_index + 1;
    
    end_index = start_index + latent_len - 1;
    assignin('base', [prefix, '_latent_data'], file_data(:, start_index:end_index));

end
% 
% % Helper function to plot full data
% function plot_full_data(data1, data2, data3, )
% 
% 
% end

%% Process each data file
process_data(dr_act_4_data, 'dr_act_4');
process_data(dr_act_6_data, 'dr_act_6');
process_data(dr_env_4_data, 'dr_env_4');
process_data(dr_env_6_data, 'dr_env_6');
process_data(dr_mod_4_data, 'dr_mod_4');
process_data(dr_mod_6_data, 'dr_mod_6');
process_data(dr_nom_4_data, 'dr_nom_4');
process_data(dr_nom_6_data, 'dr_nom_6');

process_data(erfi_50_act_4_data, 'erfi_50_act_4');
process_data(erfi_50_act_6_data, 'erfi_50_act_6');
process_data(erfi_50_env_4_data, 'erfi_50_env_4');
process_data(erfi_50_env_6_data, 'erfi_50_env_6');
process_data(erfi_50_mod_4_data, 'erfi_50_mod_4');
process_data(erfi_50_mod_6_data, 'erfi_50_mod_6');
process_data(erfi_50_nom_4_data, 'erfi_50_nom_4');
process_data(erfi_50_nom_6_data, 'erfi_50_nom_6');

process_data(nn_50_act_4_data, 'nn_50_act_4');
process_data(nn_50_act_6_data, 'nn_50_act_6');
process_data(nn_50_env_4_data, 'nn_50_env_4');
process_data(nn_50_env_6_data, 'nn_50_env_6');
process_data(nn_50_mod_4_data, 'nn_50_mod_4');
process_data(nn_50_mod_6_data, 'nn_50_mod_6');
process_data(nn_50_nom_4_data, 'nn_50_nom_4');
process_data(nn_50_nom_6_data, 'nn_50_nom_6');

process_data(nq_50_act_4_data, 'nq_50_act_4');
process_data(nq_50_act_6_data, 'nq_50_act_6');
process_data(nq_50_env_4_data, 'nq_50_env_4');
process_data(nq_50_env_6_data, 'nq_50_env_6');
process_data(nq_50_mod_4_data, 'nq_50_mod_4');
process_data(nq_50_mod_6_data, 'nq_50_mod_6');
process_data(nq_50_nom_4_data, 'nq_50_nom_4');
process_data(nq_50_nom_6_data, 'nq_50_nom_6');

process_data(q_50_act_4_data, 'q_50_act_4');
process_data(q_50_act_6_data, 'q_50_act_6');
process_data(q_50_env_4_data, 'q_50_env_4');
process_data(q_50_env_6_data, 'q_50_env_6');
process_data(q_50_mod_4_data, 'q_50_mod_4');
process_data(q_50_mod_6_data, 'q_50_mod_6');
process_data(q_50_nom_4_data, 'q_50_nom_4');
process_data(q_50_nom_6_data, 'q_50_nom_6');


process_data(erfi_30_act_4_data, 'erfi_30_act_4');
process_data(erfi_30_act_6_data, 'erfi_30_act_6');
process_data(erfi_30_env_4_data, 'erfi_30_env_4');
process_data(erfi_30_env_6_data, 'erfi_30_env_6');
process_data(erfi_30_mod_4_data, 'erfi_30_mod_4');
process_data(erfi_30_mod_6_data, 'erfi_30_mod_6');
process_data(erfi_30_nom_4_data, 'erfi_30_nom_4');
process_data(erfi_30_nom_6_data, 'erfi_30_nom_6');
% 
% process_data(nn_30_act_4_data, 'nn_30_act_4');
% process_data(nn_30_act_6_data, 'nn_30_act_6');
% process_data(nn_30_env_4_data, 'nn_30_env_4');
% process_data(nn_30_env_6_data, 'nn_30_env_6');
% process_data(nn_30_mod_4_data, 'nn_30_mod_4');
% process_data(nn_30_mod_6_data, 'nn_30_mod_6');
% process_data(nn_30_nom_4_data, 'nn_30_nom_4');
% process_data(nn_30_nom_6_data, 'nn_30_nom_6');
% 
% process_data(nq_30_act_4_data, 'nq_30_act_4');
% process_data(nq_30_act_6_data, 'nq_30_act_6');
% process_data(nq_30_env_4_data, 'nq_30_env_4');
% process_data(nq_30_env_6_data, 'nq_30_env_6');
% process_data(nq_30_mod_4_data, 'nq_30_mod_4');
% process_data(nq_30_mod_6_data, 'nq_30_mod_6');
% process_data(nq_30_nom_4_data, 'nq_30_nom_4');
% process_data(nq_30_nom_6_data, 'nq_30_nom_6');
% 
% process_data(q_30_act_4_data, 'q_30_act_4');
% process_data(q_30_act_6_data, 'q_30_act_6');
% process_data(q_30_env_4_data, 'q_30_env_4');
% process_data(q_30_env_6_data, 'q_30_env_6');
% process_data(q_30_mod_4_data, 'q_30_mod_4');
% process_data(q_30_mod_6_data, 'q_30_mod_6');
% process_data(q_30_nom_4_data, 'q_30_nom_4');
% process_data(q_30_nom_6_data, 'q_30_nom_6');


process_data(erfi_20_act_4_data, 'erfi_20_act_4');
process_data(erfi_20_act_6_data, 'erfi_20_act_6');
process_data(erfi_20_env_4_data, 'erfi_20_env_4');
process_data(erfi_20_env_6_data, 'erfi_20_env_6');
process_data(erfi_20_mod_4_data, 'erfi_20_mod_4');
process_data(erfi_20_mod_6_data, 'erfi_20_mod_6');
process_data(erfi_20_nom_4_data, 'erfi_20_nom_4');
process_data(erfi_20_nom_6_data, 'erfi_20_nom_6');

process_data(nn_20_act_4_data, 'nn_20_act_4');
process_data(nn_20_act_6_data, 'nn_20_act_6');
process_data(nn_20_env_4_data, 'nn_20_env_4');
process_data(nn_20_env_6_data, 'nn_20_env_6');
process_data(nn_20_mod_4_data, 'nn_20_mod_4');
process_data(nn_20_mod_6_data, 'nn_20_mod_6');
process_data(nn_20_nom_4_data, 'nn_20_nom_4');
process_data(nn_20_nom_6_data, 'nn_20_nom_6');

process_data(nq_20_act_4_data, 'nq_20_act_4');
process_data(nq_20_act_6_data, 'nq_20_act_6');
process_data(nq_20_env_4_data, 'nq_20_env_4');
process_data(nq_20_env_6_data, 'nq_20_env_6');
process_data(nq_20_mod_4_data, 'nq_20_mod_4');
process_data(nq_20_mod_6_data, 'nq_20_mod_6');
process_data(nq_20_nom_4_data, 'nq_20_nom_4');
process_data(nq_20_nom_6_data, 'nq_20_nom_6');

process_data(q_20_act_4_data, 'q_20_act_4');
process_data(q_20_act_6_data, 'q_20_act_6');
process_data(q_20_env_4_data, 'q_20_env_4');
process_data(q_20_env_6_data, 'q_20_env_6');
process_data(q_20_mod_4_data, 'q_20_mod_4');
process_data(q_20_mod_6_data, 'q_20_mod_6');
process_data(q_20_nom_4_data, 'q_20_nom_4');
process_data(q_20_nom_6_data, 'q_20_nom_6');



plot_start_time = 0.008;
plot_end_time = 10;
control_hz = 125;
plot_indices = control_hz*plot_start_time:control_hz*plot_end_time;
% Define figure size for single or double column
single_column_size = [4, 1.5]; % [width, height] in inches for single column
double_column_size = [7, 3]; % [width, height] in inches for double column
double_column_single_row_size = [7,1.5];
double_column_eight_row_size = [7, 12];
ylimit = [0, 0.7];
sgtitle_fontsize = 10;
subtitle_fontsize = 6;
label_fontsize = 4;
legend_fontsize = 4;
tick_fontsize = 4;

%% Compare Performance in Nominal Setting 
figure('Name', 'Nominal performance comparison, (DR, ERFI, Ours)');
set(gcf, 'Units', 'inches', 'Position', [1, 1, double_column_size]); % Choose single_column_size or double_column_size here
% sgtitle("Base Forward Velocity Tracking Performance in Nominal Settings", 'FontSize', sgtitle_fontsize);
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'normal'); % Adjust 'TileSpacing' and 'Padding' as needed
nexttile;
plot(dr_nom_4_time(plot_indices), dr_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(dr_nom_4_time(plot_indices), dr_nom_4_command_data(plot_indices, 1));
title('DR', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');

nexttile;
plot(erfi_50_nom_4_time(plot_indices), erfi_50_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(erfi_50_nom_4_time(plot_indices), erfi_50_nom_4_command_data(plot_indices, 1));
title('ERFI', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast')

nexttile;
plot(nn_50_nom_4_time(plot_indices), nn_50_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(nn_50_nom_4_time(plot_indices), nn_50_nom_4_command_data(plot_indices, 1));
title('Our method', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
lgd.Position(1) = lgd.Position(1) + 0.02;
lgd.Position(2) = lgd.Position(2) - 0.01;

% sgtitle("Base Forward Velocity Tracking Performance in Nominal Settings", 'FontSize', sgtitle_fontsize);
nexttile;
plot(dr_nom_6_time(plot_indices), dr_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(dr_nom_6_time(plot_indices), dr_nom_6_command_data(plot_indices, 1));
% title('DR', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'southeast')

nexttile;
plot(erfi_50_nom_6_time(plot_indices), erfi_50_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(erfi_50_nom_6_time(plot_indices), erfi_50_nom_6_command_data(plot_indices, 1));
% title('ERFI', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'southeast')

nexttile;
plot(nn_50_nom_6_time(plot_indices), nn_50_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(nn_50_nom_6_time(plot_indices), nn_50_nom_6_command_data(plot_indices, 1));
% title('NN', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'southeast')

han = axes(gcf, 'Visible', 'off'); % Create an invisible axes that spans the figure
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
xlabel(han, 'Time, [s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [0.5, -0.05, 0]); % Global x-axis label
ylabel(han, 'Velocity, [m/s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [-0.05, 0.5, 0]);


%% Compare DR, ERFI, and our method under actuator gap
figure('Name', 'Actuator gap, (DR, ERFI, Ours)');
set(gcf, 'Units', 'inches', 'Position', [1, 1, double_column_size]); % Choose single_column_size or double_column_size here

t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'normal'); % Adjust 'TileSpacing' and 'Padding' as needed
nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(dr_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    title('DR', 'FontSize', subtitle_fontsize);
else
    plot(dr_act_4_time(plot_indices), dr_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(dr_act_4_time(plot_indices), dr_act_4_command_data(plot_indices, 1));
    title('DR', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(erfi_50_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    title('ERFI', 'FontSize', subtitle_fontsize);
else
    plot(erfi_50_act_4_time(plot_indices), erfi_50_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_50_act_4_time(plot_indices), erfi_50_act_4_command_data(plot_indices, 1));
    title('ERFI', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(nn_50_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    title('Our method', 'FontSize', subtitle_fontsize);
else
    plot(nn_50_act_4_time(plot_indices), nn_50_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nn_50_act_4_time(plot_indices), nn_50_act_4_command_data(plot_indices, 1));
    title('Our method', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    lgd.Position(1) = lgd.Position(1) + 0.02;
    lgd.Position(2) = lgd.Position(2) - 0.01;
end

% 0.6m/s
nexttile;
if size(dr_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('DR', 'FontSize', subtitle_fontsize);
else
    plot(dr_act_6_time(plot_indices), dr_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(dr_act_6_time(plot_indices), dr_act_6_command_data(plot_indices, 1));
    % title('DR', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(erfi_50_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('ERFI', 'FontSize', subtitle_fontsize);
else
    plot(erfi_50_act_6_time(plot_indices), erfi_50_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_50_act_6_time(plot_indices), erfi_50_act_6_command_data(plot_indices, 1));
    % title('ERFI', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(nn_50_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('Our method', 'FontSize', subtitle_fontsize);
else
    plot(nn_50_act_6_time(plot_indices), nn_50_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nn_50_act_6_time(plot_indices), nn_50_act_6_command_data(plot_indices, 1));
    % title('Our method', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.02;
    % lg.Position(2) = lg.Position(2) - 0.06;
end

han = axes(gcf, 'Visible', 'off'); % Create an invisible axes that spans the figure
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
xlabel(han, 'Time, [s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [0.5, -0.05, 0]); % Global x-axis label
ylabel(han, 'Velocity, [m/s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [-0.07, 0.5, 0]);


%% Compare DR, ERFI, and our method under model gap%%%
figure('Name', 'Model gap, (DR, ERFI, Ours)');
set(gcf, 'Units', 'inches', 'Position', [1, 1, double_column_size]); % Choose single_column_size or double_column_size here

t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'normal'); % Adjust 'TileSpacing' and 'Padding' as needed
nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(dr_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    title('DR', 'FontSize', subtitle_fontsize);
else
    plot(dr_mod_4_time(plot_indices), dr_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(dr_mod_4_time(plot_indices), dr_mod_4_command_data(plot_indices, 1));
    title('DR', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(erfi_50_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    title('ERFI', 'FontSize', subtitle_fontsize);
else
    plot(erfi_50_mod_4_time(plot_indices), erfi_50_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_50_mod_4_time(plot_indices), erfi_50_mod_4_command_data(plot_indices, 1));
    title('ERFI', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(nn_50_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    title('Our method', 'FontSize', subtitle_fontsize);
else
    plot(nn_50_mod_4_time(plot_indices), nn_50_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nn_50_mod_4_time(plot_indices), nn_50_mod_4_command_data(plot_indices, 1));
    title('Our method', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    lgd.Position(1) = lgd.Position(1) + 0.02;
    lgd.Position(2) = lgd.Position(2) - 0.01;
end

% 0.6m/s
nexttile;
if size(dr_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('DR', 'FontSize', subtitle_fontsize);
else
    plot(dr_mod_6_time(plot_indices), dr_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(dr_mod_6_time(plot_indices), dr_mod_6_command_data(plot_indices, 1));
    % title('DR', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(erfi_50_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('ERFI', 'FontSize', subtitle_fontsize);
else
    plot(erfi_50_mod_6_time(plot_indices), erfi_50_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_50_mod_6_time(plot_indices), erfi_50_mod_6_command_data(plot_indices, 1));
    % title('ERFI', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(nn_50_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('Our method', 'FontSize', subtitle_fontsize);
else
    plot(nn_50_mod_6_time(plot_indices), nn_50_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nn_50_mod_6_time(plot_indices), nn_50_mod_6_command_data(plot_indices, 1));
    % title('Our method', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.02;
    % lg.Position(2) = lg.Position(2) - 0.06;
end

han = axes(gcf, 'Visible', 'off'); % Create an invisible axes that spans the figure
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
xlabel(han, 'Time, [s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [0.5, -0.05, 0]); % Global x-axis label
ylabel(han, 'Velocity, [m/s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [-0.07, 0.5, 0]);



%% Compare DR, ERFI, and our method under environment gap%%%
figure('Name', 'Environment gap, (DR, ERFI, Ours)');
set(gcf, 'Units', 'inches', 'Position', [1, 1, double_column_size]); % Choose single_column_size or double_column_size here

t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'normal'); % Adjust 'TileSpacing' and 'Padding' as needed
nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(dr_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    title('DR', 'FontSize', subtitle_fontsize);
else
    plot(dr_env_4_time(plot_indices), dr_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(dr_env_4_time(plot_indices), dr_env_4_command_data(plot_indices, 1));
    title('DR', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(erfi_50_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    title('ERFI', 'FontSize', subtitle_fontsize);
else
    plot(erfi_50_env_4_time(plot_indices), erfi_50_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_50_env_4_time(plot_indices), erfi_50_env_4_command_data(plot_indices, 1));
    title('ERFI', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(nn_50_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    title('Our method', 'FontSize', subtitle_fontsize);
else
    plot(nn_50_env_4_time(plot_indices), nn_50_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nn_50_env_4_time(plot_indices), nn_50_env_4_command_data(plot_indices, 1));
    title('Our method', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    lgd.Position(1) = lgd.Position(1) + 0.02;
    lgd.Position(2) = lgd.Position(2) - 0.01;
end

% 0.6m/s
nexttile;
if size(dr_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('DR', 'FontSize', subtitle_fontsize);
else
    plot(dr_env_6_time(plot_indices), dr_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(dr_env_6_time(plot_indices), dr_env_6_command_data(plot_indices, 1));
    % title('DR', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(erfi_50_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('ERFI', 'FontSize', subtitle_fontsize);
else
    plot(erfi_50_env_6_time(plot_indices), erfi_50_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_50_env_6_time(plot_indices), erfi_50_env_6_command_data(plot_indices, 1));
    % title('ERFI', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(nn_50_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('Our method', 'FontSize', subtitle_fontsize);
else
    plot(nn_50_env_6_time(plot_indices), nn_50_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nn_50_env_6_time(plot_indices), nn_50_env_6_command_data(plot_indices, 1));
    % title('Our method', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.02;
    % lg.Position(2) = lg.Position(2) - 0.06;
end

han = axes(gcf, 'Visible', 'off'); % Create an invisible axes that spans the figure
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
xlabel(han, 'Time, [s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [0.5, -0.05, 0]); % Global x-axis label
ylabel(han, 'Velocity, [m/s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [-0.07, 0.5, 0]);


%% Evaluate NN methods per noise lim
%% Noise lim 50, NN methods
figure('Name', 'NN method comparison at lim 50Nm, (No quad, quad, explicit');
set(gcf, 'Units', 'inches', 'Position', [1, 1, double_column_eight_row_size]); % Choose single_column_size or double_column_size here
% sgtitle("Base Forward Velocity Tracking Performance in Nominal Settings", 'FontSize', sgtitle_fontsize);
t = tiledlayout(8, 3, 'TileSpacing', 'compact', 'Padding', 'normal'); % Adjust 'TileSpacing' and 'Padding' as needed
nexttile;
plot(nq_50_nom_4_time(plot_indices), nq_50_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(nq_50_nom_4_time(plot_indices), nq_50_nom_4_command_data(plot_indices, 1));
title('No quadratic term', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');

nexttile;
plot(q_50_nom_4_time(plot_indices), q_50_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(q_50_nom_4_time(plot_indices), q_50_nom_4_command_data(plot_indices, 1));
title('Quadratic term', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast')

nexttile;
plot(nn_50_nom_4_time(plot_indices), nn_50_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(nn_50_nom_4_time(plot_indices), nn_50_nom_4_command_data(plot_indices, 1));
title('Our method', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
lgd.Position(1) = lgd.Position(1) + 0.02;
lgd.Position(2) = lgd.Position(2) - 0.01;

% sgtitle("Base Forward Velocity Tracking Performance in Nominal Settings", 'FontSize', sgtitle_fontsize);
nexttile;
plot(nq_50_nom_6_time(plot_indices), nq_50_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(nq_50_nom_6_time(plot_indices), nq_50_nom_6_command_data(plot_indices, 1));
% title('No quadratic term', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'southeast')

nexttile;
plot(q_50_nom_6_time(plot_indices), q_50_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(q_50_nom_6_time(plot_indices), q_50_nom_6_command_data(plot_indices, 1));
% title('q_50', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'southeast')

nexttile;
plot(nn_50_nom_6_time(plot_indices), nn_50_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(nn_50_nom_6_time(plot_indices), nn_50_nom_6_command_data(plot_indices, 1));
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'southeast')


% Compare No quadratic term, q_50, and our method under actuator gap
nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(nq_50_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(nq_50_act_4_time(plot_indices), nq_50_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nq_50_act_4_time(plot_indices), nq_50_act_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(q_50_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(q_50_act_4_time(plot_indices), q_50_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(q_50_act_4_time(plot_indices), q_50_act_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(nn_50_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(nn_50_act_4_time(plot_indices), nn_50_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nn_50_act_4_time(plot_indices), nn_50_act_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
end

% 0.6m/s
nexttile;
if size(nq_50_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('DR', 'FontSize', subtitle_fontsize);
else
    plot(nq_50_act_6_time(plot_indices), nq_50_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nq_50_act_6_time(plot_indices), nq_50_act_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(q_50_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('q_50', 'FontSize', subtitle_fontsize);
else
    plot(q_50_act_6_time(plot_indices), q_50_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(q_50_act_6_time(plot_indices), q_50_act_6_command_data(plot_indices, 1));
    % title('q_50', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(nn_50_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('Our method', 'FontSize', subtitle_fontsize);
else
    plot(nn_50_act_6_time(plot_indices), nn_50_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nn_50_act_6_time(plot_indices), nn_50_act_6_command_data(plot_indices, 1));
    % title('Our method', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.02;
    % lg.Position(2) = lg.Position(2) - 0.06;
end


% Compare No quadratic term, q_50, and our method under model gap%%%

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(nq_50_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(nq_50_mod_4_time(plot_indices), nq_50_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nq_50_mod_4_time(plot_indices), nq_50_mod_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(q_50_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(q_50_mod_4_time(plot_indices), q_50_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(q_50_mod_4_time(plot_indices), q_50_mod_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(nn_50_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(nn_50_mod_4_time(plot_indices), nn_50_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nn_50_mod_4_time(plot_indices), nn_50_mod_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
end

% 0.6m/s
nexttile;
if size(nq_50_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(nq_50_mod_6_time(plot_indices), nq_50_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nq_50_mod_6_time(plot_indices), nq_50_mod_6_command_data(plot_indices, 1));
    % title('DR', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(q_50_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('q_50', 'FontSize', subtitle_fontsize);
else
    plot(q_50_mod_6_time(plot_indices), q_50_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(q_50_mod_6_time(plot_indices), q_50_mod_6_command_data(plot_indices, 1));
    % title('q_50', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(nn_50_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('Our method', 'FontSize', subtitle_fontsize);
else
    plot(nn_50_mod_6_time(plot_indices), nn_50_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nn_50_mod_6_time(plot_indices), nn_50_mod_6_command_data(plot_indices, 1));
    % title('Our method', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.02;
    % lg.Position(2) = lg.Position(2) - 0.06;
end




% Compare No quadratic term, q_50, and our method under environment gap%%%
nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(nq_50_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(nq_50_env_4_time(plot_indices), nq_50_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nq_50_env_4_time(plot_indices), nq_50_env_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(q_50_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(q_50_env_4_time(plot_indices), q_50_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(q_50_env_4_time(plot_indices), q_50_env_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(nn_50_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(nn_50_env_4_time(plot_indices), nn_50_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nn_50_env_4_time(plot_indices), nn_50_env_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
end

% 0.6m/s
nexttile;
if size(nq_50_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(nq_50_env_6_time(plot_indices), nq_50_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nq_50_env_6_time(plot_indices), nq_50_env_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(q_50_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('q_50', 'FontSize', subtitle_fontsize);
else
    plot(q_50_env_6_time(plot_indices), q_50_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(q_50_env_6_time(plot_indices), q_50_env_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(nn_50_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('Our method', 'FontSize', subtitle_fontsize);
else
    plot(nn_50_env_6_time(plot_indices), nn_50_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nn_50_env_6_time(plot_indices), nn_50_env_6_command_data(plot_indices, 1));
    % title('Our method', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.02;
    % lg.Position(2) = lg.Position(2) - 0.06;
end

han = axes(gcf, 'Visible', 'off'); % Create an invisible axes that spans the figure
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
xlabel(han, 'Time, [s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [0.5, -0.05, 0]); % Global x-axis label
ylabel(han, 'Velocity, [m/s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [-0.07, 0.5, 0]);


%% Noise lim 20, NN methods
figure('Name', 'NN method comparison at lim 20Nm, (No quad, quad, explicit');
set(gcf, 'Units', 'inches', 'Position', [1, 1, double_column_eight_row_size]); % Choose single_column_size or double_column_size here
% sgtitle("Base Forward Velocity Tracking Performance in Nominal Settings", 'FontSize', sgtitle_fontsize);
t = tiledlayout(8, 3, 'TileSpacing', 'compact', 'Padding', 'normal'); % Adjust 'TileSpacing' and 'Padding' as needed
nexttile;
plot(nq_20_nom_4_time(plot_indices), nq_20_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(nq_20_nom_4_time(plot_indices), nq_20_nom_4_command_data(plot_indices, 1));
title('No quadratic term', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');

nexttile;
plot(q_20_nom_4_time(plot_indices), q_20_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(q_20_nom_4_time(plot_indices), q_20_nom_4_command_data(plot_indices, 1));
title('Quadratic term', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast')

nexttile;
plot(nn_20_nom_4_time(plot_indices), nn_20_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(nn_20_nom_4_time(plot_indices), nn_20_nom_4_command_data(plot_indices, 1));
title('Our method', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
lgd.Position(1) = lgd.Position(1) + 0.02;
lgd.Position(2) = lgd.Position(2) - 0.01;

% sgtitle("Base Forward Velocity Tracking Performance in Nominal Settings", 'FontSize', sgtitle_fontsize);
nexttile;
plot(nq_20_nom_6_time(plot_indices), nq_20_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(nq_20_nom_6_time(plot_indices), nq_20_nom_6_command_data(plot_indices, 1));
% title('No quadratic term', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'southeast')

nexttile;
plot(q_20_nom_6_time(plot_indices), q_20_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(q_20_nom_6_time(plot_indices), q_20_nom_6_command_data(plot_indices, 1));
% title('q_20', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'southeast')

nexttile;
plot(nn_20_nom_6_time(plot_indices), nn_20_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(nn_20_nom_6_time(plot_indices), nn_20_nom_6_command_data(plot_indices, 1));
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'southeast')


% Compare No quadratic term, q_20, and our method under actuator gap
nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(nq_20_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(nq_20_act_4_time(plot_indices), nq_20_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nq_20_act_4_time(plot_indices), nq_20_act_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(q_20_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(q_20_act_4_time(plot_indices), q_20_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(q_20_act_4_time(plot_indices), q_20_act_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(nn_20_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(nn_20_act_4_time(plot_indices), nn_20_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nn_20_act_4_time(plot_indices), nn_20_act_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
end

% 0.6m/s
nexttile;
if size(nq_20_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('DR', 'FontSize', subtitle_fontsize);
else
    plot(nq_20_act_6_time(plot_indices), nq_20_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nq_20_act_6_time(plot_indices), nq_20_act_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(q_20_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('q_20', 'FontSize', subtitle_fontsize);
else
    plot(q_20_act_6_time(plot_indices), q_20_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(q_20_act_6_time(plot_indices), q_20_act_6_command_data(plot_indices, 1));
    % title('q_20', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(nn_20_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('Our method', 'FontSize', subtitle_fontsize);
else
    plot(nn_20_act_6_time(plot_indices), nn_20_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nn_20_act_6_time(plot_indices), nn_20_act_6_command_data(plot_indices, 1));
    % title('Our method', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.02;
    % lg.Position(2) = lg.Position(2) - 0.06;
end


% Compare No quadratic term, q_20, and our method under model gap%%%

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(nq_20_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(nq_20_mod_4_time(plot_indices), nq_20_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nq_20_mod_4_time(plot_indices), nq_20_mod_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(q_20_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(q_20_mod_4_time(plot_indices), q_20_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(q_20_mod_4_time(plot_indices), q_20_mod_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(nn_20_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(nn_20_mod_4_time(plot_indices), nn_20_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nn_20_mod_4_time(plot_indices), nn_20_mod_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
end

% 0.6m/s
nexttile;
if size(nq_20_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(nq_20_mod_6_time(plot_indices), nq_20_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nq_20_mod_6_time(plot_indices), nq_20_mod_6_command_data(plot_indices, 1));
    % title('DR', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(q_20_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('q_20', 'FontSize', subtitle_fontsize);
else
    plot(q_20_mod_6_time(plot_indices), q_20_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(q_20_mod_6_time(plot_indices), q_20_mod_6_command_data(plot_indices, 1));
    % title('q_20', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(nn_20_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('Our method', 'FontSize', subtitle_fontsize);
else
    plot(nn_20_mod_6_time(plot_indices), nn_20_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nn_20_mod_6_time(plot_indices), nn_20_mod_6_command_data(plot_indices, 1));
    % title('Our method', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.02;
    % lg.Position(2) = lg.Position(2) - 0.06;
end




% Compare No quadratic term, q_20, and our method under environment gap%%%
nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(nq_20_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(nq_20_env_4_time(plot_indices), nq_20_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nq_20_env_4_time(plot_indices), nq_20_env_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(q_20_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(q_20_env_4_time(plot_indices), q_20_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(q_20_env_4_time(plot_indices), q_20_env_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(nn_20_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(nn_20_env_4_time(plot_indices), nn_20_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nn_20_env_4_time(plot_indices), nn_20_env_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
end

% 0.6m/s
nexttile;
if size(nq_20_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(nq_20_env_6_time(plot_indices), nq_20_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nq_20_env_6_time(plot_indices), nq_20_env_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(q_20_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('q_20', 'FontSize', subtitle_fontsize);
else
    plot(q_20_env_6_time(plot_indices), q_20_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(q_20_env_6_time(plot_indices), q_20_env_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(nn_20_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('Our method', 'FontSize', subtitle_fontsize);
else
    plot(nn_20_env_6_time(plot_indices), nn_20_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nn_20_env_6_time(plot_indices), nn_20_env_6_command_data(plot_indices, 1));
    % title('Our method', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.02;
    % lg.Position(2) = lg.Position(2) - 0.06;
end

han = axes(gcf, 'Visible', 'off'); % Create an invisible axes that spans the figure
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
xlabel(han, 'Time, [s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [0.5, -0.05, 0]); % Global x-axis label
ylabel(han, 'Velocity, [m/s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [-0.07, 0.5, 0]);

%% Evaluate ERFI varying noise lim
figure('Name', 'ERFI comparison, (20Nm, 30Nm, 50Nm');
set(gcf, 'Units', 'inches', 'Position', [1, 1, double_column_eight_row_size]); % Choose single_column_size or double_column_size here
% sgtitle("Base Forward Velocity Tracking Performance in Nominal Settings", 'FontSize', sgtitle_fontsize);
t = tiledlayout(8, 3, 'TileSpacing', 'compact', 'Padding', 'normal'); % Adjust 'TileSpacing' and 'Padding' as needed
nexttile;
plot(erfi_20_nom_4_time(plot_indices), erfi_20_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(erfi_20_nom_4_time(plot_indices), erfi_20_nom_4_command_data(plot_indices, 1));
title('ERFI, 20Nm', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');

nexttile;
plot(erfi_30_nom_4_time(plot_indices), erfi_30_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(erfi_30_nom_4_time(plot_indices), erfi_30_nom_4_command_data(plot_indices, 1));
title('ERFI, 30Nm', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast')

nexttile;
plot(erfi_50_nom_4_time(plot_indices), erfi_50_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(erfi_50_nom_4_time(plot_indices), erfi_50_nom_4_command_data(plot_indices, 1));
title('ERFI, 50Nm', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
lgd.Position(1) = lgd.Position(1) + 0.02;
lgd.Position(2) = lgd.Position(2) - 0.01;

% sgtitle("Base Forward Velocity Tracking Performance in Nominal Settings", 'FontSize', sgtitle_fontsize);
nexttile;
plot(erfi_20_nom_6_time(plot_indices), erfi_20_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(erfi_20_nom_6_time(plot_indices), erfi_20_nom_6_command_data(plot_indices, 1));
% title('No quadratic term', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'southeast')

nexttile;
plot(erfi_30_nom_6_time(plot_indices), erfi_30_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(erfi_30_nom_6_time(plot_indices), erfi_30_nom_6_command_data(plot_indices, 1));
% title('erfi_30', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'southeast')

nexttile;
plot(erfi_50_nom_6_time(plot_indices), erfi_50_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(erfi_50_nom_6_time(plot_indices), erfi_50_nom_6_command_data(plot_indices, 1));
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'southeast')


% Compare No quadratic term, erfi_30, and our method under actuator gap
nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(erfi_20_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(erfi_20_act_4_time(plot_indices), erfi_20_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_20_act_4_time(plot_indices), erfi_20_act_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(erfi_30_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(erfi_30_act_4_time(plot_indices), erfi_30_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_30_act_4_time(plot_indices), erfi_30_act_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(erfi_50_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(erfi_50_act_4_time(plot_indices), erfi_50_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_50_act_4_time(plot_indices), erfi_50_act_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
end

% 0.6m/s
nexttile;
if size(erfi_20_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('DR', 'FontSize', subtitle_fontsize);
else
    plot(erfi_20_act_6_time(plot_indices), erfi_20_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_20_act_6_time(plot_indices), erfi_20_act_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(erfi_30_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('erfi_30', 'FontSize', subtitle_fontsize);
else
    plot(erfi_30_act_6_time(plot_indices), erfi_30_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_30_act_6_time(plot_indices), erfi_30_act_6_command_data(plot_indices, 1));
    % title('erfi_30', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(erfi_50_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('Our method', 'FontSize', subtitle_fontsize);
else
    plot(erfi_50_act_6_time(plot_indices), erfi_50_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_50_act_6_time(plot_indices), erfi_50_act_6_command_data(plot_indices, 1));
    % title('Our method', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.02;
    % lg.Position(2) = lg.Position(2) - 0.06;
end


% Compare No quadratic term, erfi_30, and our method under model gap%%%

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(erfi_20_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(erfi_20_mod_4_time(plot_indices), erfi_20_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_20_mod_4_time(plot_indices), erfi_20_mod_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(erfi_30_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(erfi_30_mod_4_time(plot_indices), erfi_30_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_30_mod_4_time(plot_indices), erfi_30_mod_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(erfi_50_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(erfi_50_mod_4_time(plot_indices), erfi_50_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_50_mod_4_time(plot_indices), erfi_50_mod_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
end

% 0.6m/s
nexttile;
if size(erfi_20_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(erfi_20_mod_6_time(plot_indices), erfi_20_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_20_mod_6_time(plot_indices), erfi_20_mod_6_command_data(plot_indices, 1));
    % title('DR', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(erfi_30_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('erfi_30', 'FontSize', subtitle_fontsize);
else
    plot(erfi_30_mod_6_time(plot_indices), erfi_30_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_30_mod_6_time(plot_indices), erfi_30_mod_6_command_data(plot_indices, 1));
    % title('erfi_30', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(erfi_50_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('Our method', 'FontSize', subtitle_fontsize);
else
    plot(erfi_50_mod_6_time(plot_indices), erfi_50_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_50_mod_6_time(plot_indices), erfi_50_mod_6_command_data(plot_indices, 1));
    % title('Our method', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.02;
    % lg.Position(2) = lg.Position(2) - 0.06;
end




% Compare No quadratic term, erfi_30, and our method under environment gap%%%
nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(erfi_20_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(erfi_20_env_4_time(plot_indices), erfi_20_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_20_env_4_time(plot_indices), erfi_20_env_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(erfi_30_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(erfi_30_env_4_time(plot_indices), erfi_30_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_30_env_4_time(plot_indices), erfi_30_env_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(erfi_50_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(erfi_50_env_4_time(plot_indices), erfi_50_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_50_env_4_time(plot_indices), erfi_50_env_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
end

% 0.6m/s
nexttile;
if size(erfi_20_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
else
    plot(erfi_20_env_6_time(plot_indices), erfi_20_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_20_env_6_time(plot_indices), erfi_20_env_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(erfi_30_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('erfi_30', 'FontSize', subtitle_fontsize);
else
    plot(erfi_30_env_6_time(plot_indices), erfi_30_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_30_env_6_time(plot_indices), erfi_30_env_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(erfi_50_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    % title('Our method', 'FontSize', subtitle_fontsize);
else
    plot(erfi_50_env_6_time(plot_indices), erfi_50_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_50_env_6_time(plot_indices), erfi_50_env_6_command_data(plot_indices, 1));
    % title('Our method', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.02;
    % lg.Position(2) = lg.Position(2) - 0.06;
end

han = axes(gcf, 'Visible', 'off'); % Create an invisible axes that spans the figure
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
xlabel(han, 'Time, [s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [0.5, -0.05, 0]); % Global x-axis label
ylabel(han, 'Velocity, [m/s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [-0.07, 0.5, 0]);

