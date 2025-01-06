close all; clc;

% %% Load Data. Do only once
% %% Load DR data
% clear all;
% dr_act_4_data = load('DR/actuator_gap/0.4ms/data.csv');
% dr_act_6_data = load('DR/actuator_gap/0.6ms/data.csv');
% dr_env_4_data = load('DR/environment_gap/0.4ms/data.csv');
% dr_env_6_data = load('DR/environment_gap/0.6ms/data.csv');
% dr_mod_4_data = load('DR/model_gap/0.4ms/data.csv');
% dr_mod_6_data = load('DR/model_gap/0.6ms/data.csv');
% dr_nom_4_data = load('DR/nominal/0.4ms/data.csv');
% dr_nom_6_data = load('DR/nominal/0.6ms/data.csv');
% %%
% %% Load noise data
% erfi_50_act_4_data = load('ERFI_50_uniform/actuator_gap/0.4/data.csv');
% erfi_50_act_6_data = load('ERFI_50_uniform/actuator_gap/0.6/data.csv');
% erfi_50_env_4_data = load('ERFI_50_uniform/environment_gap/0.4/data.csv');
% erfi_50_env_6_data = load('ERFI_50_uniform/environment_gap/0.6/data.csv');
% erfi_50_mod_4_data = load('ERFI_50_uniform/model_gap/0.4/data.csv');
% erfi_50_mod_6_data = load('ERFI_50_uniform/model_gap/0.6/data.csv');
% erfi_50_nom_4_data = load('ERFI_50_uniform/nominal/0.4/data.csv');
% erfi_50_nom_6_data = load('ERFI_50_uniform/nominal/0.6/data.csv');
% 
% erfi_beta_50_act_4_data = load('ERFI_50_beta/actuator_gap/0.4/data.csv');
% erfi_beta_50_act_6_data = load('ERFI_50_beta/actuator_gap/0.6/data.csv');
% erfi_beta_50_env_4_data = load('ERFI_50_beta/environment_gap/0.4/data.csv');
% erfi_beta_50_env_6_data = load('ERFI_50_beta/environment_gap/0.6/data.csv');
% erfi_beta_50_mod_4_data = load('ERFI_50_beta/model_gap/0.4/data.csv');
% erfi_beta_50_mod_6_data = load('ERFI_50_beta/model_gap/0.6/data.csv');
% erfi_beta_50_nom_4_data = load('ERFI_50_beta/nominal/0.4/data.csv');
% erfi_beta_50_nom_6_data = load('ERFI_50_beta/nominal/0.6/data.csv');
% 
% erfi_80_act_4_data = load('ERFI_80_uniform/actuator_gap/0.4/data.csv');
% erfi_80_act_6_data = load('ERFI_80_uniform/actuator_gap/0.6/data.csv');
% erfi_80_env_4_data = load('ERFI_80_uniform/environment_gap/0.4/data.csv');
% erfi_80_env_6_data = load('ERFI_80_uniform/environment_gap/0.6/data.csv');
% erfi_80_mod_4_data = load('ERFI_80_uniform/model_gap/0.4/data.csv');
% erfi_80_mod_6_data = load('ERFI_80_uniform/model_gap/0.6/data.csv');
% erfi_80_nom_4_data = load('ERFI_80_uniform/nominal/0.4/data.csv');
% erfi_80_nom_6_data = load('ERFI_80_uniform/nominal/0.6/data.csv');
% 
% erfi_beta_80_act_4_data = load('ERFI_80_beta/actuator_gap/0.4/data.csv');
% erfi_beta_80_act_6_data = load('ERFI_80_beta/actuator_gap/0.6/data.csv');
% erfi_beta_80_env_4_data = load('ERFI_80_beta/environment_gap/0.4/data.csv');
% erfi_beta_80_env_6_data = load('ERFI_80_beta/environment_gap/0.6/data.csv');
% erfi_beta_80_mod_4_data = load('ERFI_80_beta/model_gap/0.4/data.csv');
% erfi_beta_80_mod_6_data = load('ERFI_80_beta/model_gap/0.6/data.csv');
% erfi_beta_80_nom_4_data = load('ERFI_80_beta/nominal/0.4/data.csv');
% erfi_beta_80_nom_6_data = load('ERFI_80_beta/nominal/0.6/data.csv');
% 
% nq_50_act_4_data = load('NQ_50_gaussian_1.5/actuator_gap/0.4/data.csv');
% nq_50_act_6_data = load('NQ_50_gaussian_1.5/actuator_gap/0.6/data.csv');
% nq_50_env_4_data = load('NQ_50_gaussian_1.5/environment_gap/0.4/data.csv');
% nq_50_env_6_data = load('NQ_50_gaussian_1.5/environment_gap/0.6/data.csv');
% nq_50_mod_4_data = load('NQ_50_gaussian_1.5/model_gap/0.4/data.csv');
% nq_50_mod_6_data = load('NQ_50_gaussian_1.5/model_gap/0.6/data.csv');
% nq_50_nom_4_data = load('NQ_50_gaussian_1.5/nominal/0.4/data.csv');
% nq_50_nom_6_data = load('NQ_50_gaussian_1.5/nominal/0.6/data.csv');
% 
% nq_80_act_4_data = load('NQ_80_gaussian_1.5/actuator_gap/0.4/data.csv');
% nq_80_act_6_data = load('NQ_80_gaussian_1.5/actuator_gap/0.6/data.csv');
% nq_80_env_4_data = load('NQ_80_gaussian_1.5/environment_gap/0.4/data.csv');
% nq_80_env_6_data = load('NQ_80_gaussian_1.5/environment_gap/0.6/data.csv');
% nq_80_mod_4_data = load('NQ_80_gaussian_1.5/model_gap/0.4/data.csv');
% nq_80_mod_6_data = load('NQ_80_gaussian_1.5/model_gap/0.6/data.csv');
% nq_80_nom_4_data = load('NQ_80_gaussian_1.5/nominal/0.4/data.csv');
% nq_80_nom_6_data = load('NQ_80_gaussian_1.5/nominal/0.6/data.csv');
% 
% 
% erfi_20_act_4_data = load('ERFI_20_uniform/actuator_gap/0.4/data.csv');
% erfi_20_act_6_data = load('ERFI_20_uniform/actuator_gap/0.6/data.csv');
% erfi_20_env_4_data = load('ERFI_20_uniform/environment_gap/0.4/data.csv');
% erfi_20_env_6_data = load('ERFI_20_uniform/environment_gap/0.6/data.csv');
% erfi_20_mod_4_data = load('ERFI_20_uniform/model_gap/0.4/data.csv');
% erfi_20_mod_6_data = load('ERFI_20_uniform/model_gap/0.6/data.csv');
% erfi_20_nom_4_data = load('ERFI_20_uniform/nominal/0.4/data.csv');
% erfi_20_nom_6_data = load('ERFI_20_uniform/nominal/0.6/data.csv');
% 
% erfi_beta_20_act_4_data = load('ERFI_20_uniform/actuator_gap/0.4/data.csv');
% erfi_beta_20_act_6_data = load('ERFI_20_uniform/actuator_gap/0.6/data.csv');
% erfi_beta_20_env_4_data = load('ERFI_20_uniform/environment_gap/0.4/data.csv');
% erfi_beta_20_env_6_data = load('ERFI_20_uniform/environment_gap/0.6/data.csv');
% erfi_beta_20_mod_4_data = load('ERFI_20_uniform/model_gap/0.4/data.csv');
% erfi_beta_20_mod_6_data = load('ERFI_20_uniform/model_gap/0.6/data.csv');
% erfi_beta_20_nom_4_data = load('ERFI_20_uniform/nominal/0.4/data.csv');
% erfi_beta_20_nom_6_data = load('ERFI_20_uniform/nominal/0.6/data.csv');
% 
% nq_20_act_4_data = load('NQ_20_gaussian_1.5/actuator_gap/0.4/data.csv');
% nq_20_act_6_data = load('NQ_20_gaussian_1.5/actuator_gap/0.6/data.csv');
% nq_20_env_4_data = load('NQ_20_gaussian_1.5/environment_gap/0.4/data.csv');
% nq_20_env_6_data = load('NQ_20_gaussian_1.5/environment_gap/0.6/data.csv');
% nq_20_mod_4_data = load('NQ_20_gaussian_1.5/model_gap/0.4/data.csv');
% nq_20_mod_6_data = load('NQ_20_gaussian_1.5/model_gap/0.6/data.csv');
% nq_20_nom_4_data = load('NQ_20_gaussian_1.5/nominal/0.4/data.csv');
% nq_20_nom_6_data = load('NQ_20_gaussian_1.5/nominal/0.6/data.csv');
% 
% 
% hy_50_act_4_data = load('Hybrid_50/actuator_gap/0.4/data.csv');
% hy_50_act_6_data = load('Hybrid_50/actuator_gap/0.6/data.csv');
% hy_50_env_4_data = load('Hybrid_50/environment_gap/0.4/data.csv');
% hy_50_env_6_data = load('Hybrid_50/environment_gap/0.6/data.csv');
% hy_50_mod_4_data = load('Hybrid_50/model_gap/0.4/data.csv');
% hy_50_mod_6_data = load('Hybrid_50/model_gap/0.6/data.csv');
% hy_50_nom_4_data = load('Hybrid_50/nominal/0.4/data.csv');
% hy_50_nom_6_data = load('Hybrid_50/nominal/0.6/data.csv');
% 
% hy_100_act_4_data = load('Hybrid_100/actuator_gap/0.4/data.csv');
% hy_100_act_6_data = load('Hybrid_100/actuator_gap/0.6/data.csv');
% hy_100_env_4_data = load('Hybrid_100/environment_gap/0.4/data.csv');
% hy_100_env_6_data = load('Hybrid_100/environment_gap/0.6/data.csv');
% hy_100_mod_4_data = load('Hybrid_100/model_gap/0.4/data.csv');
% hy_100_mod_6_data = load('Hybrid_100/model_gap/0.6/data.csv');
% hy_100_nom_4_data = load('Hybrid_100/nominal/0.4/data.csv');
% hy_100_nom_6_data = load('Hybrid_100/nominal/0.6/data.csv');
% %%
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
%%
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

process_data(erfi_beta_50_act_4_data, 'erfi_beta_50_act_4');
process_data(erfi_beta_50_act_6_data, 'erfi_beta_50_act_6');
process_data(erfi_beta_50_env_4_data, 'erfi_beta_50_env_4');
process_data(erfi_beta_50_env_6_data, 'erfi_beta_50_env_6');
process_data(erfi_beta_50_mod_4_data, 'erfi_beta_50_mod_4');
process_data(erfi_beta_50_mod_6_data, 'erfi_beta_50_mod_6');
process_data(erfi_beta_50_nom_4_data, 'erfi_beta_50_nom_4');
process_data(erfi_beta_50_nom_6_data, 'erfi_beta_50_nom_6');

process_data(nq_50_act_4_data, 'nq_50_act_4');
process_data(nq_50_act_6_data, 'nq_50_act_6');
process_data(nq_50_env_4_data, 'nq_50_env_4');
process_data(nq_50_env_6_data, 'nq_50_env_6');
process_data(nq_50_mod_4_data, 'nq_50_mod_4');
process_data(nq_50_mod_6_data, 'nq_50_mod_6');
process_data(nq_50_nom_4_data, 'nq_50_nom_4');
process_data(nq_50_nom_6_data, 'nq_50_nom_6');


process_data(erfi_80_act_4_data, 'erfi_80_act_4');
process_data(erfi_80_act_6_data, 'erfi_80_act_6');
process_data(erfi_80_env_4_data, 'erfi_80_env_4');
process_data(erfi_80_env_6_data, 'erfi_80_env_6');
process_data(erfi_80_mod_4_data, 'erfi_80_mod_4');
process_data(erfi_80_mod_6_data, 'erfi_80_mod_6');
process_data(erfi_80_nom_4_data, 'erfi_80_nom_4');
process_data(erfi_80_nom_6_data, 'erfi_80_nom_6');

process_data(erfi_beta_80_act_4_data, 'erfi_beta_80_act_4');
process_data(erfi_beta_80_act_6_data, 'erfi_beta_80_act_6');
process_data(erfi_beta_80_env_4_data, 'erfi_beta_80_env_4');
process_data(erfi_beta_80_env_6_data, 'erfi_beta_80_env_6');
process_data(erfi_beta_80_mod_4_data, 'erfi_beta_80_mod_4');
process_data(erfi_beta_80_mod_6_data, 'erfi_beta_80_mod_6');
process_data(erfi_beta_80_nom_4_data, 'erfi_beta_80_nom_4');
process_data(erfi_beta_80_nom_6_data, 'erfi_beta_80_nom_6');

process_data(nq_80_act_4_data, 'nq_80_act_4');
process_data(nq_80_act_6_data, 'nq_80_act_6');
process_data(nq_80_env_4_data, 'nq_80_env_4');
process_data(nq_80_env_6_data, 'nq_80_env_6');
process_data(nq_80_mod_4_data, 'nq_80_mod_4');
process_data(nq_80_mod_6_data, 'nq_80_mod_6');
process_data(nq_80_nom_4_data, 'nq_80_nom_4');
process_data(nq_80_nom_6_data, 'nq_80_nom_6');

process_data(erfi_20_act_4_data, 'erfi_20_act_4');
process_data(erfi_20_act_6_data, 'erfi_20_act_6');
process_data(erfi_20_env_4_data, 'erfi_20_env_4');
process_data(erfi_20_env_6_data, 'erfi_20_env_6');
process_data(erfi_20_mod_4_data, 'erfi_20_mod_4');
process_data(erfi_20_mod_6_data, 'erfi_20_mod_6');
process_data(erfi_20_nom_4_data, 'erfi_20_nom_4');
process_data(erfi_20_nom_6_data, 'erfi_20_nom_6');

process_data(erfi_beta_20_act_4_data, 'erfi_beta_20_act_4');
process_data(erfi_beta_20_act_6_data, 'erfi_beta_20_act_6');
process_data(erfi_beta_20_env_4_data, 'erfi_beta_20_env_4');
process_data(erfi_beta_20_env_6_data, 'erfi_beta_20_env_6');
process_data(erfi_beta_20_mod_4_data, 'erfi_beta_20_mod_4');
process_data(erfi_beta_20_mod_6_data, 'erfi_beta_20_mod_6');
process_data(erfi_beta_20_nom_4_data, 'erfi_beta_20_nom_4');
process_data(erfi_beta_20_nom_6_data, 'erfi_beta_20_nom_6');

process_data(nq_20_act_4_data, 'nq_20_act_4');
process_data(nq_20_act_6_data, 'nq_20_act_6');
process_data(nq_20_env_4_data, 'nq_20_env_4');
process_data(nq_20_env_6_data, 'nq_20_env_6');
process_data(nq_20_mod_4_data, 'nq_20_mod_4');
process_data(nq_20_mod_6_data, 'nq_20_mod_6');
process_data(nq_20_nom_4_data, 'nq_20_nom_4');
process_data(nq_20_nom_6_data, 'nq_20_nom_6');

process_data(hy_50_act_4_data, 'hy_50_act_4');
process_data(hy_50_act_6_data, 'hy_50_act_6');
process_data(hy_50_env_4_data, 'hy_50_env_4');
process_data(hy_50_env_6_data, 'hy_50_env_6');
process_data(hy_50_mod_4_data, 'hy_50_mod_4');
process_data(hy_50_mod_6_data, 'hy_50_mod_6');
process_data(hy_50_nom_4_data, 'hy_50_nom_4');
process_data(hy_50_nom_6_data, 'hy_50_nom_6');

process_data(hy_100_act_4_data, 'hy_100_act_4');
process_data(hy_100_act_6_data, 'hy_100_act_6');
process_data(hy_100_env_4_data, 'hy_100_env_4');
process_data(hy_100_env_6_data, 'hy_100_env_6');
process_data(hy_100_mod_4_data, 'hy_100_mod_4');
process_data(hy_100_mod_6_data, 'hy_100_mod_6');
process_data(hy_100_nom_4_data, 'hy_100_nom_4');
process_data(hy_100_nom_6_data, 'hy_100_nom_6');

%%
%% Final baselines

process_data(dr_act_4_data, 'DR_act_4');
process_data(dr_act_6_data, 'DR_act_6');
process_data(dr_env_4_data, 'DR_env_4');
process_data(dr_env_6_data, 'DR_env_6');
process_data(dr_mod_4_data, 'DR_mod_4');
process_data(dr_mod_6_data, 'DR_mod_6');
process_data(dr_nom_4_data, 'DR_nom_4');
process_data(dr_nom_6_data, 'DR_nom_6');

process_data(erfi_beta_50_act_4_data, 'ERFI_act_4');
process_data(erfi_beta_50_act_6_data, 'ERFI_act_6');
process_data(erfi_beta_50_env_4_data, 'ERFI_env_4');
process_data(erfi_beta_50_env_6_data, 'ERFI_env_6');
process_data(erfi_beta_50_mod_4_data, 'ERFI_mod_4');
process_data(erfi_beta_50_mod_6_data, 'ERFI_mod_6');
process_data(erfi_beta_50_nom_4_data, 'ERFI_nom_4');
process_data(erfi_beta_50_nom_6_data, 'ERFI_nom_6');

process_data(nq_50_act_4_data, 'OUR_act_4');
process_data(nq_50_act_6_data, 'OUR_act_6');
process_data(nq_50_env_4_data, 'OUR_env_4');
process_data(nq_50_env_6_data, 'OUR_env_6');
process_data(nq_50_mod_4_data, 'OUR_mod_4');
process_data(nq_50_mod_6_data, 'OUR_mod_6');
process_data(nq_50_nom_4_data, 'OUR_nom_4');
process_data(nq_50_nom_6_data, 'OUR_nom_6');
%%
%% Set plotting params

plot_start_time = 0.008;
plot_end_time = 10;
control_hz = 125;
plot_indices = control_hz*plot_start_time:control_hz*plot_end_time;
% Define figure size for single or double column
single_column_size = [4, 1.5]; % [width, height] in inches for single column
double_column_size = [7, 3]; % [width, height] in inches for double column
double_column_single_row_size = [7,1.5];
double_column_eight_row_size = [7, 12];
single_column_eight_row_size = [4, 12];
ylimit = [0, 0.7];
sgtitle_fontsize = 10;
subtitle_fontsize = 6;
label_fontsize = 4;
legend_fontsize = 4;
tick_fontsize = 4;
%%
%% Evaluate ERFI per noise lim
figure('Name', 'ERFI uniform');
set(gcf, 'Units', 'inches', 'Position', [1, 1, double_column_eight_row_size]); % Choose single_column_size or double_column_size here
% sgtitle("Base Forward Velocity Tracking Performance in Nominal Settings", 'FontSize', sgtitle_fontsize);
t = tiledlayout(8, 3, 'TileSpacing', 'compact', 'Padding', 'normal'); % Adjust 'TileSpacing' and 'Padding' as needed
nexttile;
plot(erfi_20_nom_4_time(plot_indices), erfi_20_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(erfi_20_nom_4_time(plot_indices), erfi_20_nom_4_command_data(plot_indices, 1));
title('ERFI 20', 'FontSize', subtitle_fontsize);
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
title('ERFI 50', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');

nexttile;
plot(erfi_80_nom_4_time(plot_indices), erfi_80_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(erfi_80_nom_4_time(plot_indices), erfi_80_nom_4_command_data(plot_indices, 1));
title('ERFI 80', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
lgd.Position(1) = lgd.Position(1) + 0.02;
lgd.Position(2) = lgd.Position(2) - 0.01;

% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast')

nexttile;
plot(erfi_20_nom_6_time(plot_indices), erfi_20_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(erfi_20_nom_6_time(plot_indices), erfi_20_nom_6_command_data(plot_indices, 1));
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;

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

% sgtitle("Base Forward Velocity Tracking Performance in Nominal Settings", 'FontSize', sgtitle_fontsize);
nexttile;
plot(erfi_80_nom_6_time(plot_indices), erfi_80_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(erfi_80_nom_6_time(plot_indices), erfi_80_nom_6_command_data(plot_indices, 1));
% title('No quadratic term', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'southeast')



% Compare No quadratic term, erfi_20, and our method under actuator gap
nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(erfi_20_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
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

% 0.4m/s
if size(erfi_50_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
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
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(erfi_80_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_80_act_4_time(plot_indices), erfi_80_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_80_act_4_time(plot_indices), erfi_80_act_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

% 0.6m/s
nexttile;
if size(erfi_20_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
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

% 0.4m/s
if size(erfi_50_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_50_act_6_time(plot_indices), erfi_50_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_50_act_6_time(plot_indices), erfi_50_act_6_command_data(plot_indices, 1));
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

if size(erfi_80_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    % title('erfi_20', 'FontSize', subtitle_fontsize);
else
    plot(erfi_80_act_6_time(plot_indices), erfi_80_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_80_act_6_time(plot_indices), erfi_80_act_6_command_data(plot_indices, 1));
    % title('erfi_20', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end


% Compare No quadratic term, erfi_20, and our method under model gap%%%

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(erfi_20_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
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

% 0.4m/s
if size(erfi_50_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
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
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end


nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(erfi_80_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_80_mod_4_time(plot_indices), erfi_80_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_80_mod_4_time(plot_indices), erfi_80_mod_4_command_data(plot_indices, 1));
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

if size(erfi_20_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_20_mod_6_time(plot_indices), erfi_20_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_20_mod_6_time(plot_indices), erfi_20_mod_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
end


nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(erfi_50_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_50_mod_6_time(plot_indices), erfi_50_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_50_mod_6_time(plot_indices), erfi_50_mod_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;
if size(erfi_80_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_80_mod_6_time(plot_indices), erfi_80_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_80_mod_6_time(plot_indices), erfi_80_mod_6_command_data(plot_indices, 1));
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




% Compare No quadratic term, erfi_20, and our method under environment gap%%%
nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(erfi_20_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
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

% 0.4m/s
if size(erfi_50_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
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
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

if size(erfi_80_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_80_env_4_time(plot_indices), erfi_80_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_80_env_4_time(plot_indices), erfi_80_env_4_command_data(plot_indices, 1));
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

if size(erfi_20_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
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
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(erfi_50_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_50_env_6_time(plot_indices), erfi_50_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_50_env_6_time(plot_indices), erfi_50_env_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end


nexttile;
if size(erfi_80_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_80_env_6_time(plot_indices), erfi_80_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_80_env_6_time(plot_indices), erfi_80_env_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end


han = axes(gcf, 'Visible', 'off'); % Create an invisible axes that spans the figure
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
xlabel(han, 'Time, [s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [0.5, -0.05, 0]); % Global x-axis label
ylabel(han, 'Velocity, [m/s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [-0.07, 0.5, 0]);
%%
%% Evaluate ERFI-beta per noise lim
figure('Name', 'ERFI beta');
set(gcf, 'Units', 'inches', 'Position', [1, 1, double_column_eight_row_size]); % Choose single_column_size or double_column_size here
% sgtitle("Base Forward Velocity Tracking Performance in Nominal Settings", 'FontSize', sgtitle_fontsize);
t = tiledlayout(8, 3, 'TileSpacing', 'compact', 'Padding', 'normal'); % Adjust 'TileSpacing' and 'Padding' as needed
nexttile;
plot(erfi_beta_20_nom_4_time(plot_indices), erfi_beta_20_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(erfi_beta_20_nom_4_time(plot_indices), erfi_beta_20_nom_4_command_data(plot_indices, 1));
title('erfi beta 20', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
nexttile;
plot(erfi_beta_50_nom_4_time(plot_indices), erfi_beta_50_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(erfi_beta_50_nom_4_time(plot_indices), erfi_beta_50_nom_4_command_data(plot_indices, 1));
title('erfi beta 50', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');

nexttile;
plot(erfi_beta_80_nom_4_time(plot_indices), erfi_beta_80_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(erfi_beta_80_nom_4_time(plot_indices), erfi_beta_80_nom_4_command_data(plot_indices, 1));
title('erfi beta 80', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
lgd.Position(1) = lgd.Position(1) + 0.02;
lgd.Position(2) = lgd.Position(2) - 0.01;

% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast')

nexttile;
plot(erfi_beta_20_nom_6_time(plot_indices), erfi_beta_20_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(erfi_beta_20_nom_6_time(plot_indices), erfi_beta_20_nom_6_command_data(plot_indices, 1));
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;

nexttile;
plot(erfi_beta_50_nom_6_time(plot_indices), erfi_beta_50_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(erfi_beta_50_nom_6_time(plot_indices), erfi_beta_50_nom_6_command_data(plot_indices, 1));
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;

% sgtitle("Base Forward Velocity Tracking Performance in Nominal Settings", 'FontSize', sgtitle_fontsize);
nexttile;
plot(erfi_beta_80_nom_6_time(plot_indices), erfi_beta_80_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(erfi_beta_80_nom_6_time(plot_indices), erfi_beta_80_nom_6_command_data(plot_indices, 1));
% title('No quadratic term', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'southeast')



% Compare No quadratic term, erfi_beta_20, and our method under actuator gap
nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(erfi_beta_20_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_beta_20_act_4_time(plot_indices), erfi_beta_20_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_beta_20_act_4_time(plot_indices), erfi_beta_20_act_4_command_data(plot_indices, 1));
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

% 0.4m/s
if size(erfi_beta_50_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_beta_50_act_4_time(plot_indices), erfi_beta_50_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_beta_50_act_4_time(plot_indices), erfi_beta_50_act_4_command_data(plot_indices, 1));
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

if size(erfi_beta_80_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_beta_80_act_4_time(plot_indices), erfi_beta_80_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_beta_80_act_4_time(plot_indices), erfi_beta_80_act_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

% 0.6m/s
nexttile;
if size(erfi_beta_20_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    % title('DR', 'FontSize', subtitle_fontsize);
else
    plot(erfi_beta_20_act_6_time(plot_indices), erfi_beta_20_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_beta_20_act_6_time(plot_indices), erfi_beta_20_act_6_command_data(plot_indices, 1));
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

% 0.4m/s
if size(erfi_beta_50_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_beta_50_act_6_time(plot_indices), erfi_beta_50_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_beta_50_act_6_time(plot_indices), erfi_beta_50_act_6_command_data(plot_indices, 1));
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

if size(erfi_beta_80_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    % title('erfi_beta_20', 'FontSize', subtitle_fontsize);
else
    plot(erfi_beta_80_act_6_time(plot_indices), erfi_beta_80_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_beta_80_act_6_time(plot_indices), erfi_beta_80_act_6_command_data(plot_indices, 1));
    % title('erfi_beta_20', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end


% Compare No quadratic term, erfi_beta_20, and our method under model gap%%%

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(erfi_beta_20_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_beta_20_mod_4_time(plot_indices), erfi_beta_20_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_beta_20_mod_4_time(plot_indices), erfi_beta_20_mod_4_command_data(plot_indices, 1));
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

% 0.4m/s
if size(erfi_beta_50_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_beta_50_mod_4_time(plot_indices), erfi_beta_50_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_beta_50_mod_4_time(plot_indices), erfi_beta_50_mod_4_command_data(plot_indices, 1));
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

if size(erfi_beta_80_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_beta_80_mod_4_time(plot_indices), erfi_beta_80_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_beta_80_mod_4_time(plot_indices), erfi_beta_80_mod_4_command_data(plot_indices, 1));
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

if size(erfi_beta_20_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_beta_20_mod_6_time(plot_indices), erfi_beta_20_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_beta_20_mod_6_time(plot_indices), erfi_beta_20_mod_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
end


nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(erfi_beta_50_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_beta_50_mod_6_time(plot_indices), erfi_beta_50_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_beta_50_mod_6_time(plot_indices), erfi_beta_50_mod_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;
if size(erfi_beta_80_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_beta_80_mod_6_time(plot_indices), erfi_beta_80_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_beta_80_mod_6_time(plot_indices), erfi_beta_80_mod_6_command_data(plot_indices, 1));
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




% Compare No quadratic term, erfi_beta_20, and our method under environment gap%%%
nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(erfi_beta_20_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_beta_20_env_4_time(plot_indices), erfi_beta_20_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_beta_20_env_4_time(plot_indices), erfi_beta_20_env_4_command_data(plot_indices, 1));
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

% 0.4m/s
if size(erfi_beta_50_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_beta_50_env_4_time(plot_indices), erfi_beta_50_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_beta_50_env_4_time(plot_indices), erfi_beta_50_env_4_command_data(plot_indices, 1));
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

if size(erfi_beta_80_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_beta_80_env_4_time(plot_indices), erfi_beta_80_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_beta_80_env_4_time(plot_indices), erfi_beta_80_env_4_command_data(plot_indices, 1));
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

if size(erfi_beta_20_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_beta_20_env_6_time(plot_indices), erfi_beta_20_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_beta_20_env_6_time(plot_indices), erfi_beta_20_env_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(erfi_beta_50_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_beta_50_env_6_time(plot_indices), erfi_beta_50_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_beta_50_env_6_time(plot_indices), erfi_beta_50_env_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end


nexttile;
if size(erfi_beta_80_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(erfi_beta_80_env_6_time(plot_indices), erfi_beta_80_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(erfi_beta_80_env_6_time(plot_indices), erfi_beta_80_env_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end


han = axes(gcf, 'Visible', 'off'); % Create an invisible axes that spans the figure
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
xlabel(han, 'Time, [s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [0.5, -0.05, 0]); % Global x-axis label
ylabel(han, 'Velocity, [m/s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [-0.07, 0.5, 0]);
%%
%% Evaluate NQ per noise lim
figure('Name', 'NQ');
set(gcf, 'Units', 'inches', 'Position', [1, 1, double_column_eight_row_size]); % Choose single_column_size or double_column_size here
% sgtitle("Base Forward Velocity Tracking Performance in Nominal Settings", 'FontSize', sgtitle_fontsize);
t = tiledlayout(8, 3, 'TileSpacing', 'compact', 'Padding', 'normal'); % Adjust 'TileSpacing' and 'Padding' as needed
nexttile;
plot(nq_20_nom_4_time(plot_indices), nq_20_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(nq_20_nom_4_time(plot_indices), nq_20_nom_4_command_data(plot_indices, 1));
title('nq 20', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
nexttile;
plot(nq_50_nom_4_time(plot_indices), nq_50_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(nq_50_nom_4_time(plot_indices), nq_50_nom_4_command_data(plot_indices, 1));
title('nq 50', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');

nexttile;
plot(nq_80_nom_4_time(plot_indices), nq_80_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(nq_80_nom_4_time(plot_indices), nq_80_nom_4_command_data(plot_indices, 1));
title('nq 80', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
lgd.Position(1) = lgd.Position(1) + 0.02;
lgd.Position(2) = lgd.Position(2) - 0.01;

% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast')

nexttile;
plot(nq_20_nom_6_time(plot_indices), nq_20_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(nq_20_nom_6_time(plot_indices), nq_20_nom_6_command_data(plot_indices, 1));
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;

nexttile;
plot(nq_50_nom_6_time(plot_indices), nq_50_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(nq_50_nom_6_time(plot_indices), nq_50_nom_6_command_data(plot_indices, 1));
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;

% sgtitle("Base Forward Velocity Tracking Performance in Nominal Settings", 'FontSize', sgtitle_fontsize);
nexttile;
plot(nq_80_nom_6_time(plot_indices), nq_80_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(nq_80_nom_6_time(plot_indices), nq_80_nom_6_command_data(plot_indices, 1));
% title('No quadratic term', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'southeast')



% Compare No quadratic term, nq_20, and our method under actuator gap
nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(nq_20_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
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

% 0.4m/s
if size(nq_50_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
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

if size(nq_80_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(nq_80_act_4_time(plot_indices), nq_80_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nq_80_act_4_time(plot_indices), nq_80_act_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

% 0.6m/s
nexttile;
if size(nq_20_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
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

% 0.4m/s
if size(nq_50_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
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

if size(nq_80_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    % title('nq_20', 'FontSize', subtitle_fontsize);
else
    plot(nq_80_act_6_time(plot_indices), nq_80_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nq_80_act_6_time(plot_indices), nq_80_act_6_command_data(plot_indices, 1));
    % title('nq_20', 'FontSize', subtitle_fontsize);
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end


% Compare No quadratic term, nq_20, and our method under model gap%%%

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(nq_20_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
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

% 0.4m/s
if size(nq_50_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
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

if size(nq_80_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(nq_80_mod_4_time(plot_indices), nq_80_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nq_80_mod_4_time(plot_indices), nq_80_mod_4_command_data(plot_indices, 1));
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

if size(nq_20_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(nq_20_mod_6_time(plot_indices), nq_20_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nq_20_mod_6_time(plot_indices), nq_20_mod_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
end


nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(nq_50_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(nq_50_mod_6_time(plot_indices), nq_50_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nq_50_mod_6_time(plot_indices), nq_50_mod_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

nexttile;
if size(nq_80_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(nq_80_mod_6_time(plot_indices), nq_80_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nq_80_mod_6_time(plot_indices), nq_80_mod_6_command_data(plot_indices, 1));
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




% Compare No quadratic term, nq_20, and our method under environment gap%%%
nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(nq_20_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
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

% 0.4m/s
if size(nq_50_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
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

if size(nq_80_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(nq_80_env_4_time(plot_indices), nq_80_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nq_80_env_4_time(plot_indices), nq_80_env_4_command_data(plot_indices, 1));
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

if size(nq_20_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(nq_20_env_6_time(plot_indices), nq_60_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nq_20_env_6_time(plot_indices), nq_60_env_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(nq_50_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
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


nexttile;
if size(nq_80_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(nq_80_env_6_time(plot_indices), nq_80_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(nq_80_env_6_time(plot_indices), nq_80_env_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end


han = axes(gcf, 'Visible', 'off'); % Create an invisible axes that spans the figure
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
xlabel(han, 'Time, [s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [0.5, -0.05, 0]); % Global x-axis label
ylabel(han, 'Velocity, [m/s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [-0.07, 0.5, 0]);
%%
%% Evaluate Hybrid per noise lim
figure('Name', 'Hybrid');
set(gcf, 'Units', 'inches', 'Position', [1, 1, single_column_eight_row_size]); % Choose single_column_size or double_column_size here
% sgtitle("Base Forward Velocity Tracking Performance in Nominal Settings", 'FontSize', sgtitle_fontsize);
t = tiledlayout(8, 2, 'TileSpacing', 'compact', 'Padding', 'normal'); % Adjust 'TileSpacing' and 'Padding' as needed
nexttile;
plot(hy_50_nom_4_time(plot_indices), hy_50_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(hy_50_nom_4_time(plot_indices), hy_50_nom_4_command_data(plot_indices, 1));
title('hy 50', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
nexttile;
plot(hy_100_nom_4_time(plot_indices), hy_100_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(hy_100_nom_4_time(plot_indices), hy_100_nom_4_command_data(plot_indices, 1));
title('hy 100', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');

nexttile;
plot(hy_50_nom_6_time(plot_indices), hy_50_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(hy_50_nom_6_time(plot_indices), hy_50_nom_6_command_data(plot_indices, 1));
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;

nexttile;
plot(hy_100_nom_6_time(plot_indices), hy_100_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(hy_100_nom_6_time(plot_indices), hy_100_nom_6_command_data(plot_indices, 1));
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;


% Compare No quadratic term, nq_20, and our method under actuator gap
nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(hy_50_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(hy_50_act_4_time(plot_indices), hy_50_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(hy_50_act_4_time(plot_indices), hy_50_act_4_command_data(plot_indices, 1));
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

% 0.4m/s
if size(hy_100_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(hy_100_act_4_time(plot_indices), hy_100_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(hy_100_act_4_time(plot_indices), hy_100_act_4_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end

% 0.6m/s
nexttile;
if size(hy_50_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    % title('DR', 'FontSize', subtitle_fontsize);
else
    plot(hy_50_act_6_time(plot_indices), hy_50_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(hy_50_act_6_time(plot_indices), hy_50_act_6_command_data(plot_indices, 1));
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

% 0.4m/s
if size(hy_100_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(hy_100_act_6_time(plot_indices), hy_100_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(hy_100_act_6_time(plot_indices), hy_100_act_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end




% Compare No quadratic term, hy_50, and our method under model gap%%%

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(hy_50_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(hy_50_mod_4_time(plot_indices), hy_50_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(hy_50_mod_4_time(plot_indices), hy_50_mod_4_command_data(plot_indices, 1));
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

% 0.4m/s
if size(hy_100_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(hy_100_mod_4_time(plot_indices), hy_100_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(hy_100_mod_4_time(plot_indices), hy_100_mod_4_command_data(plot_indices, 1));
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

if size(hy_50_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(hy_50_mod_6_time(plot_indices), hy_50_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(hy_50_mod_6_time(plot_indices), hy_50_mod_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
end


nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(hy_100_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(hy_100_mod_6_time(plot_indices), hy_100_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(hy_100_mod_6_time(plot_indices), hy_100_mod_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end


% Compare No quadratic term, hy_50, and our method under environment gap%%%
nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(hy_50_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(hy_50_env_4_time(plot_indices), hy_50_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(hy_50_env_4_time(plot_indices), hy_50_env_4_command_data(plot_indices, 1));
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

% 0.4m/s
if size(hy_100_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(hy_100_env_4_time(plot_indices), hy_100_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(hy_100_env_4_time(plot_indices), hy_100_env_4_command_data(plot_indices, 1));
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

if size(hy_50_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(hy_50_env_6_time(plot_indices), hy_50_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(hy_50_env_6_time(plot_indices), hy_50_env_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
end

nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(hy_100_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
else
    plot(hy_100_env_6_time(plot_indices), hy_100_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(hy_100_env_6_time(plot_indices), hy_100_env_6_command_data(plot_indices, 1));
    % xlabel('Time, [s]', 'FontSize', label_fontsize);
    % ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
    xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
    ylim(ylimit); % Adjust these values to your desired y-axis range
    set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
    grid on;
    % lg = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');
    % lg.Position(1) = lg.Position(1) + 0.05;
end



han = axes(gcf, 'Visible', 'off'); % Create an invisible axes that spans the figure
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
xlabel(han, 'Time, [s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [0.5, -0.05, 0]); % Global x-axis label
ylabel(han, 'Velocity, [m/s]', 'FontSize', label_fontsize, 'Units', 'normalized', 'Position', [-0.07, 0.5, 0]);
%%
%% Compare Performance in Nominal Setting 
figure('Name', 'Nominal performance comparison, (DR, ERFI, Ours)');
set(gcf, 'Units', 'inches', 'Position', [1, 1, double_column_size]); % Choose single_column_size or double_column_size here
% sgtitle("Base Forward Velocity Tracking Performance in Nominal Settings", 'FontSize', sgtitle_fontsize);
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'normal'); % Adjust 'TileSpacing' and 'Padding' as needed
nexttile;
plot(DR_nom_4_time(plot_indices), DR_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(DR_nom_4_time(plot_indices), DR_nom_4_command_data(plot_indices, 1));
title('DR', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% lgd = legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast');

nexttile;
plot(ERFI_nom_4_time(plot_indices), ERFI_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(ERFI_nom_4_time(plot_indices), ERFI_nom_4_command_data(plot_indices, 1));
title('ERFI', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'northeast')

nexttile;
plot(OUR_nom_4_time(plot_indices), OUR_nom_4_q_dot_virtual(plot_indices, 1));
hold on;
plot(OUR_nom_4_time(plot_indices), OUR_nom_4_command_data(plot_indices, 1));
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
plot(DR_nom_6_time(plot_indices), DR_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(DR_nom_6_time(plot_indices), DR_nom_6_command_data(plot_indices, 1));
% title('DR', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'southeast')

nexttile;
plot(ERFI_nom_6_time(plot_indices), ERFI_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(ERFI_nom_6_time(plot_indices), ERFI_nom_6_command_data(plot_indices, 1));
% title('ERFI', 'FontSize', subtitle_fontsize);
% xlabel('Time, [s]', 'FontSize', label_fontsize);
% ylabel('Velocity, [m/s]', 'FontSize', label_fontsize);
xlim([plot_start_time, plot_end_time]); % Adjust these values to your desired x-axis range
ylim(ylimit); % Adjust these values to your desired y-axis range
set(gca, 'FontSize', tick_fontsize); % This sets the font size for axis tick labels
grid on;
% legend('Measured', 'Commanded', 'FontSize', legend_fontsize, 'Location', 'southeast')

nexttile;
plot(OUR_nom_6_time(plot_indices), OUR_nom_6_q_dot_virtual(plot_indices, 1));
hold on;
plot(OUR_nom_6_time(plot_indices), OUR_nom_6_command_data(plot_indices, 1));
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
%%
%% Compare DR, ERFI, and our method under actuator gap
figure('Name', 'Actuator gap, (DR, ERFI, Ours)');
set(gcf, 'Units', 'inches', 'Position', [1, 1, double_column_size]); % Choose single_column_size or double_column_size here

t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'normal'); % Adjust 'TileSpacing' and 'Padding' as needed
nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(DR_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    title('DR', 'FontSize', subtitle_fontsize);
else
    plot(DR_act_4_time(plot_indices), DR_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(DR_act_4_time(plot_indices), DR_act_4_command_data(plot_indices, 1));
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

if size(ERFI_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    title('ERFI', 'FontSize', subtitle_fontsize);
else
    plot(ERFI_act_4_time(plot_indices), ERFI_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(ERFI_act_4_time(plot_indices), ERFI_act_4_command_data(plot_indices, 1));
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

if size(OUR_act_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    title('Our method', 'FontSize', subtitle_fontsize);
else
    plot(OUR_act_4_time(plot_indices), OUR_act_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(OUR_act_4_time(plot_indices), OUR_act_4_command_data(plot_indices, 1));
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
if size(DR_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    % title('DR', 'FontSize', subtitle_fontsize);
else
    plot(DR_act_6_time(plot_indices), DR_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(DR_act_6_time(plot_indices), DR_act_6_command_data(plot_indices, 1));
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

if size(ERFI_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    % title('ERFI', 'FontSize', subtitle_fontsize);
else
    plot(ERFI_act_6_time(plot_indices), ERFI_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(ERFI_act_6_time(plot_indices), ERFI_act_6_command_data(plot_indices, 1));
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

if size(OUR_act_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    % title('Our method', 'FontSize', subtitle_fontsize);
else
    plot(OUR_act_6_time(plot_indices), OUR_act_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(OUR_act_6_time(plot_indices), OUR_act_6_command_data(plot_indices, 1));
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

%%
%% Compare DR, ERFI, and our method under model gap%%%
figure('Name', 'Model gap, (DR, ERFI, Ours)');
set(gcf, 'Units', 'inches', 'Position', [1, 1, double_column_size]); % Choose single_column_size or double_column_size here

t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'normal'); % Adjust 'TileSpacing' and 'Padding' as needed
nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(DR_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    title('DR', 'FontSize', subtitle_fontsize);
else
    plot(DR_mod_4_time(plot_indices), DR_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(DR_mod_4_time(plot_indices), DR_mod_4_command_data(plot_indices, 1));
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

if size(ERFI_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    title('ERFI', 'FontSize', subtitle_fontsize);
else
    plot(ERFI_mod_4_time(plot_indices), ERFI_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(ERFI_mod_4_time(plot_indices), ERFI_mod_4_command_data(plot_indices, 1));
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

if size(OUR_mod_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    title('Our method', 'FontSize', subtitle_fontsize);
else
    plot(OUR_mod_4_time(plot_indices), OUR_mod_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(OUR_mod_4_time(plot_indices), OUR_mod_4_command_data(plot_indices, 1));
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
if size(DR_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    % title('DR', 'FontSize', subtitle_fontsize);
else
    plot(DR_mod_6_time(plot_indices), DR_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(DR_mod_6_time(plot_indices), DR_mod_6_command_data(plot_indices, 1));
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

if size(ERFI_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    % title('ERFI', 'FontSize', subtitle_fontsize);
else
    plot(ERFI_mod_6_time(plot_indices), ERFI_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(ERFI_mod_6_time(plot_indices), ERFI_mod_6_command_data(plot_indices, 1));
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

if size(OUR_mod_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    % title('Our method', 'FontSize', subtitle_fontsize);
else
    plot(OUR_mod_6_time(plot_indices), OUR_mod_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(OUR_mod_6_time(plot_indices), OUR_mod_6_command_data(plot_indices, 1));
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
%%
%% Compare DR, ERFI, and our method under environment gap%%%
figure('Name', 'Environment gap, (DR, ERFI, Ours)');
set(gcf, 'Units', 'inches', 'Position', [1, 1, double_column_size]); % Choose single_column_size or double_column_size here

t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'normal'); % Adjust 'TileSpacing' and 'Padding' as needed
nexttile;% sgtitle('Base Forward Velocity Tracking Performance under Actuator Gap', 'FontSize', sgtitle_fontsize)

% 0.4m/s
if size(DR_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    title('DR', 'FontSize', subtitle_fontsize);
else
    plot(DR_env_4_time(plot_indices), DR_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(DR_env_4_time(plot_indices), DR_env_4_command_data(plot_indices, 1));
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

if size(ERFI_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    title('ERFI', 'FontSize', subtitle_fontsize);
else
    plot(ERFI_env_4_time(plot_indices), ERFI_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(ERFI_env_4_time(plot_indices), ERFI_env_4_command_data(plot_indices, 1));
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

if size(OUR_env_4_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    title('Our method', 'FontSize', subtitle_fontsize);
else
    plot(OUR_env_4_time(plot_indices), OUR_env_4_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(OUR_env_4_time(plot_indices), OUR_env_4_command_data(plot_indices, 1));
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
if size(DR_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    % title('DR', 'FontSize', subtitle_fontsize);
else
    plot(DR_env_6_time(plot_indices), DR_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(DR_env_6_time(plot_indices), DR_env_6_command_data(plot_indices, 1));
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

if size(ERFI_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    % title('ERFI', 'FontSize', subtitle_fontsize);
else
    plot(ERFI_env_6_time(plot_indices), ERFI_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(ERFI_env_6_time(plot_indices), ERFI_env_6_command_data(plot_indices, 1));
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

if size(OUR_env_6_q_dot_virtual, 1) < plot_end_time * control_hz
    text(0.5, 0.5, 'FAILED', 'HorizontalAlignment', 'center', 'FontSize', subtitle_fontsize);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', tick_fontsize); % Hide axes
    % title('Our method', 'FontSize', subtitle_fontsize);
else
    plot(OUR_env_6_time(plot_indices), OUR_env_6_q_dot_virtual(plot_indices, 1));
    hold on;
    plot(OUR_env_6_time(plot_indices), OUR_env_6_command_data(plot_indices, 1));
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
%%