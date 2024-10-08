#include "tocabi_lib/robot_data.h"
#include "wholebody_functions.h"
#include <random>
#include <cmath>

#include <ros/ros.h>
#include <sensor_msgs/Joy.h>
#include <tocabi_msgs/RLCommand.h>

class CustomController
{
public:
    CustomController(RobotData &rd);
    Eigen::VectorQd getControl();

    //void taskCommandToCC(TaskCommand tc_);
    
    void computeSlow();
    void computeFast();
    void computePlanner();
    void copyRobotData(RobotData &rd_l);

    RobotData &rd_;
    RobotData rd_cc_;

    void loadNetwork();
    void processNoise();
    void processBias();
    void initBias();
    Eigen::Matrix<double, MODEL_DOF, 1> q_bias_;
    void processObservation();
    void feedforwardPolicy();
    MatrixXd conv_layer(const MatrixXd &input, const MatrixXd &weights, const MatrixXd &biases, int kernel_size, int stride);
    void feedforwardEncoder();
    void initVariable();

    Eigen::Vector3d mat2euler(Eigen::Matrix3d mat);

    static const int num_action = 12;
    static const int num_actuator_action = 12;
    static const int num_cur_state = 48;
    static const int num_cur_internal_state = 36;
    static const int num_state_skip = 2;
    static const int num_state_hist = 10;
    static const int num_state = num_cur_state * num_state_hist; // num_cur_internal_state*num_state_hist+num_action*(num_state_hist-1);
    static const int num_hidden1 = 512;
    static const int num_hidden2 = 512;
    static const int num_hidden3 = 128;

    Eigen::MatrixXd policy_net_w0_;
    Eigen::MatrixXd policy_net_b0_;
    Eigen::MatrixXd policy_net_w2_;
    Eigen::MatrixXd policy_net_b2_;
    Eigen::MatrixXd policy_net_w4_;
    Eigen::MatrixXd policy_net_b4_;
    Eigen::MatrixXd action_net_w_;
    Eigen::MatrixXd action_net_b_;
    Eigen::MatrixXd hidden_layer1_;
    Eigen::MatrixXd hidden_layer2_;
    Eigen::MatrixXd hidden_layer3_;
    Eigen::MatrixXd rl_action_;

    Eigen::MatrixXd value_net_w0_;
    Eigen::MatrixXd value_net_b0_;
    Eigen::MatrixXd value_net_w2_;
    Eigen::MatrixXd value_net_b2_;
    Eigen::MatrixXd value_net_w4_;
    Eigen::MatrixXd value_net_b4_;
    Eigen::MatrixXd value_net_w_;
    Eigen::MatrixXd value_net_b_;
    Eigen::MatrixXd value_hidden_layer1_;
    Eigen::MatrixXd value_hidden_layer2_;
    Eigen::MatrixXd value_hidden_layer3_;
    double value_;

    Eigen::MatrixXd morph_net_w0_;
    Eigen::MatrixXd morph_net_b0_;
    Eigen::MatrixXd morph_net_w2_;
    Eigen::MatrixXd morph_net_b2_;
    Eigen::MatrixXd morph_net_w_;
    Eigen::MatrixXd morph_net_b_;
    Eigen::MatrixXd morph_hidden_layer1_;
    Eigen::MatrixXd morph_hidden_layer2_;
    Eigen::MatrixXd morphnet_input_;
    Eigen::MatrixXd morphnet_output_;    
    static const int morph_num_hidden_1 = 512;
    static const int morph_num_hidden_2 = 256;   
    static const int morph_params_dim = 4;
    static const int morph_history_len_ = 5;
    static const int morph_history_skip_ = 1;     
    static const bool morphnet = true;
    void loadMorphnet();
    void feedforwardMorphnet();

    bool stop_by_value_thres_ = false;
    Eigen::Matrix<double, MODEL_DOF, 1> q_stop_;
    float stop_start_time_;
    
    Eigen::MatrixXd state_;
    Eigen::MatrixXd state_cur_;
    Eigen::MatrixXd state_buffer_;
    Eigen::MatrixXd state_mean_;
    Eigen::MatrixXd state_var_;

    // encoder
    static const bool use_encoder_ = true;
    bool encoder_initialized = false;
    
    Eigen::MatrixXd state_history_;
    Eigen::MatrixXd encoder_input_;
    Eigen::MatrixXd encoder_output_;
    static const int history_len_ = 50;
    static const int history_skip_ = 5;
    static const int encoder_dim_ = 8;
    static const int encoder_conv1_kernel_ = 6;
    static const int encoder_conv1_stride_ = 3;
    static const int encoder_conv1_output_channel_ = 32;
    static const int encoder_conv1_input_channel_ = num_cur_state;
    static const int encoder_conv1_output_size_ = (history_len_-encoder_conv1_kernel_) / encoder_conv1_stride_ + 1;

    static const int encoder_conv2_kernel_ = 4;
    static const int encoder_conv2_stride_ = 2;
    static const int encoder_conv2_output_channel_ = 16;
    static const int encoder_conv2_input_channel_ = encoder_conv1_output_channel_;
    static const int encoder_conv2_output_size_ = (encoder_conv1_output_size_-encoder_conv2_kernel_) / encoder_conv2_stride_ + 1;

    static const int encoder_fc_input_ = 16*encoder_conv2_output_size_;

    Eigen::MatrixXd encoder_conv1_w_;
    Eigen::MatrixXd encoder_conv1_b_;
    Eigen::MatrixXd encoder_hidden_layer1_;
    Eigen::MatrixXd encoder_conv2_w_;
    Eigen::MatrixXd encoder_conv2_b_;
    Eigen::MatrixXd encoder_hidden_layer2_;
    Eigen::MatrixXd encoder_fc_w_;
    Eigen::MatrixXd encoder_fc_b_;
    void loadEncoderNetwork();

    Eigen::MatrixXd policy_input_;
    static const int policy_input_dim_ = (morphnet) 
                                     ? num_state + morph_params_dim
                                     : (use_encoder_) 
                                       ? num_state + encoder_dim_
                                       : num_state;
    std::ofstream writeFile;
    std::ofstream actuator_data_file;
    bool actuator_net_log = true;

    float phase_ = 0.0;

    bool is_on_robot_ = false;
    bool is_write_file_ = true;
    Eigen::Matrix<double, MODEL_DOF, 1> q_dot_lpf_;

    Eigen::Matrix<double, MODEL_DOF, 1> q_init_;
    Eigen::Matrix<double, MODEL_DOF, 1> q_noise_;
    Eigen::Matrix<double, MODEL_DOF, 1> q_noise_pre_;
    Eigen::Matrix<double, MODEL_DOF, 1> q_vel_noise_;

    Eigen::Matrix<double, MODEL_DOF, 1> torque_init_;
    Eigen::Matrix<double, MODEL_DOF, 1> torque_spline_;
    Eigen::Matrix<double, MODEL_DOF, 1> torque_rl_;
    Eigen::Matrix<double, MODEL_DOF, 1> torque_bound_;

    Eigen::Matrix<double, MODEL_DOF, MODEL_DOF> kp_;
    Eigen::Matrix<double, MODEL_DOF, MODEL_DOF> kv_;

    float start_time_;
    float time_inference_pre_ = 0.0;
    float time_write_pre_ = 0.0;

    double time_cur_;
    double time_pre_;
    double action_dt_accumulate_ = 0.0;

    Eigen::Vector3d euler_angle_;

    // Joystick
    ros::NodeHandle nh_;

    void joyCallback(const sensor_msgs::Joy::ConstPtr& joy);
    ros::Subscriber joy_sub_;

    // slider command
    void rlcommandCallback(const tocabi_msgs::RLCommand::ConstPtr& command);
    ros::Subscriber rl_command_sub_;


    double target_vel_x_ = 0.0;
    double pre_target_vel_x_ = 0.0;
    double target_vel_y_ = 0.0;
    double pre_target_vel_y_ = 0.0;
    double target_heading_ = 0.0;
    double pre_target_vel_yaw_ = 0.0;
    double heading = 0.0;

private:
    Eigen::VectorQd ControlVal_;
};