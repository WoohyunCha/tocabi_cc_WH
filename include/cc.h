#include "tocabi_lib/robot_data.h"
#include "wholebody_functions.h"
#include <random>
#include <cmath>

#include <ros/ros.h>
#include <sensor_msgs/Joy.h>

class CustomController
{
public:
    CustomController(RobotData &rd);
    Eigen::VectorQd getControl();

    //void taskCommandToCC(TaskCommand tc_);

    const double hz_ =125.;
    const double pd_hz_ = 2000;
    double del_t = 1 / hz_;


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
    static const int num_cur_state = 73;
    static const int num_cur_internal_state = 73;
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


    bool stop_by_value_thres_ = false;
    Eigen::Matrix<double, MODEL_DOF, 1> q_stop_;
    float stop_start_time_;
    
    Eigen::MatrixXd state_;
    Eigen::MatrixXd state_cur_;
    Eigen::MatrixXd state_buffer_;
    Eigen::MatrixXd state_mean_;
    Eigen::MatrixXd state_var_;

    // encoder
    static const bool use_encoder_ = false;
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
    static const int policy_input_dim_ =  (use_encoder_) 
                                       ? num_state + encoder_dim_
                                       : num_state;
    std::ofstream writeFile;
    std::ofstream actuator_data_file;
    bool actuator_net_log = false;

    float phase_ = 0.0;

    bool is_on_robot_ = true;
    bool is_write_file_ = true;
    Eigen::Matrix<double, MODEL_DOF, 1> q_dot_lpf_;

    Eigen::Matrix<double, MODEL_DOF, 1> q_init_;
    Eigen::Matrix<double, MODEL_DOF, 1> q_noise_;
    Eigen::Matrix<double, MODEL_DOF, 1> q_noise_pre_;
    Eigen::Matrix<double, MODEL_DOF, 1> q_vel_noise_;
    Eigen::Vector12d q_leg_desired_;

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
    Vector3_t base_lin_vel, base_ang_vel;

    double target_vel_x_ = 0.0;
    double pre_target_vel_x_ = 0.0;
    double target_vel_y_ = 0.0;
    double pre_target_vel_y_ = 0.0;
    double target_heading_ = 0.0;
    double pre_target_vel_yaw_ = 0.0;
    double heading = 0.0;

    std::string base_path = "";
    std::string loadPathFromConfig(const std::string &config_file);
    std::string loadCommand(const std::string &command_file);


    // BIPED WALKING PARAMETER
    void walkingParameterSetting();
    const int number_of_foot_step = 3;
    Eigen::Vector3d phase_indicator_;
    Eigen::Vector3d step_length_x_;
    Eigen::Vector3d step_length_y_;
    Eigen::Vector3d step_yaw_;
    Eigen::Vector3d t_dsp_;
    Eigen::Vector3d t_ssp_;
    Eigen::Vector3d t_dsp_seconds;
    Eigen::Vector3d t_ssp_seconds;
    Eigen::Vector3d foot_height_;
    Eigen::Vector3d t_total_;
    int first_stance_foot_ = 0; // 1 means right foot stance, 0 means left foot stance
    const double com_height_ = 0.68;
    int current_step_num = 0;

    double t_last_;
    double t_start_;
    double t_temp_;  
    double t_double1_;
    double t_double2_;
    double zmp_offset = 0.;




    // Policy input
    Eigen::VectorXd target_com_state_stance_frame_;
    Eigen::VectorXd target_swing_state_stance_frame_;


    Eigen::MatrixXd foot_step_;
    Eigen::MatrixXd foot_step_support_frame_;
    Eigen::MatrixXd foot_step_support_frame_offset_;

    Eigen::Isometry3d pelv_support_start_;
    // Eigen::Isometry3d pelv_support_init_;
    Eigen::Isometry3d pelv_support_init_yaw_;

    Eigen::Vector3d com_desired_;
    Eigen::Vector3d com_desired_dot_;

    Eigen::Vector3d com_support_init_yaw_;
    Eigen::Vector3d com_support_init_dot_yaw_;

    Eigen::Vector3d com_support_current_;
    Eigen::Vector3d com_support_current_dot_;
    Eigen::Vector3d com_support_current_dot_prev_;
    Eigen::Vector3d com_global_current_;
    Eigen::Vector3d com_global_current_dot_;

    Eigen::Vector3d pelv_rpy_current_;
    Eigen::Vector3d rfoot_rpy_current_;
    Eigen::Vector3d lfoot_rpy_current_;
    Eigen::Isometry3d pelv_yaw_rot_current_from_global_;
    Eigen::Isometry3d rfoot_roll_rot_;
    Eigen::Isometry3d lfoot_roll_rot_;
    Eigen::Isometry3d rfoot_pitch_rot_;
    Eigen::Isometry3d lfoot_pitch_rot_;
    Eigen::Isometry3d rfoot_yaw_rot_;
    Eigen::Isometry3d lfoot_yaw_rot_;

    Eigen::Isometry3d pelv_global_current_;
    Eigen::Isometry3d lfoot_global_current_;
    Eigen::Isometry3d rfoot_global_current_;
    Eigen::Isometry3d pelv_global_init_;
    Eigen::Isometry3d lfoot_global_init_;
    Eigen::Isometry3d rfoot_global_init_;

    Eigen::Isometry3d pelv_trajectory_support_; //local frame
    Eigen::Isometry3d pelv_trajectory_support_fast_; //local frame
    Eigen::Isometry3d pelv_trajectory_support_slow_; //local frame
    
    Eigen::Isometry3d rfoot_trajectory_support_;  //local frame
    Eigen::Isometry3d lfoot_trajectory_support_;
    Eigen::Isometry3d lfoot_trajectory_support_fast_;
    Eigen::Isometry3d lfoot_trajectory_support_slow_;

    Eigen::Vector3d rfoot_trajectory_euler_support_;
    Eigen::Vector3d lfoot_trajectory_euler_support_;

    Eigen::Isometry3d pelv_trajectory_float_; //pelvis frame

    Eigen::Vector3d com_trajectory_float_;

    Eigen::Isometry3d lfoot_trajectory_float_;
    Eigen::Isometry3d lfoot_trajectory_float_fast_;
    Eigen::Isometry3d lfoot_trajectory_float_slow_;

    Eigen::Isometry3d rfoot_trajectory_float_;
    Eigen::Isometry3d rfoot_trajectory_float_fast_;
    Eigen::Isometry3d rfoot_trajectory_float_slow_;

    Eigen::Vector3d pelv_support_euler_init_;
    Eigen::Vector3d pelv_support_euler_init_yaw_;
    Eigen::Vector3d lfoot_support_euler_init_;
    Eigen::Vector3d lfoot_support_euler_init_yaw_;
    Eigen::Vector3d rfoot_support_euler_init_;
    Eigen::Vector3d rfoot_support_euler_init_yaw_;
    double wn = sqrt(GRAVITY / com_height_);

    double walking_end_flag = 0;
    
    Eigen::Isometry3d swingfoot_global_current_; 
    Eigen::Isometry3d supportfoot_global_current_; // Only yaw is considered

    Eigen::Isometry3d pelv_support_current_;
    Eigen::Isometry3d lfoot_support_current_;
    Eigen::Isometry3d rfoot_support_current_;

    Eigen::Isometry3d lfoot_support_init_;
    Eigen::Isometry3d lfoot_support_init_yaw_;
    Eigen::Isometry3d rfoot_support_init_;
    Eigen::Isometry3d rfoot_support_init_yaw_;
    
    Eigen::Isometry3d target_com_state_float_frame_, target_lfoot_state_float_frame_, target_rfoot_state_float_frame_;
    
    Eigen::MatrixXd ref_zmp_;
    Eigen::MatrixXd ref_zmp_container;
    Eigen::MatrixXd ref_zmp_thread3;
    Eigen::VectorXd ref_com_yaw_;
    Eigen::VectorXd ref_com_yawvel_;
    // PREVIEW CONTROL
    Eigen::Vector3d x_preview_;
    Eigen::Vector3d y_preview_;

    Eigen::Vector3d xs_preview_;
    Eigen::Vector3d ys_preview_;
    Eigen::Vector3d xd_preview_;
    Eigen::Vector3d yd_preview_; 

    Eigen::MatrixXd Gi_preview_;
    Eigen::MatrixXd Gx_preview_;
    Eigen::VectorXd Gd_preview_;
    Eigen::MatrixXd A_preview_;
    Eigen::VectorXd B_preview_;
    Eigen::MatrixXd C_preview_;
    double UX_preview_, UY_preview_, EX_preview_, EY_preview_;


    Eigen::Vector6d l_ft_;
    Eigen::Vector6d r_ft_;
    Eigen::Vector6d l_ft_LPF;
    Eigen::Vector6d r_ft_LPF;
    Eigen::Vector2d zmp_measured_mj_;
    Eigen::Vector2d zmp_measured_LPF_;

    double P_angle_i = 0;
    double P_angle = 0;
    double P_angle_input_dot = 0;
    double P_angle_input = 0;
    double R_angle = 0;
    double R_angle_input_dot = 0;
    double R_angle_input = 0;
    double aa = 0; 
    double Y_angle_input = 0;

    // BOOLEAN OPERATOR
    bool walking_enable_;
    bool is_lfoot_support = false;
    bool is_rfoot_support = false;
    bool is_dsp1 = false;
    bool is_ssp  = false;
    bool is_dsp2 = false;
    bool is_preview_ctrl_init = true;

    // Logging
    Eigen::Isometry3d supportfoot_global_init_; // Used just to compute init_yaw_
    Eigen::Isometry3d supportfoot_global_init_yaw_;
    Eigen::VectorXd swing_state_stance_frame_;
    Eigen::VectorXd com_state_stance_frame_;

    // USER COMMAND
    double Lcommand_step_length_x_ = 0.1;
    double Lcommand_step_length_y_ = 0.22;
    double Lcommand_step_yaw_ = 0.;
    double Lcommand_t_dsp_ = 0.1;
    double Lcommand_t_ssp_ = 0.8;
    double Lcommand_foot_height_ = 0.08;

    double Rcommand_step_length_x_ = 0.1;
    double Rcommand_step_length_y_ = 0.22;
    double Rcommand_step_yaw_ = 0.;
    double Rcommand_t_dsp_ = 0.1;
    double Rcommand_t_ssp_ = 0.8;
    double Rcommand_foot_height_ = 0.08;

    // PREVIEW CONTROL
    void updateInitialState();
    void updateFootstepCommand();
    void getRobotState();
    void walkingStateMachine();
    void calculateFootStepTotal();
    void supportToFloatPattern();
    void floatToSupportFootstep();
    void updateNextStepTime();
    void resetPreviewState();
    void computeIkControl(const Eigen::Isometry3d &float_trunk_transform, const Eigen::Isometry3d &float_lleg_transform, const Eigen::Isometry3d &float_rleg_transform, Eigen::Vector12d &q_des);

    void getZmpTrajectory();
    void addZmpOffset();
    void zmpGenerator(const unsigned int norm_size);
    void onestepZmp(unsigned int current_step_number, Eigen::VectorXd &temp_px, Eigen::VectorXd &temp_py, Eigen::VectorXd& temp_yaw, Eigen::VectorXd &temp_yawvel);

    void getComTrajectory(); 
    void previewcontroller(double dt, int NL, int tick, 
                           Eigen::Vector3d &x_k, Eigen::Vector3d &y_k, double &UX, double &UY,
                           const Eigen::MatrixXd &Gi, const Eigen::VectorXd &Gd, const Eigen::MatrixXd &Gx, 
                           const Eigen::MatrixXd &A,  const Eigen::VectorXd &B,  const Eigen::MatrixXd &C);
    void preview_Parameter(double dt, int NL, Eigen::MatrixXd& Gi, Eigen::VectorXd& Gd, Eigen::MatrixXd& Gx, Eigen::MatrixXd& A, Eigen::VectorXd& B, Eigen::MatrixXd& C);
    void getComTrajectory_mpc();
    void getFootTrajectory(); 
    void getTargetState(); 

private:
    Eigen::VectorQd ControlVal_;
    unsigned int walking_tick = 0;
    unsigned int walking_tick_container = 0;

};