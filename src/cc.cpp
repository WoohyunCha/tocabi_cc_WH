#include "cc.h"

using namespace TOCABI;

CustomController::CustomController(RobotData &rd) : rd_(rd) //, wbc_(dc.wbc_)
{
    ControlVal_.setZero();

    if (is_write_file_)
    {
        if (is_on_robot_)
        {
            writeFile.open("/home/dyros/catkin_ws/src/tocabi_cc/result/data.csv", std::ofstream::out | std::ofstream::app);
            actuator_net_log = false;
        }
        else
        {
            writeFile.open("/home/cha/catkin_ws/src/tocabi_cc/result/data.csv", std::ofstream::out);
            actuator_data_file.open("/home/cha/catkin_ws/src/tocabi_cc/result/actuator_data.csv", std::ofstream::app);
        }
        writeFile << std::fixed << std::setprecision(8);
        actuator_data_file << std::fixed << std::setprecision(8);
        if (actuator_net_log) actuator_data_file << "START OF EPISODE\n";  // Mark the end of the data stream
    }
    initVariable();
    loadNetwork();

    joy_sub_ = nh_.subscribe<sensor_msgs::Joy>("joy", 10, &CustomController::joyCallback, this);
    rl_command_sub_ = nh_.subscribe<tocabi_msgs::RLCommand>("/tocabi/rlcommand", 10, &CustomController::rlcommandCallback, this);
}

Eigen::VectorQd CustomController::getControl()
{
    return ControlVal_;
}

void CustomController::loadNetwork()
{
    state_.setZero();
    rl_action_.setZero();


    string cur_path = "/home/cha/catkin_ws/src/tocabi_cc/";

    if (is_on_robot_)
    {
        cur_path = "/home/dyros/catkin_ws/src/tocabi_cc/";
    }
    std::ifstream file[14];
    // file[0].open(cur_path+"weight/a2c_network_actor_mlp_0_weight.txt", std::ios::in);
    // file[1].open(cur_path+"weight/a2c_network_actor_mlp_0_bias.txt", std::ios::in);
    // file[2].open(cur_path+"weight/a2c_network_actor_mlp_2_weight.txt", std::ios::in);
    // file[3].open(cur_path+"weight/a2c_network_actor_mlp_2_bias.txt", std::ios::in);
    // file[4].open(cur_path+"weight/a2c_network_mu_weight.txt", std::ios::in);
    // file[5].open(cur_path+"weight/a2c_network_mu_bias.txt", std::ios::in);
    // file[6].open(cur_path+"weight/obs_mean_fixed.txt", std::ios::in);
    // file[7].open(cur_path+"weight/obs_variance_fixed.txt", std::ios::in);
    // file[8].open(cur_path+"weight/a2c_network_critic_mlp_0_weight.txt", std::ios::in);
    // file[9].open(cur_path+"weight/a2c_network_critic_mlp_0_bias.txt", std::ios::in);
    // file[10].open(cur_path+"weight/a2c_network_critic_mlp_2_weight.txt", std::ios::in);
    // file[11].open(cur_path+"weight/a2c_network_critic_mlp_2_bias.txt", std::ios::in);
    // file[12].open(cur_path+"weight/a2c_network_value_weight.txt", std::ios::in);
    // file[13].open(cur_path+"weight/a2c_network_value_bias.txt", std::ios::in);
    file[0].open("/home/cha/isaac_ws/AMP_for_hardware/logs/policy/0_weight.txt", std::ios::in);
    file[1].open("/home/cha/isaac_ws/AMP_for_hardware/logs/policy/0_bias.txt", std::ios::in);
    file[2].open("/home/cha/isaac_ws/AMP_for_hardware/logs/policy/2_weight.txt", std::ios::in);
    file[3].open("/home/cha/isaac_ws/AMP_for_hardware/logs/policy/2_bias.txt", std::ios::in);
    file[4].open("/home/cha/isaac_ws/AMP_for_hardware/logs/policy/4_weight.txt", std::ios::in);
    file[5].open("/home/cha/isaac_ws/AMP_for_hardware/logs/policy/4_bias.txt", std::ios::in);
    file[6].open("/home/cha/isaac_ws/AMP_for_hardware/logs/normalizer/running_mean.txt", std::ios::in);
    file[7].open("/home/cha/isaac_ws/AMP_for_hardware/logs/normalizer/running_var.txt", std::ios::in);
    file[8].open("/home/cha/isaac_ws/AMP_for_hardware/logs/critic/0_weight.txt", std::ios::in);
    file[9].open("/home/cha/isaac_ws/AMP_for_hardware/logs/critic/0_bias.txt", std::ios::in);
    file[10].open("/home/cha/isaac_ws/AMP_for_hardware/logs/critic/2_weight.txt", std::ios::in);
    file[11].open("/home/cha/isaac_ws/AMP_for_hardware/logs/critic/2_bias.txt", std::ios::in);
    file[12].open("/home/cha/isaac_ws/AMP_for_hardware/logs/critic/4_weight.txt", std::ios::in);
    file[13].open("/home/cha/isaac_ws/AMP_for_hardware/logs/critic/4_bias.txt", std::ios::in);

    if(!file[0].is_open())
    {
        std::cout<<"Can not find the weight file"<<std::endl;
    }

    float temp;
    int row = 0;
    int col = 0;

    while(!file[0].eof() && row != policy_net_w0_.rows())
    {
        file[0] >> temp;
        if(temp != '\n')
        {
            policy_net_w0_(row, col) = temp;
            col ++;
            if (col == policy_net_w0_.cols())
            {
                col = 0;
                row ++;
            }
        }
    }
    row = 0;
    col = 0;
    while(!file[1].eof() && row != policy_net_b0_.rows())
    {
        file[1] >> temp;
        if(temp != '\n')
        {
            policy_net_b0_(row, col) = temp;
            col ++;
            if (col == policy_net_b0_.cols())
            {
                col = 0;
                row ++;
            }
        }
    }
    row = 0;
    col = 0;
    while(!file[2].eof() && row != policy_net_w2_.rows())
    {
        file[2] >> temp;
        if(temp != '\n')
        {
            policy_net_w2_(row, col) = temp;
            col ++;
            if (col == policy_net_w2_.cols())
            {
                col = 0;
                row ++;
            }
        }
    }
    row = 0;
    col = 0;
    while(!file[3].eof() && row != policy_net_b2_.rows())
    {
        file[3] >> temp;
        if(temp != '\n')
        {
            policy_net_b2_(row, col) = temp;
            col ++;
            if (col == policy_net_b2_.cols())
            {
                col = 0;
                row ++;
            }
        }
    }
    row = 0;
    col = 0;
    while(!file[4].eof() && row != action_net_w_.rows())
    {
        file[4] >> temp;
        if(temp != '\n')
        {
            action_net_w_(row, col) = temp;
            col ++;
            if (col == action_net_w_.cols())
            {
                col = 0;
                row ++;
            }
        }
    }
    row = 0;
    col = 0;
    while(!file[5].eof() && row != action_net_b_.rows())
    {
        file[5] >> temp;
        if(temp != '\n')
        {
            action_net_b_(row, col) = temp;
            col ++;
            if (col == action_net_b_.cols())
            {
                col = 0;
                row ++;
            }
        }
    }
    row = 0;
    col = 0;
    while(!file[6].eof() && row != state_mean_.rows())
    {
        file[6] >> temp;
        if(temp != '\n')
        {
            state_mean_(row, col) = temp;
            col ++;
            if (col == state_mean_.cols())
            {
                col = 0;
                row ++;
            }
        }
    }
    row = 0;
    col = 0;
    while(!file[7].eof() && row != state_var_.rows())
    {
        file[7] >> temp;
        if(temp != '\n')
        {
            state_var_(row, col) = temp + 1.e-4;
            col ++;
            if (col == state_var_.cols())
            {
                col = 0;
                row ++;
            }
        }
    }
    row = 0;
    col = 0;
    while(!file[8].eof() && row != value_net_w0_.rows())
    {
        file[8] >> temp;
        if(temp != '\n')
        {
            value_net_w0_(row, col) = temp;
            col ++;
            if (col == value_net_w0_.cols())
            {
                col = 0;
                row ++;
            }
        }
    }
    row = 0;
    col = 0;
    while(!file[9].eof() && row != value_net_b0_.rows())
    {
        file[9] >> temp;
        if(temp != '\n')
        {
            value_net_b0_(row, col) = temp;
            col ++;
            if (col == value_net_b0_.cols())
            {
                col = 0;
                row ++;
            }
        }
    }
    row = 0;
    col = 0;
    while(!file[10].eof() && row != value_net_w2_.rows())
    {
        file[10] >> temp;
        if(temp != '\n')
        {
            value_net_w2_(row, col) = temp;
            col ++;
            if (col == value_net_w2_.cols())
            {
                col = 0;
                row ++;
            }
        }
    }
    row = 0;
    col = 0;
    while(!file[11].eof() && row != value_net_b2_.rows())
    {
        file[11] >> temp;
        if(temp != '\n')
        {
            value_net_b2_(row, col) = temp;
            col ++;
            if (col == value_net_b2_.cols())
            {
                col = 0;
                row ++;
            }
        }
    }
    row = 0;
    col = 0;
    while(!file[12].eof() && row != value_net_w_.rows())
    {
        file[12] >> temp;
        if(temp != '\n')
        {
            value_net_w_(row, col) = temp;
            col ++;
            if (col == value_net_w_.cols())
            {
                col = 0;
                row ++;
            }
        }
    }
    row = 0;
    col = 0;
    while(!file[13].eof() && row != value_net_b_.rows())
    {
        file[13] >> temp;
        if(temp != '\n')
        {
            value_net_b_(row, col) = temp;
            col ++;
            if (col == value_net_b_.cols())
            {
                col = 0;
                row ++;
            }
        }
    }
    if (morphnet) loadMorphnet();
    else if (use_encoder_) loadEncoderNetwork();
}

void CustomController::initVariable()
{    
    
    policy_net_w0_.resize(num_hidden1, policy_input_dim_);
    policy_net_b0_.resize(num_hidden1, 1);
    policy_net_w2_.resize(num_hidden2, num_hidden1);
    policy_net_b2_.resize(num_hidden2, 1);
    action_net_w_.resize(num_action, num_hidden2);
    action_net_b_.resize(num_action, 1);
    hidden_layer1_.resize(num_hidden1, 1);
    hidden_layer2_.resize(num_hidden2, 1);
    rl_action_.resize(num_action, 1);

    value_net_w0_.resize(num_hidden1, policy_input_dim_);
    value_net_b0_.resize(num_hidden1, 1);
    value_net_w2_.resize(num_hidden2, num_hidden1);
    value_net_b2_.resize(num_hidden2, 1);
    value_net_w_.resize(1, num_hidden2);
    value_net_b_.resize(1, 1);
    value_hidden_layer1_.resize(num_hidden1, 1);
    value_hidden_layer2_.resize(num_hidden2, 1);
    
    // state_cur_.resize(num_cur_state, 1);
    state_cur_ = MatrixXd::Zero(num_cur_state, 1);
    state_.resize(num_state, 1);
    state_buffer_.resize(num_cur_state*num_state_skip*num_state_hist, 1);
    state_mean_.resize(num_cur_state, 1);
    state_var_.resize(num_cur_state, 1);

    // Encoder
    if (morphnet){
        state_history_.resize(num_cur_state, morph_history_len_ * morph_history_skip_);
        morphnet_input_.resize(num_cur_state * morph_history_len_, 1);
        morph_net_w0_.resize( morph_num_hidden_1, num_cur_state * morph_history_len_);
        morph_net_b0_.resize(morph_num_hidden_1, 1);

        morph_net_w2_.resize( morph_num_hidden_2, morph_num_hidden_1);
        morph_net_b2_.resize(morph_num_hidden_2, 1);

        morph_net_w_.resize(morph_params_dim, morph_num_hidden_2);
        morph_net_b_.resize(morph_params_dim, 1);

        morph_hidden_layer1_.resize(morph_num_hidden_1, 1);
        morph_hidden_layer2_.resize(morph_num_hidden_2, 1);

        morphnet_output_ = MatrixXd::Zero(morph_params_dim, 1);

    }
    else{
        state_history_.resize(num_cur_state, history_len_ * history_skip_);

        encoder_input_.resize(num_cur_state, history_len_);
        encoder_conv1_w_.resize(encoder_conv1_output_channel_, encoder_conv1_input_channel_*encoder_conv1_kernel_);
        encoder_conv1_b_.resize(encoder_conv1_output_channel_, 1);

        encoder_conv2_w_.resize(encoder_conv2_output_channel_, encoder_conv2_input_channel_*encoder_conv2_kernel_);
        encoder_conv2_b_.resize(encoder_conv2_output_channel_, 1);

        encoder_fc_w_.resize(encoder_dim_, encoder_conv2_output_channel_*encoder_conv2_output_size_);
        encoder_fc_b_.resize(encoder_dim_, 1);

        encoder_hidden_layer1_ = MatrixXd::Zero(encoder_conv1_output_channel_, encoder_conv1_output_size_);
        encoder_hidden_layer2_ = MatrixXd::Zero(encoder_conv2_output_channel_, encoder_conv2_output_size_);
        encoder_output_ = MatrixXd::Zero(encoder_dim_, 1);
    }

    if (morphnet) policy_input_.resize(num_state + morph_params_dim, 1);
    else if (use_encoder_) policy_input_.resize(num_state + encoder_dim_, 1);
    else policy_input_.resize(num_state, 1);

    q_dot_lpf_.setZero();

    torque_bound_ << 333, 232, 263, 289, 222, 166,
                    333, 232, 263, 289, 222, 166,
                    303, 303, 303, 
                    64, 64, 64, 64, 23, 23, 10, 10,
                    10, 10,
                    64, 64, 64, 64, 23, 23, 10, 10;  
                    
    q_init_ << 0.0, 0.0, -0.28, 0.6, -0.32, 0.0,
                0.0, 0.0, -0.28, 0.6, -0.32, 0.0,
                0.0, 0.0, 0.0,
                0.3, 0.3, 1.5, -1.27, -1.0, 0.0, -1.0, 0.0,
                0.0, 0.0,
                -0.3, -0.3, -1.5, 1.27, 1.0, 0.0, 1.0, 0.0;

    kp_.setZero();
    kv_.setZero();
    kp_.diagonal() <<   2000.0, 5000.0, 4000.0, 3700.0, 3200.0, 3200.0,
                        2000.0, 5000.0, 4000.0, 3700.0, 3200.0, 3200.0,
                        6000.0, 10000.0, 10000.0,
                        400.0, 1000.0, 400.0, 400.0, 400.0, 400.0, 100.0, 100.0,
                        100.0, 100.0,
                        400.0, 1000.0, 400.0, 400.0, 400.0, 400.0, 100.0, 100.0;
    kp_.diagonal() /= 9.0;  
    kv_.diagonal() << 15.0, 50.0, 20.0, 25.0, 24.0, 24.0,
                        15.0, 50.0, 20.0, 25.0, 24.0, 24.0,
                        200.0, 100.0, 100.0,
                        10.0, 28.0, 10.0, 10.0, 10.0, 10.0, 3.0, 3.0,
                        2.0, 2.0,
                        10.0, 28.0, 10.0, 10.0, 10.0, 10.0, 3.0, 3.0;
    kv_.diagonal() /= 3.0;

    // Woohyun
    initBias();
}

Eigen::Vector3d CustomController::mat2euler(Eigen::Matrix3d mat)
{
    Eigen::Vector3d euler;

    double cy = std::sqrt(mat(2, 2) * mat(2, 2) + mat(1, 2) * mat(1, 2));
    if (cy > std::numeric_limits<double>::epsilon())
    {
        euler(2) = -atan2(mat(0, 1), mat(0, 0));
        euler(1) =  -atan2(-mat(0, 2), cy);
        euler(0) = -atan2(mat(1, 2), mat(2, 2));
    }
    else
    {
        euler(2) = -atan2(-mat(1, 0), mat(1, 1));
        euler(1) =  -atan2(-mat(0, 2), cy);
        euler(0) = 0.0;
    }
    return euler;
}

void CustomController::initBias()
{
    q_bias_.setZero();
    if (~is_on_robot_){
        std::random_device rd;  
        std::mt19937 gen(rd());
        float bias_std = 0.;
        std::uniform_real_distribution<> dis(-bias_std, bias_std);
        q_bias_(2) = dis(gen);
        q_bias_(3) = dis(gen);
        q_bias_(4) = dis(gen);
        q_bias_(8) = dis(gen);
        q_bias_(9) = dis(gen);
        q_bias_(10) = dis(gen);
        // for (int i = 0; i < num_actuator_action; i++){
        //     q_bias_(i) = dis(gen);

        // }
    }
}

void CustomController::processBias()
{
    for (int i = 0; i < MODEL_DOF; i++){
        q_noise_(i) += q_bias_(i);
    }
}

void CustomController::processNoise()
{
    time_cur_ = rd_cc_.control_time_us_ / 1e6;
    if (is_on_robot_)
    {
        q_vel_noise_ = rd_cc_.q_dot_virtual_.segment(6,MODEL_DOF);
        q_noise_= rd_cc_.q_virtual_.segment(6,MODEL_DOF);
        if (time_cur_ - time_pre_ > 0.0)
        {
            q_dot_lpf_ = DyrosMath::lpf<MODEL_DOF>(q_vel_noise_, q_dot_lpf_, 1/(time_cur_ - time_pre_), 4.0);
        }
        else
        {
            q_dot_lpf_ = q_dot_lpf_;
        }
    }
    else
    {
        std::random_device rd;  
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(-0.00001, 0.00001);
        for (int i = 0; i < MODEL_DOF; i++) {
            q_noise_(i) = rd_cc_.q_virtual_(6+i) + dis(gen);
        }
        if (time_cur_ - time_pre_ > 0.0)
        {
            q_vel_noise_ = (q_noise_ - q_noise_pre_) / (time_cur_ - time_pre_);
            q_dot_lpf_ = DyrosMath::lpf<MODEL_DOF>(q_vel_noise_, q_dot_lpf_, 1/(time_cur_ - time_pre_), 4.0);
        }
        else
        {
            q_vel_noise_ = q_vel_noise_;
            q_dot_lpf_ = q_dot_lpf_;
        }
        q_noise_pre_ = q_noise_;
    }
    time_pre_ = time_cur_;
}

void CustomController::processObservation() // [linvel, angvel, proj_grav, commands, dof_pos, dof_vel, actions]
{

    int data_idx = 0;

    Eigen::Quaterniond q;
    q.x() = rd_cc_.q_virtual_(3);
    q.y() = rd_cc_.q_virtual_(4);
    q.z() = rd_cc_.q_virtual_(5);
    q.w() = rd_cc_.q_virtual_(MODEL_DOF_QVIRTUAL-1);   

    Vector3_t base_lin_vel, base_ang_vel;
    base_lin_vel = q.conjugate()*(rd_cc_.q_dot_virtual_.segment(0,3));
    base_ang_vel = q.conjugate()*(rd_cc_.q_dot_virtual_.segment(3,3));

    // for (int i=0; i<6; i++)
    // {
    //     state_cur_(data_idx) = rd_cc_.q_dot_virtual_(i);
    //     data_idx++;
    // }

    for (int i = 0; i < 3; i++){
        state_cur_(data_idx) = base_lin_vel(i);
        data_idx++;
    }

    for (int i = 0; i < 3; i++){
        state_cur_(data_idx) = base_ang_vel(i);
        data_idx++;
    }
 
    Vector3_t grav, projected_grav, forward_vec;
    grav << 0, 0, -1.;
    forward_vec << 1., 0, 0;
    projected_grav = q.conjugate()*grav;
    
    // euler_angle_ = DyrosMath::rot2Euler_tf(q.toRotationMatrix());

    state_cur_(data_idx) = projected_grav(0);
    data_idx++;

    state_cur_(data_idx) = projected_grav(1);
    data_idx++;

    state_cur_(data_idx) = projected_grav(2);
    data_idx++;

    // std::cout << "start time : " << start_time_ << std::endl;
    state_cur_(data_idx) = DyrosMath::cubic(rd_cc_.control_time_us_, start_time_, start_time_ + 3e6, pre_target_vel_x_, target_vel_x_, 0.0, 0.0);// .5;//target_vel_x_;
    // std::cout << "command : " << state_cur_(data_idx) << std::endl;
    // state_cur_(data_idx) = 0.5;//target_vel_x_;
    data_idx++;

    state_cur_(data_idx) = DyrosMath::cubic(rd_cc_.control_time_us_, start_time_, start_time_ + 3e6, pre_target_vel_y_, target_vel_y_, 0.0, 0.0); //target_vel_y_; //0.0;//target_vel_y_;
    data_idx++;

    Vector3_t forward = q * forward_vec;
    heading = atan2(forward(1), forward(0));
    // state_cur_(data_idx) = DyrosMath::minmax_cut(0.5*DyrosMath::wrap_to_pi(target_heading_ - heading), -1., 1.);
    state_cur_(data_idx) = DyrosMath::cubic(rd_cc_.control_time_us_, start_time_, start_time_ + 3e6, pre_target_vel_yaw_, DyrosMath::minmax_cut(0.5*DyrosMath::wrap_to_pi(target_heading_ - heading), -1., 1.), 0.0, 0.0);
    // ROS_INFO("Current heading : %f\n", heading);
    // ROS_INFO("Target heading : %f\n", target_heading_);
    // ROS_INFO("Target yaw vel : %f\n", state_cur_(data_idx));
    // state_cur_(data_idx) = 0.0;//target_vel_yaw;
    data_idx++;

    for (int i = 0; i < num_actuator_action; i++)
    {
        state_cur_(data_idx) = q_noise_(i) - q_init_(i);
        data_idx++;
    }

    for (int i = 0; i < num_actuator_action; i++)
    {
        if (is_on_robot_)
        {
            state_cur_(data_idx) = q_vel_noise_(i);
        }
        else
        {
            state_cur_(data_idx) = q_vel_noise_(i); //rd_cc_.q_dot_virtual_(i+6);
        }
        data_idx++;
    }

    // float squat_duration = 1.7995;
    // phase_ = std::fmod((rd_cc_.control_time_us_-start_time_)/1e6 + action_dt_accumulate_, squat_duration) / squat_duration;

    // state_cur_(data_idx) = sin(2*M_PI*phase_);
    // data_idx++;
    // state_cur_(data_idx) = cos(2*M_PI*phase_);
    // data_idx++;
    




    // state_cur_(data_idx) = -rd_cc_.LF_FT(2);
    // data_idx++;

    // state_cur_(data_idx) = -rd_cc_.RF_FT(2);
    // data_idx++;

    // state_cur_(data_idx) = rd_cc_.LF_FT(3);
    // data_idx++;

    // state_cur_(data_idx) = rd_cc_.RF_FT(3);
    // data_idx++;

    // state_cur_(data_idx) = rd_cc_.LF_FT(4);
    // data_idx++;

    // state_cur_(data_idx) = rd_cc_.RF_FT(4);
    // data_idx++;

    for (int i = 0; i <num_actuator_action; i++) 
    {
        state_cur_(data_idx) = DyrosMath::minmax_cut(rl_action_(i), -1.0, 1.0);
        data_idx++;
    }
    // state_cur_(data_idx) = DyrosMath::minmax_cut(rl_action_(num_actuator_action), 0.0, 1.0);
    // data_idx++;
    state_buffer_.block(0, 0, num_cur_state*(num_state_skip*num_state_hist-1),1) = state_buffer_.block(num_cur_state, 0, num_cur_state*(num_state_skip*num_state_hist-1),1);
    if (morphnet) state_history_.block(0, 0, num_cur_state, morph_history_skip_*morph_history_len_-1) = state_history_.block(0, 1, num_cur_state,morph_history_skip_*morph_history_len_-1);
    else state_history_.block(0, 0, num_cur_state, history_skip_*history_len_-1) = state_history_.block(0, 1, num_cur_state,history_skip_*history_len_-1);

    MatrixXd tmp = (state_cur_ - state_mean_).array() / state_var_.cwiseSqrt().array();
    for (int i = 0; i < num_cur_state; i++){
        tmp(i) = DyrosMath::minmax_cut(tmp(i), -10., 10.);
    }
    state_buffer_.block(num_cur_state*(num_state_skip*num_state_hist-1), 0, num_cur_state,1) = tmp;
    if (morphnet) state_history_.block(0,morph_history_skip_*morph_history_len_-1, num_cur_state, 1) = tmp;
    else state_history_.block(0,history_skip_*history_len_-1, num_cur_state, 1) = tmp;

    if (!encoder_initialized){
        if (morphnet){
            for (int i = 0; i < morph_history_skip_*morph_history_len_; i++) 
                state_history_.block(0,i, num_cur_state, 1) = tmp;
        }
        else{
            for (int i = 0; i < history_skip_*history_len_; i++) 
                state_history_.block(0,i, num_cur_state, 1) = tmp;
        }
        encoder_initialized = true;
    }

    // // Internal State First
    // for (int i = 0; i < num_state_hist; i++)
    // {
    //     state_.block(num_cur_internal_state*i, 0, num_cur_internal_state, 1) = state_buffer_.block(num_cur_state*(num_state_skip*(i+1)-1), 0, num_cur_internal_state, 1);
    // }
    // // Action History Second
    // for (int i = 0; i < num_state_hist-1; i++)
    // {
    //     state_.block(num_state_hist*num_cur_internal_state + num_action*i, 0, num_action, 1) = state_buffer_.block(num_cur_state*(num_state_skip*(i+1)) + num_cur_internal_state, 0, num_action, 1);
    // }
    for (int i = 0; i < num_state_hist; i++){
        state_.block(num_cur_state*i, 0, num_cur_state, 1) = state_buffer_.block(num_cur_state*(num_state_skip*(i+1)-1), 0, num_cur_state, 1);
    }
    if (morphnet){
        for (int i = 0; i < morph_history_len_; i++){
            morphnet_input_.block(num_cur_state * i, 0, num_cur_state, 1) = state_history_.block(0, morph_history_skip_*(i+1)-1, num_cur_state, 1);
        }
    }
    else if (use_encoder_){
        for (int i = 0; i < history_len_; i++){
            encoder_input_.block(0, i, num_cur_state, 1) = state_history_.block(0, history_skip_*(i+1)-1, num_cur_state, 1);
       }
    }



}

// RELU
// void CustomController::feedforwardPolicy()
// {
//     hidden_layer1_ = policy_net_w0_ * state_ + policy_net_b0_;
//     for (int i = 0; i < num_hidden; i++) 
//     {
//         if (hidden_layer1_(i) < 0)
//             hidden_layer1_(i) = 0.0;
//     }

//     hidden_layer2_ = policy_net_w2_ * hidden_layer1_ + policy_net_b2_;
//     for (int i = 0; i < num_hidden; i++) 
//     {
//         if (hidden_layer2_(i) < 0)
//             hidden_layer2_(i) = 0.0;
//     }

//     rl_action_ = action_net_w_ * hidden_layer2_ + action_net_b_;

//     value_hidden_layer1_ = value_net_w0_ * state_ + value_net_b0_;
//     for (int i = 0; i < num_hidden; i++) 
//     {
//         if (value_hidden_layer1_(i) < 0)
//             value_hidden_layer1_(i) = 0.0;
//     }

//     value_hidden_layer2_ = value_net_w2_ * value_hidden_layer1_ + value_net_b2_;
//     for (int i = 0; i < num_hidden; i++) 
//     {
//         if (value_hidden_layer2_(i) < 0)
//             value_hidden_layer2_(i) = 0.0;
//     }

//     value_ = (value_net_w_ * value_hidden_layer2_ + value_net_b_)(0);
    
// }

// ELU VERSION
void CustomController::feedforwardPolicy()
{
    if (morphnet){
        feedforwardMorphnet();
        policy_input_.block(0, 0, num_state, 1) = state_;
        policy_input_.block(num_state, 0, morph_params_dim, 1) = morphnet_output_; 
    }
    else if (use_encoder_){
        feedforwardEncoder();
        policy_input_.block(0, 0, num_state, 1) = state_;
        policy_input_.block(num_state, 0, encoder_dim_, 1) = encoder_output_;
        
    }
    else policy_input_ = state_;    
    // First hidden layer for policy network
    hidden_layer1_ = policy_net_w0_ * policy_input_ + policy_net_b0_;
    std::cout << "hidden layer 1" << std::endl;
    for (int i = 0; i < num_hidden1; i++) 
    {
        if (hidden_layer1_(i) < 0)
            hidden_layer1_(i) = std::exp(hidden_layer1_(i)) - 1.0;
    }

    // Second hidden layer for policy network
    hidden_layer2_ = policy_net_w2_ * hidden_layer1_ + policy_net_b2_;
    std::cout << "hidden layer 2" << std::endl;
    for (int i = 0; i < num_hidden2; i++) 
    {
        if (hidden_layer2_(i) < 0)
            hidden_layer2_(i) = std::exp(hidden_layer2_(i)) - 1.0;
    }

    // Output layer for policy network
    rl_action_ = action_net_w_ * hidden_layer2_ + action_net_b_;
    std::cout << "action" << std::endl;

    // First hidden layer for value network
    value_hidden_layer1_ = value_net_w0_ * policy_input_ + value_net_b0_;
    std::cout << "value layer 1" << std::endl;
    for (int i = 0; i < num_hidden1; i++) 
    {
        if (value_hidden_layer1_(i) < 0)
            value_hidden_layer1_(i) = std::exp(value_hidden_layer1_(i)) - 1.0;
    }

    // Second hidden layer for value network
    value_hidden_layer2_ = value_net_w2_ * value_hidden_layer1_ + value_net_b2_;
    std::cout << "value layer 2" << std::endl;

    for (int i = 0; i < num_hidden2; i++) 
    {
        if (value_hidden_layer2_(i) < 0)
            value_hidden_layer2_(i) = std::exp(value_hidden_layer2_(i)) - 1.0;
    }

    // Output layer for value network
    value_ = (value_net_w_ * value_hidden_layer2_ + value_net_b_)(0);
    std::cout << "value " << std::endl;

}

void CustomController::feedforwardEncoder()
{
    encoder_hidden_layer1_ = conv_layer(encoder_input_, encoder_conv1_w_, encoder_conv1_b_, encoder_conv1_kernel_, encoder_conv1_stride_);
    // encoder_hidden_layer1_ = encoder_hidden_layer1_.cwiseMax(0.);  // ReLU
    encoder_hidden_layer1_ = encoder_hidden_layer1_.array().tanh();  // Tanh

    encoder_hidden_layer2_ = conv_layer(encoder_hidden_layer1_, encoder_conv2_w_, encoder_conv2_b_, encoder_conv2_kernel_, encoder_conv2_stride_);
    // encoder_hidden_layer2_ = encoder_hidden_layer2_.cwiseMax(0.);  // ReLU
    encoder_hidden_layer2_ = encoder_hidden_layer2_.array().tanh();  // Tanh

    Eigen::MatrixXd flattened = Eigen::Map<const Eigen::MatrixXd>(encoder_hidden_layer2_.data(), encoder_hidden_layer2_.size(), 1);
    encoder_output_ = (encoder_fc_w_ * flattened) + encoder_fc_b_;
    encoder_output_ = encoder_output_.array().tanh();
}

void CustomController::feedforwardMorphnet(){
    // First hidden layer for policy network
    morph_hidden_layer1_ = morph_net_w0_ * morphnet_input_ + morph_net_b0_;
    morph_hidden_layer1_ = morph_hidden_layer1_.array().tanh();

    // Second hidden layer for policy network
    morph_hidden_layer2_ = morph_net_w2_ * morph_hidden_layer1_ + morph_net_b2_;
    morph_hidden_layer2_ = morph_hidden_layer2_.array().tanh();

    // Output layer for policy network
    morphnet_output_ = morph_net_w_ * morph_hidden_layer2_ + morph_net_b_;
    std::cout << "Morph params : " << morphnet_output_.transpose() << std::endl;
}

Eigen::MatrixXd CustomController::conv_layer(const MatrixXd &input, const MatrixXd &weights, const MatrixXd &biases, int kernel_size, int stride)
{
    // Implement 1D convolution manually
    // Here, input is of size (d, T), weights of size (out_channels, in_channels * kernel_size)
    int out_channels = biases.rows();
    int T_out = (input.cols() - kernel_size) / stride + 1;
    Eigen::MatrixXd output(out_channels, T_out);

    for (int i = 0; i < out_channels; ++i) {
        for (int t = 0; t < T_out; ++t) {
            int start_idx = t * stride;
            Eigen::VectorXd segment = Eigen::Map<const Eigen::VectorXd>(
                input.data() + start_idx * input.rows(), input.rows() * kernel_size);
            output(i, t) = segment.dot(weights.row(i)) + biases(i, 0);
        }
    }

    return output;
}

void CustomController::loadEncoderNetwork()
{
    string cur_path = "/home/cha/catkin_ws/src/tocabi_cc/";

    if (is_on_robot_)
    {
        cur_path = "/home/dyros/catkin_ws/src/tocabi_cc/";
    }

    std::ifstream file[6];
    file[0].open("/home/cha/isaac_ws/AMP_for_hardware/logs/encoder/conv1_weight.txt", std::ios::in);
    file[1].open("/home/cha/isaac_ws/AMP_for_hardware/logs/encoder/conv1_bias.txt", std::ios::in);
    file[2].open("/home/cha/isaac_ws/AMP_for_hardware/logs/encoder/conv2_weight.txt", std::ios::in);
    file[3].open("/home/cha/isaac_ws/AMP_for_hardware/logs/encoder/conv2_bias.txt", std::ios::in);
    file[4].open("/home/cha/isaac_ws/AMP_for_hardware/logs/encoder/fc2_weight.txt", std::ios::in);
    file[5].open("/home/cha/isaac_ws/AMP_for_hardware/logs/encoder/fc2_bias.txt", std::ios::in);

    if (!file[0].is_open()) {
        std::cout << "Cannot find the weight file" << std::endl;
    }

    float temp;
    int row = 0;
    int col = 0;

    // Load conv1 weights
    while (!file[0].eof() && row != encoder_conv1_w_.rows())
    {
        file[0] >> temp;
        encoder_conv1_w_(row, col) = temp;
        col++;
        if (col == encoder_conv1_w_.cols())
        {
            col = 0;
            row++;
        }
    }

    // Load conv1 biases
    row = 0;
    col = 0;
    while (!file[1].eof() && row != encoder_conv1_b_.rows())
    {
        file[1] >> temp;
        encoder_conv1_b_(row, col) = temp;
        col++;
        if (col == encoder_conv1_b_.cols())
        {
            col = 0;
            row++;
        }
    }

    // Load conv2 weights
    row = 0;
    col = 0;
    while (!file[2].eof() && row != encoder_conv2_w_.rows())
    {
        file[2] >> temp;
        encoder_conv2_w_(row, col) = temp;
        col++;
        if (col == encoder_conv2_w_.cols())
        {
            col = 0;
            row++;
        }
    }

    // Load conv2 biases
    row = 0;
    col = 0;
    while (!file[3].eof() && row != encoder_conv2_b_.rows())
    {
        file[3] >> temp;
        encoder_conv2_b_(row, col) = temp;
        col++;
        if (col == encoder_conv2_b_.cols())
        {
            col = 0;
            row++;
        }
    }

    // Load fully connected layer weights
    row = 0;
    col = 0;
    while (!file[4].eof() && row != encoder_fc_w_.rows())
    {
        file[4] >> temp;
        encoder_fc_w_(row, col) = temp;
        col++;
        if (col == encoder_fc_w_.cols())
        {
            col = 0;
            row++;
        }
    }

    // Load fully connected layer biases
    row = 0;
    col = 0;
    while (!file[5].eof() && row != encoder_fc_b_.rows())
    {
        file[5] >> temp;
        encoder_fc_b_(row, col) = temp;
        col++;
        if (col == encoder_fc_b_.cols())
        {
            col = 0;
            row++;
        }
    }
}

void CustomController::loadMorphnet()
{
    morphnet_input_.setZero();
    morphnet_output_.setZero();


    string cur_path = "/home/cha/catkin_ws/src/tocabi_cc/";

    if (is_on_robot_)
    {
        cur_path = "/home/dyros/catkin_ws/src/tocabi_cc/";
    }
    std::ifstream file[6];
    file[0].open("/home/cha/isaac_ws/AMP_for_hardware/logs/encoder/linear1_weight.txt", std::ios::in);
    file[1].open("/home/cha/isaac_ws/AMP_for_hardware/logs/encoder/linear1_bias.txt", std::ios::in);
    file[2].open("/home/cha/isaac_ws/AMP_for_hardware/logs/encoder/linear2_weight.txt", std::ios::in);
    file[3].open("/home/cha/isaac_ws/AMP_for_hardware/logs/encoder/linear2_bias.txt", std::ios::in);
    file[4].open("/home/cha/isaac_ws/AMP_for_hardware/logs/encoder/head_weight.txt", std::ios::in);
    file[5].open("/home/cha/isaac_ws/AMP_for_hardware/logs/encoder/head_bias.txt", std::ios::in);


    if(!file[0].is_open())
    {
        std::cout<<"Can not find the weight file"<<std::endl;
        exit(1);
    }

    float temp;
    int row = 0;
    int col = 0;

    while(!file[0].eof() && row != morph_net_w0_.rows())
    {
        file[0] >> temp;
        if(temp != '\n')
        {
            morph_net_w0_(row, col) = temp;
            col ++;
            if (col == morph_net_w0_.cols())
            {
                col = 0;
                row ++;
            }
        }
    }
    row = 0;
    col = 0;
    while(!file[1].eof() && row != morph_net_b0_.rows())
    {
        file[1] >> temp;
        if(temp != '\n')
        {
            morph_net_b0_(row, col) = temp;
            col ++;
            if (col == morph_net_b0_.cols())
            {
                col = 0;
                row ++;
            }
        }
    }
    row = 0;
    col = 0;
    while(!file[2].eof() && row != morph_net_w2_.rows())
    {
        file[2] >> temp;
        if(temp != '\n')
        {
            morph_net_w2_(row, col) = temp;
            col ++;
            if (col == morph_net_w2_.cols())
            {
                col = 0;
                row ++;
            }
        }
    }
    row = 0;
    col = 0;
    while(!file[3].eof() && row != morph_net_b2_.rows())
    {
        file[3] >> temp;
        if(temp != '\n')
        {
            morph_net_b2_(row, col) = temp;
            col ++;
            if (col == morph_net_b2_.cols())
            {
                col = 0;
                row ++;
            }
        }
    }
    row = 0;
    col = 0;
    while(!file[4].eof() && row != morph_net_w_.rows())
    {
        file[4] >> temp;
        if(temp != '\n')
        {
            morph_net_w_(row, col) = temp;
            col ++;
            if (col == morph_net_w_.cols())
            {
                col = 0;
                row ++;
            }
        }
    }
    row = 0;
    col = 0;
    while(!file[5].eof() && row != morph_net_b_.rows())
    {
        file[5] >> temp;
        if(temp != '\n')
        {
            morph_net_b_(row, col) = temp;
            col ++;
            if (col == morph_net_b_.cols())
            {
                col = 0;
                row ++;
            }
        }
    }
}

void CustomController::computeSlow()
{
    copyRobotData(rd_);
    if (rd_cc_.tc_.mode == 7)
    {
        if (rd_cc_.tc_init)
        {
            //Initialize settings for Task Control! 
            start_time_ = rd_cc_.control_time_us_;
            q_noise_pre_ = q_noise_ = q_init_ = rd_cc_.q_virtual_.segment(6,MODEL_DOF);
            time_cur_ = start_time_ / 1e6;
            time_pre_ = time_cur_ - 0.005;
            // time_inference_pre_ = rd_cc_.control_time_us_ - (1/249.9)*1e6;
            time_inference_pre_ = rd_cc_.control_time_us_ - (1/124.9)*1e6;

            rd_.tc_init = false;
            std::cout<<"cc mode 7"<<std::endl;
            torque_init_ = rd_cc_.torque_desired;

            processNoise();
            // Woohyun
            processBias();
            processObservation();
            for (int i = 0; i < num_state_skip*num_state_hist; i++) 
            {
                state_buffer_.block(num_cur_state*i, 0, num_cur_state, 1) = (state_cur_ - state_mean_).array() / state_var_.cwiseSqrt().array();
                // state_buffer_.block(num_cur_state*i, 0, num_cur_state, 1).setZero();
            }
        }

        processNoise();
        // Woohyun
        processBias();

        // processObservation and feedforwardPolicy mean time: 15 us, max 53 us
        // With encoder, 
        if ((rd_cc_.control_time_us_ - time_inference_pre_)/1.0e6 >= 1/125.0 - 4/10000.0) // 125 is the control frequency
        // if ((rd_cc_.control_time_us_ - time_inference_pre_)/1.0e6 >= 1/250.0 - 1/10000.0) // 250 is the control frequency
        {
            // auto start_time = std::chrono::high_resolution_clock::now();

            // Call the functions you want to measure
            processObservation();
            feedforwardPolicy();

            // End time measurement
            // auto end_time = std::chrono::high_resolution_clock::now();

            // // Calculate the duration in microseconds
            // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();

            // // Output the time taken
            // std::cout << "processObservation and feedforwardPolicy took " << duration << " us" << std::endl;

            
            // action_dt_accumulate_ += DyrosMath::minmax_cut(rl_action_(num_action-1)*5/250.0, 0.0, 5/250.0);
            action_dt_accumulate_ += DyrosMath::minmax_cut(rl_action_(num_action-1)*5/125.0, 0.0, 5/125.0);
            std::cout << "Value : " << value_ << std::endl;
            if (value_ < 60.0)
            {
                if (stop_by_value_thres_ == false)
                {
                    stop_by_value_thres_ = true;
                    stop_start_time_ = rd_cc_.control_time_us_;
                    q_stop_ = q_noise_;
                    std::cout << "Stop by Value Function" << std::endl;
                }
            }

            if (is_write_file_)
            {
                    writeFile << (rd_cc_.control_time_us_ - time_inference_pre_)/1e6 << "\t";
                    // writeFile << DyrosMath::minmax_cut(rl_action_(num_action-1)*1/100.0, 0.0, 1/100.0) << "\t";

                    // writeFile << rd_cc_.LF_FT.transpose() << "\t";
                    // writeFile << rd_cc_.RF_FT.transpose() << "\t";
                    writeFile << rd_cc_.LF_CF_FT.transpose() << "\t";
                    writeFile << rd_cc_.RF_CF_FT.transpose() << "\t";

                    writeFile << rd_cc_.torque_desired.transpose()  << "\t";
                    writeFile << q_noise_.transpose() << "\t";
                    writeFile << q_dot_lpf_.transpose() << "\t";
                    writeFile << rd_cc_.q_dot_virtual_.transpose() << "\t";
                    writeFile << rd_cc_.q_virtual_.transpose() << "\t";
                    writeFile << heading << "\t";

                    writeFile << value_ << "\t" << stop_by_value_thres_ << "\t";
                    writeFile << state_cur_(9) << "\t" << state_cur_(11) << "\t" << target_heading_ << "\t";
                    if (morphnet) writeFile << morphnet_output_.transpose() << "\t";
                    writeFile << std::endl;

                    time_write_pre_ = rd_cc_.control_time_us_;

                    // Data for actuator net training
                    if (actuator_net_log){
                        actuator_data_file << rd_cc_.torque_elmo_.segment(0,12).transpose() << "\t";
                        actuator_data_file << rd_cc_.q_dot_virtual_.segment(6, 12).transpose() << "\t";
                        actuator_data_file << rd_cc_.q_virtual_.segment(7,12).transpose() << "\t";
                        actuator_data_file << rd_cc_.torque_input_.segment(0,12).transpose() << "\t";
                        actuator_data_file << std::endl;
                    }

            }
            
            time_inference_pre_ = rd_cc_.control_time_us_;
        }

        for (int i = 0; i < num_actuator_action; i++)
        {
            // WH
            torque_rl_(i) = DyrosMath::minmax_cut(DyrosMath::minmax_cut(rl_action_(i), -1., 1.) *torque_bound_(i), -torque_bound_(i), torque_bound_(i));
        }
        for (int i = num_actuator_action; i < MODEL_DOF; i++)
        {
            torque_rl_(i) = kp_(i,i) * (q_init_(i) - q_noise_(i)) - kv_(i,i)*q_vel_noise_(i);
        }
        
        if (rd_cc_.control_time_us_ < start_time_ + 0.1e6)
        {

            for (int i = 0; i <MODEL_DOF; i++)
            {
                torque_spline_(i) = DyrosMath::cubic(rd_cc_.control_time_us_, start_time_, start_time_ + 0.1e6, torque_init_(i), torque_rl_(i), 0.0, 0.0);
            }
            rd_.torque_desired = torque_spline_;    
        }
        else
        {
             rd_.torque_desired = torque_rl_;
            
        }

        if (stop_by_value_thres_)
        {
            // rd_.torque_desired = kp_ * (q_stop_ - q_noise_) - kv_*q_vel_noise_;
        }

    }
}

void CustomController::computeFast()
{
    // if (tc.mode == 10)
    // {
    // }
    // else if (tc.mode == 11)
    // {
    // }
}

void CustomController::computePlanner()
{
}

void CustomController::copyRobotData(RobotData &rd_l)
{
    std::memcpy(&rd_cc_, &rd_l, sizeof(RobotData));
}

void CustomController::joyCallback(const sensor_msgs::Joy::ConstPtr& joy)
{
    target_vel_x_ = DyrosMath::minmax_cut(0.5*joy->axes[1], -0.2, 0.5);
    target_vel_y_ = DyrosMath::minmax_cut(0.5*joy->axes[0], -0.2, 0.2);
}

void CustomController::rlcommandCallback(const tocabi_msgs::RLCommand::ConstPtr& command){
    target_vel_x_ = DyrosMath::minmax_cut(command->forward, 0., 1.);
    target_vel_y_ = DyrosMath::minmax_cut(command->lateral * 0.5, -.5, .5);
    target_heading_ = DyrosMath::minmax_cut(command->heading * 3.14, -3.14, 3.14);
    pre_target_vel_x_ = state_cur_(9);
    pre_target_vel_y_ = state_cur_(10);
    pre_target_vel_yaw_ = state_cur_(11);


}