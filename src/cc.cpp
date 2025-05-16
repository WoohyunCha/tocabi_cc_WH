#include "cc.h"

using namespace TOCABI;

CustomController::CustomController(RobotData &rd) : rd_(rd) //, wbc_(dc.wbc_)
{
    ControlVal_.setZero();

    if (is_write_file_)
    {
        if (is_on_robot_)
        {
            writeFile.open("/home/dyros/catkin_ws/src/tocabi_cc/result/data.csv", std::ofstream::out);
        }
        else
        {
            writeFile.open("/home/cha/catkin_ws/src/tocabi_cc/result/data.csv", std::ofstream::out);
            evalFile.open("/home/cha/catkin_ws/src/tocabi_cc/result/eval_data.csv", std::ofstream::out);
        }
        writeFile << std::fixed << std::setprecision(8);
        evalFile << std::fixed << std::setprecision(8);
    }
    initVariable();
    std::cout << "Load network start\n" << std::endl;

    loadNetwork();
    std::cout << "Load network end\n" << std::endl;

    joy_sub_ = nh_.subscribe<sensor_msgs::Joy>("joy", 10, &CustomController::joyCallback, this);
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

    base_path = loadPathFromConfig(cur_path + "weight_directory.txt");


    std::ifstream file[14];
    file[0].open(base_path + "policy/0_weight.txt", std::ios::in);
    file[1].open(base_path + "policy/0_bias.txt", std::ios::in);
    file[2].open(base_path + "policy/2_weight.txt", std::ios::in);
    file[3].open(base_path + "policy/2_bias.txt", std::ios::in);
    file[4].open(base_path + "policy/4_weight.txt", std::ios::in);
    file[5].open(base_path + "policy/4_bias.txt", std::ios::in);
    file[6].open(base_path + "normalizer/running_mean.txt", std::ios::in);
    file[7].open(base_path + "normalizer/running_var.txt", std::ios::in);
    file[8].open(base_path + "critic/0_weight.txt", std::ios::in);
    file[9].open(base_path + "critic/0_bias.txt", std::ios::in);
    file[10].open(base_path + "critic/2_weight.txt", std::ios::in);
    file[11].open(base_path + "critic/2_bias.txt", std::ios::in);
    file[12].open(base_path + "critic/4_weight.txt", std::ios::in);
    file[13].open(base_path + "critic/4_bias.txt", std::ios::in);

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
    if (use_encoder_) loadEncoderNetwork();
}

void CustomController::initVariable()
{    
    // Load the path from the configuration file

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

    if (use_encoder_) policy_input_.resize(num_state + encoder_dim_, 1);
    else policy_input_.resize(num_state, 1);

    q_dot_lpf_.setZero();

    torque_bound_ << 333, 232, 263, 289, 222, 166,
                    333, 232, 263, 289, 222, 166,
                    303, 303, 303, 
                    64, 64, 64, 64, 23, 23, 10, 10,
                    10, 10,
                    64, 64, 64, 64, 23, 23, 10, 10;  
                    
    if (com_height_ == 0.68){
        q_init_ << 0.0, 0.0, -0.5, 1.06, -0.56, 0.0,
                    0.0, 0.0, -0.5, 1.06, -0.56, 0.0,
                    0.0, 0.0, 0.0,
                    0.3, 0.1, 1.5, -1.27, -0.3, 0.0, -0.0, 0.0,
                    0.0, 0.0,
                    -0.3, -0.1, -1.5, 1.27, 0.3, 0.0, 0.0, 0.0;
    }
    else if (com_height_ == 0.728){
        q_init_ << 0.0, 0.0, -0.24, 0.6, -0.36, 0.0,
                    0.0, 0.0, -0.24, 0.6, -0.36, 0.0,
                    0.0, 0.0, 0.0,
                    0.3, 0.3, 1.5, -1.27, -1.0, 0.0, -1.0, 0.0,
                    0.0, 0.0,
                    -0.3, -0.3, -1.5, 1.27, 1.0, 0.0, 1.0, 0.0;

    }
    else {
        std::cout << "WRONG COM HEIGHT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        exit(0);
    }

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
    base_lin_vel.setZero();
    base_ang_vel.setZero();
    swing_state_stance_frame_.setZero(13);
    com_state_stance_frame_.setZero(13);
    q_leg_desired_ = q_init_.segment(0, num_actuator_action);
    com_support_current_dot_prev_.setZero(3);

    string cur_path = "/home/cha/catkin_ws/src/tocabi_cc/";

    if (is_on_robot_)
    {
        cur_path = "/home/dyros/catkin_ws/src/tocabi_cc/";
    }


    if (ctrl_mode){
        foothold_x_planned.setZero(planned_step_number);
        foothold_y_planned.setZero(planned_step_number);
        foothold_yaw_planned.setZero(planned_step_number);
        t_dsp_planned.setZero(planned_step_number);
        t_ssp_planned.setZero(planned_step_number);
        foot_height_planned.setZero(planned_step_number);
        lfoot_global_state.setZero(3);
        rfoot_global_state.setZero(3);
        loadCommand(cur_path + "commands.txt");
    }
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
    
    base_lin_vel = q.conjugate()*(rd_cc_.q_dot_virtual_.segment(0,3));
    base_ang_vel = (rd_cc_.q_dot_virtual_.segment(3,3));

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

    Vector3_t forward = q * forward_vec;
    double heading = atan2(forward(1), forward(0));
    heading_error = target_heading - heading;

    state_cur_(data_idx) = q.x();
    data_idx++;
    state_cur_(data_idx) = q.y();
    data_idx++;
    state_cur_(data_idx) = q.z();
    data_idx++;
    state_cur_(data_idx) = q.w();
    data_idx++;

    for (int i = 0; i < num_actuator_action; i++)
    {
        state_cur_(data_idx) = q_noise_(i) - q_init_(i);
        data_idx++;
    }

    state_cur_(data_idx) = q_noise_(15) - q_init_(15);
    data_idx++;
    state_cur_(data_idx) = q_noise_(16) - q_init_(16);
    data_idx++;
    // state_cur_(data_idx) = q_noise_(17) - q_init_(17);
    // data_idx++;
    state_cur_(data_idx) = q_noise_(19) - q_init_(19);
    data_idx++;


    state_cur_(data_idx) = q_noise_(25) - q_init_(25);
    data_idx++;
    state_cur_(data_idx) = q_noise_(26) - q_init_(26);
    data_idx++;
    // state_cur_(data_idx) = q_noise_(27) - q_init_(27);
    // data_idx++;
    state_cur_(data_idx) = q_noise_(29) - q_init_(29);
    data_idx++;
    

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

    state_cur_(data_idx) = q_vel_noise_(15); //rd_cc_.q_dot_virtual_(i+6);
    data_idx++;
    state_cur_(data_idx) = q_vel_noise_(16); //rd_cc_.q_dot_virtual_(i+6);
    data_idx++;
    // state_cur_(data_idx) = q_vel_noise_(17); //rd_cc_.q_dot_virtual_(i+6);
    // data_idx++;
    state_cur_(data_idx) = q_vel_noise_(19); //rd_cc_.q_dot_virtual_(i+6);
    data_idx++;

    state_cur_(data_idx) = q_vel_noise_(25); //rd_cc_.q_dot_virtual_(i+6);
    data_idx++;
    state_cur_(data_idx) = q_vel_noise_(26); //rd_cc_.q_dot_virtual_(i+6);
    data_idx++;
    // state_cur_(data_idx) = q_vel_noise_(27); //rd_cc_.q_dot_virtual_(i+6);
    // data_idx++;
    state_cur_(data_idx) = q_vel_noise_(29); //rd_cc_.q_dot_virtual_(i+6);
    data_idx++;


    for (int i = 0; i < num_actuator_action; i++)
    {
        // state_cur_(data_idx) = q_leg_desired_(i) - q_noise_(i);
        state_cur_(data_idx) = q_leg_desired_(i);
        data_idx++;
    }

    // state_cur_(data_idx) = com_state_stance_frame_(0);
    // data_idx++;
    // state_cur_(data_idx) = com_state_stance_frame_(1);
    // data_idx++;
    // state_cur_(data_idx) = com_state_stance_frame_(2);
    // data_idx++;
    // state_cur_(data_idx) = com_state_stance_frame_(3);
    // data_idx++;
    // state_cur_(data_idx) = com_state_stance_frame_(4);
    // data_idx++;
    // state_cur_(data_idx) = com_state_stance_frame_(5);
    // data_idx++;
    // state_cur_(data_idx) = com_state_stance_frame_(6);
    // data_idx++;
    // state_cur_(data_idx) = swing_state_stance_frame_(0);
    // data_idx++;
    // state_cur_(data_idx) = swing_state_stance_frame_(1);
    // data_idx++;
    // state_cur_(data_idx) = swing_state_stance_frame_(2);
    // data_idx++;
    // state_cur_(data_idx) = swing_state_stance_frame_(3);
    // data_idx++;
    // state_cur_(data_idx) = swing_state_stance_frame_(4);
    // data_idx++;
    // state_cur_(data_idx) = swing_state_stance_frame_(5);
    // data_idx++;
    // state_cur_(data_idx) = swing_state_stance_frame_(6);
    // data_idx++;



    state_cur_(data_idx) = cos(float(walking_tick) / float(t_total_(0)) * 2 * M_PI);
    data_idx++;
    state_cur_(data_idx) = sin(float(walking_tick) / float(t_total_(0)) * 2 * M_PI);
    data_idx++;

    state_cur_(data_idx) = step_length_x_(0);
    data_idx++;
    state_cur_(data_idx) = step_length_y_(0);
    data_idx++;
    state_cur_(data_idx) = step_yaw_(0);
    data_idx++;
    state_cur_(data_idx) = t_dsp_seconds(0);
    data_idx++;
    state_cur_(data_idx) = t_ssp_seconds(0);
    data_idx++;
    state_cur_(data_idx) = foot_height_(0);
    data_idx++;

    for (int i = 0; i <num_action; i++) 
    {
        state_cur_(data_idx) = DyrosMath::minmax_cut(rl_action_(i), -1.0, 1.0);
        data_idx++;
    }

    state_buffer_.block(0, 0, num_cur_state*(num_state_skip*num_state_hist-1),1) = state_buffer_.block(num_cur_state, 0, num_cur_state*(num_state_skip*num_state_hist-1),1);
    state_history_.block(0, 0, num_cur_state, history_skip_*history_len_-1) = state_history_.block(0, 1, num_cur_state,history_skip_*history_len_-1);

    MatrixXd tmp = (state_cur_ - state_mean_).array() / state_var_.cwiseSqrt().array();
    for (int i = 0; i < num_cur_state; i++){
        tmp(i) = DyrosMath::minmax_cut(tmp(i), -10., 10.);
    }
    state_buffer_.block(num_cur_state*(num_state_skip*num_state_hist-1), 0, num_cur_state,1) = tmp;
    state_history_.block(0,history_skip_*history_len_-1, num_cur_state, 1) = tmp;

    if (!encoder_initialized){
        for (int i = 0; i < history_skip_*history_len_; i++) 
            state_history_.block(0,i, num_cur_state, 1) = tmp;
        
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
    if (use_encoder_){
        for (int i = 0; i < history_len_; i++){
            encoder_input_.block(0, i, num_cur_state, 1) = state_history_.block(0, history_skip_*(i+1)-1, num_cur_state, 1);
       }
    }



}



// // RELU
// void CustomController::feedforwardPolicy()
// {
//     hidden_layer1_ = policy_net_w0_ * state_ + policy_net_b0_;
//     for (int i = 0; i < num_hidden1; i++) 
//     {
//         if (hidden_layer1_(i) < 0)
//             hidden_layer1_(i) = 0.0;
//     }

//     hidden_layer2_ = policy_net_w2_ * hidden_layer1_ + policy_net_b2_;
//     for (int i = 0; i < num_hidden2; i++) 
//     {
//         if (hidden_layer2_(i) < 0)
//             hidden_layer2_(i) = 0.0;
//     }

//     rl_action_ = action_net_w_ * hidden_layer2_ + action_net_b_;

//     value_hidden_layer1_ = value_net_w0_ * state_ + value_net_b0_;
//     for (int i = 0; i < num_hidden1; i++) 
//     {
//         if (value_hidden_layer1_(i) < 0)
//             value_hidden_layer1_(i) = 0.0;
//     }

//     value_hidden_layer2_ = value_net_w2_ * value_hidden_layer1_ + value_net_b2_;
//     for (int i = 0; i < num_hidden2; i++) 
//     {
//         if (value_hidden_layer2_(i) < 0)
//             value_hidden_layer2_(i) = 0.0;
//     }

//     value_ = (value_net_w_ * value_hidden_layer2_ + value_net_b_)(0);
    
// }

// ELU VERSION
void CustomController::feedforwardPolicy()
{
    if (use_encoder_){
        feedforwardEncoder();
        policy_input_.block(0, 0, num_state, 1) = state_;
        policy_input_.block(num_state, 0, encoder_dim_, 1) = encoder_output_;
        
    }
    else policy_input_ = state_;    
    // First hidden layer for policy network
    hidden_layer1_ = policy_net_w0_ * policy_input_ + policy_net_b0_;
    for (int i = 0; i < num_hidden1; i++) 
    {
        if (hidden_layer1_(i) < 0)
            hidden_layer1_(i) = std::exp(hidden_layer1_(i)) - 1.0;
    }

    // Second hidden layer for policy network
    hidden_layer2_ = policy_net_w2_ * hidden_layer1_ + policy_net_b2_;
    for (int i = 0; i < num_hidden2; i++) 
    {
        if (hidden_layer2_(i) < 0)
            hidden_layer2_(i) = std::exp(hidden_layer2_(i)) - 1.0;
    }

    // Output layer for policy network
    rl_action_ = action_net_w_ * hidden_layer2_ + action_net_b_;

    // First hidden layer for value network
    value_hidden_layer1_ = value_net_w0_ * policy_input_ + value_net_b0_;
    for (int i = 0; i < num_hidden1; i++) 
    {
        if (value_hidden_layer1_(i) < 0)
            value_hidden_layer1_(i) = std::exp(value_hidden_layer1_(i)) - 1.0;
    }

    // Second hidden layer for value network
    value_hidden_layer2_ = value_net_w2_ * value_hidden_layer1_ + value_net_b2_;

    for (int i = 0; i < num_hidden2; i++) 
    {
        if (value_hidden_layer2_(i) < 0)
            value_hidden_layer2_(i) = std::exp(value_hidden_layer2_(i)) - 1.0;
    }

    // Output layer for value network
    value_ = (value_net_w_ * value_hidden_layer2_ + value_net_b_)(0);

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
    base_path = loadPathFromConfig(cur_path + "weight_directory.txt");

    std::ifstream file[6];
    file[0].open(base_path + "encoder/conv1_weight.txt", std::ios::in);
    file[1].open(base_path + "encoder/conv1_bias.txt", std::ios::in);
    file[2].open(base_path + "encoder/conv2_weight.txt", std::ios::in);
    file[3].open(base_path + "encoder/conv2_bias.txt", std::ios::in);
    file[4].open(base_path + "encoder/fc2_weight.txt", std::ios::in);
    file[5].open(base_path + "encoder/fc2_bias.txt", std::ios::in);

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



void CustomController::joyCallback(const sensor_msgs::Joy::ConstPtr& joy)
{   

    last_buttons_0.push_back(joy->buttons[0]);
    last_buttons_1.push_back(joy->buttons[1]);
    last_buttons_2.push_back(joy->buttons[2]);
    last_buttons_3.push_back(joy->buttons[3]);
    last_buttons_7.push_back(joy->axes[7]);
    // joy_length_previous_vec.push_back(joy_length);

    if (last_buttons_0.size() > 10) {
        last_buttons_0.erase(last_buttons_0.begin());
    }
    if (last_buttons_1.size() > 10) {
        last_buttons_1.erase(last_buttons_1.begin()); 
    }
    if (last_buttons_2.size() > 10) {
        last_buttons_2.erase(last_buttons_2.begin());
    }
    if (last_buttons_3.size() > 10) {
        last_buttons_3.erase(last_buttons_3.begin()); 
    }
    if (last_buttons_7.size() > 10) {
        last_buttons_7.erase(last_buttons_7.begin());
    }
    if(iter_x_l>1){
        joy_length_previous = joy_length_command;
        iter_x_l =0;
    }
    if(iter_x_r>1){
        joy_length_previous = joy_length_command;
        iter_x_r =0;
    }
    if(iter_end_x>3){
        joy_length_previous = joy_length_command;
        iter_end_x =0;
    }



    if(iter_y_l_ls>1){
        joy_length_y_l_previous = joy_length_y_l_command;
        iter_y_l_ls =0;
    }
    if(iter_y_l_rs>2){
        joy_length_y_l_previous = joy_length_y_l_command;
        iter_y_l_rs =0;
    }
    if(iter_end_y_l_ls>1){
        joy_length_y_l_previous = joy_length_y_l_command;
        iter_end_y_l_ls =0;
    }
    if(iter_end_y_l_rs>2){
        joy_length_y_l_previous = joy_length_y_l_command;
        iter_end_y_l_rs =0;
    }



    if(iter_y_r_ls>2){
        joy_length_y_r_previous = joy_length_y_r_command;
        iter_y_r_ls =0;
    }
    if(iter_y_r_rs>1){
        joy_length_y_r_previous = joy_length_y_r_command;
        iter_y_r_rs =0;
    }
    if(iter_end_y_r_ls>2){
        joy_length_y_r_previous = joy_length_y_r_command;
        iter_end_y_r_ls =0;
    }
    if(iter_end_y_r_rs>1){
        joy_length_y_r_previous = joy_length_y_r_command;
        iter_end_y_r_rs =0;
    }


    
    if(fabs(joy->axes[0]) < 0.1 && fabs(joy->axes[1]) < 0.1){ 
        
        joy_length =0.0;
        joy_length_y_r_temp =0.21;
        joy_length_y_l_temp =0.21;
        // if(joy_length = 0.0){
        //     joy_length_previous =0.0;
        // }
        // Lcommand_step_length_x_ = joy_length_l;
        Lcommand_step_length_x_ = 0.;
        // std::cout<<joy_length_r<<std::endl;
        // Lcommand_step_length_y_ = joy_length_y_l;
        Lcommand_step_length_y_ = 0.21;
        // Lcommand_step_yaw_ = joy_yaw_l;
        // std::cout << "joy yaw l :" << joy_yaw_l << std::endl;
        // Lcommand_step_yaw_ = 0.;
       // Lcommand_t_dsp_ = 0.2;
        // Lcommand_t_ssp_ = 1.0;
        // Lcommand_foot_height_ = 0.1; 
        // Rcommand_step_length_x_ = joy_length_r;
        Rcommand_step_length_x_ = 0.;
        // Rcommand_step_length_y_ = joy_length_y_r;
        Rcommand_step_length_y_ = 0.21;
        // Rcommand_step_yaw_ = joy_yaw_r;
        // std::cout << "joy yaw r :" << joy_yaw_r << std::endl;
        // Rcommand_step_yaw_ = 0.;
        // Rcommand_t_dsp_ = 0.2;
        // Rcommand_t_ssp_ = 1.0;
        // Rcommand_foot_height_ = 0.1;

    }else{

        joy_x = joy->axes[1];
        joy_y = joy->axes[0];
        norm = sqrt(pow(joy_x,2) + pow(joy_y,2));

        joy_x /=norm;
        joy_y /=norm;

        
        if(joy_x > 0.0){
            joy_length = joy_length_temp*joy_x;
            // joy_length_l =joy_length_temp*joy_x;
            // joy_length_r =joy_length_temp*joy_x;
            // std::cout<<joy_length<<std::endl;
            Lcommand_step_length_x_ = joy_length;
            Rcommand_step_length_x_ = joy_length;
        }else if(joy_x < 0.0){
            joy_length = DyrosMath::minmax_cut(joy_length_temp*joy_x, -0.3, 0.0);
            // joy_length_l =joy_length_temp*-joy_x;
            // joy_length_r =joy_length_temp*-joy_x;
                // std::cout<<joy_length<<std::endl;
            Lcommand_step_length_x_ = joy_length;
            Rcommand_step_length_x_ = joy_length;
        }

        if(joy_y > 0.0){

            joy_length_y_l_temp = 0.21 + (joy_length_y_temp - 0.21)*joy_y; 
            joy_length_y_r_temp =0.21;
            Lcommand_step_length_y_= joy_length_y_l_temp;
            Rcommand_step_length_y_= joy_length_y_r_temp;
            // std::cout<<joy_length_y_r<<std::endl;
            // std::cout<<joy_y<<std::endl;
            // std::cout<<Lcommand_step_length_y_<<std::endl;
        }else if(joy_y < 0.0){

            joy_length_y_r_temp = 0.21 +  (joy_length_y_temp - 0.21)*-joy_y;
            joy_length_y_l_temp = 0.21;
            Lcommand_step_length_y_= joy_length_y_l_temp;             
            Rcommand_step_length_y_= joy_length_y_r_temp;
            // std::cout<<Rcommand_step_length_y_<<std::endl;
        }
    }
    //  std::cout<<joy_length_y_r<<std::endl;
    

    if(joy->axes[7] == -1.0 && Lcommand_t_dsp_ < 0.2 && Lcommand_t_ssp_<1.0){ //long step
        if(last_buttons_7[last_buttons_7.size()-2]==0.0){
            Lcommand_t_dsp_ = DyrosMath::minmax_cut(Lcommand_t_dsp_+0.01, 0.02, 0.2);
            Lcommand_t_ssp_ = DyrosMath::minmax_cut(Lcommand_t_ssp_+0.05, 0.5, 1.0);
            Rcommand_t_dsp_ = Lcommand_t_dsp_;
            Rcommand_t_ssp_ = Lcommand_t_ssp_;
            ROS_INFO("%f",Lcommand_t_dsp_);
            ROS_INFO("%f",Lcommand_t_ssp_);
        }

    }else if(joy->axes[7] == 1.0 && Lcommand_t_dsp_> 0.02 && Lcommand_t_ssp_> 0.5){//short step
        if(last_buttons_7[last_buttons_7.size()-2]==0.0){
            Lcommand_t_dsp_ = DyrosMath::minmax_cut(Lcommand_t_dsp_-0.01, 0.02, 0.2);
            Lcommand_t_ssp_ = DyrosMath::minmax_cut(Lcommand_t_ssp_-0.05, 0.5, 1.0);
            Rcommand_t_dsp_ = Lcommand_t_dsp_;
            Rcommand_t_ssp_ = Lcommand_t_ssp_;
            ROS_INFO("%f",Lcommand_t_dsp_);
            ROS_INFO("%f",Lcommand_t_ssp_);
            }
    }else{
        Lcommand_t_dsp_ = Lcommand_t_dsp_;
        Lcommand_t_ssp_ = Lcommand_t_ssp_;

        Rcommand_t_dsp_ = Lcommand_t_dsp_;
        Rcommand_t_ssp_ = Lcommand_t_ssp_;
    }
    
    if(joy->buttons[1] == 1 && joy_length_temp < 0.6 && joy_length_y_temp < 0.5){ 
        if(last_buttons_1[last_buttons_1.size()-2]==0){

            joy_length_temp+=0.08;//보폭
            joy_length_y_temp +=0.02;
            joy_length_temp = DyrosMath::minmax_cut(joy_length_temp, 0.1, 0.6);
            joy_length_y_temp = DyrosMath::minmax_cut(joy_length_y_temp, 0.25, 0.5);

            ROS_INFO("step_length_x_ : %f",joy_length_temp);
            ROS_INFO("step_length_y_ : %f",joy_length_y_temp);
        }            
    }
    if(joy->buttons[0] == 1 &&  joy_length_temp > 0.1 && joy_length_y_temp > 0.25){ 
        if(last_buttons_0[last_buttons_0.size()-2]==0){
            // 0.06 0.12 0.18 0.24 0.30 0.36
            // 0.25 0.29 0.33 0.37 0.41 0.45         

            joy_length_temp-=0.08;//보폭
            joy_length_y_temp-=0.02;//보폭
            joy_length_temp = DyrosMath::minmax_cut(joy_length_temp, 0.1, 0.6);
            joy_length_y_temp = DyrosMath::minmax_cut(joy_length_y_temp, 0.25, 0.5);
            ROS_INFO("step_length_x_ : %f",joy_length_temp);
            ROS_INFO("step_length_y_ : %f",joy_length_y_temp);
            
        }            
    }

    if(joy->buttons[3]==1 && Lcommand_foot_height_<0.2){
        if(last_buttons_3[last_buttons_3.size()-2]==0){
            joy_height= DyrosMath::minmax_cut(Lcommand_foot_height_+0.02, 0.05, 0.2);
            Lcommand_foot_height_= Rcommand_foot_height_ = joy_height;
            ROS_INFO("%f",joy_height);
        }

    }else if(joy->buttons[2]==1 && Lcommand_foot_height_>0.05){
        if(last_buttons_2[last_buttons_2.size()-2]==0){
            joy_height= DyrosMath::minmax_cut(Lcommand_foot_height_-0.02, 0.05, 0.2);
        Lcommand_foot_height_= Rcommand_foot_height_ = joy_height;
        ROS_INFO("%f",joy_height);
        }
    }else{
        Lcommand_foot_height_= joy_height;
        Rcommand_foot_height_= joy_height;
    }

    if(joy->buttons[4] == 1){
        joy_yaw_l_command =0.3;
        joy_yaw_r_command =0.0;

        Lcommand_step_yaw_ = joy_yaw_l_command;           
        Rcommand_step_yaw_ = joy_yaw_r_command;
    }
    if(joy->buttons[5] == 1){
        joy_yaw_l_command =0.0;
        joy_yaw_r_command =0.3;

        Lcommand_step_yaw_ = joy_yaw_l_command;           
        Rcommand_step_yaw_ = joy_yaw_r_command;

    }
    if(joy->buttons[5] != 1 && joy->buttons[4] != 1){
        joy_yaw_l_command =0.0;
        joy_yaw_r_command =0.0;
        
        Lcommand_step_yaw_ = joy_yaw_l_command;           
        Rcommand_step_yaw_ = joy_yaw_r_command;
    }

    
    // target_vel_x_ = DyrosMath::minmax_cut(0.5*joy->axes[1], -0.2, 0.5);
    // target_vel_y_ = DyrosMath::minmax_cut(0.5*joy->axes[0], -0.2, 0.2);
    // std::cout << "Rcommand_step_yaw_ :" << Rcommand_step_yaw_ << std::endl;
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

            q_leg_desired_ = rd_cc_.q_.segment(0,12);

            time_cur_ = start_time_ / 1e6;

            time_pre_ = time_cur_ - 0.005;

            // time_inference_pre_ = rd_cc_.control_time_us_ - (1/249.9)*1e6;

            time_inference_pre_ = rd_cc_.control_time_us_ - (1/(hz_))*1e6;



            rd_.tc_init = false;

            std::cout<<"cc mode 7"<<std::endl;

            torque_init_ = rd_cc_.torque_desired;

            target_com_state_stance_frame_.setZero(13);

            target_swing_state_stance_frame_.setZero(13);

            updateFootstepCommand();

            getRobotState();

            walkingStateMachine();

            getComTrajectory(); 

            getFootTrajectory();



            getTargetState();

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

        if ((rd_cc_.control_time_us_ - time_inference_pre_)/1.0e6 >= 1/hz_) // 125 is the control frequency

        {

            auto start_time = std::chrono::high_resolution_clock::now();



            // Call the functions you want to measure

            updateFootstepCommand();

            getRobotState();

            walkingStateMachine();

            getComTrajectory(); 

            getFootTrajectory();



            getTargetState();



            processObservation();

            feedforwardPolicy();


            
            updateNextStepTime();



            // End time measurement

            auto end_time = std::chrono::high_resolution_clock::now();



            // // Calculate the duration in microseconds

            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();



            // // Output the time taken

            // std::cout << "processObservation and feedforwardPolicy took " << duration << " us" << std::endl;



            

            // action_dt_accumulate_ += DyrosMath::minmax_cut(rl_action_(num_action-1)*5/250.0, 0.0, 5/250.0);

            action_dt_accumulate_ += DyrosMath::minmax_cut(rl_action_(num_action-1)*5/hz_, 0.0, 5/hz_);


            if (value_ < 60.0)

            {

                if (stop_by_value_thres_ == false)

                {



                    stop_by_value_thres_ = true;

                    stop_start_time_ = rd_cc_.control_time_us_;

                    q_stop_ = q_noise_;

                    std::cout << "Stop by Value Function : " << walking_tick << ", Value : " << value_ << std::endl;

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

                    writeFile << base_lin_vel.transpose() << "\t" << base_ang_vel.transpose() << "\t" << rd_cc_.q_dot_virtual_.segment(6,33).transpose() << "\t";

                    writeFile << rd_cc_.q_virtual_.transpose() << "\t";

                    writeFile << heading << "\t";



                    writeFile << value_ << "\t" << stop_by_value_thres_ << "\t";

                    writeFile << target_swing_state_stance_frame_.transpose() << "\t";

                    writeFile << target_com_state_stance_frame_.transpose() << "\t";

                    writeFile << swing_state_stance_frame_.transpose() << "\t";

                    writeFile << com_state_stance_frame_.transpose() << "\t";

                    writeFile << q_leg_desired_.transpose() << "\t";

                    writeFile << ref_zmp_(walking_tick,0) << "\t";

                    writeFile << ref_zmp_(walking_tick, 1) << "\t";

                    
                    writeFile << std::endl;



                    time_write_pre_ = rd_cc_.control_time_us_;




            }

            

            time_inference_pre_ = rd_cc_.control_time_us_;

        }







        for (int i = 0; i < num_actuator_action; i++)

        {

            // WH

            torque_rl_(i) = DyrosMath::minmax_cut(rl_action_(i), -1., 1.) *torque_bound_(i) ;

            // torque_rl_(i) = DyrosMath::minmax_cut(torque_rl_(i)  + 

            // DyrosMath::minmax_cut(kp_(i,i) * (q_leg_desired_(i) - q_noise_(i)) - kv_(i,i)*q_vel_noise_(i),-torque_bound_(i), torque_bound_(i) ),

            // -torque_bound_(i), torque_bound_(i));

        }

        for (int i = num_actuator_action; i < MODEL_DOF; i++)

        {

            torque_rl_(i) = kp_(i,i) * (q_init_(i) - q_noise_(i)) - kv_(i,i)*q_vel_noise_(i);

        }

        // torque_rl_(15) = DyrosMath::minmax_cut(rl_action_(12), -1., 1.) *torque_bound_(15);
        // torque_rl_(16) = DyrosMath::minmax_cut(rl_action_(13), -1., 1.) *torque_bound_(16);
        // // torque_rl_(17) = DyrosMath::minmax_cut(rl_action_(14), -1., 1.) *torque_bound_(17);
        // torque_rl_(19) = DyrosMath::minmax_cut(rl_action_(15), -1., 1.) *torque_bound_(19);

        // torque_rl_(25) = DyrosMath::minmax_cut(rl_action_(16), -1., 1.) *torque_bound_(25);
        // torque_rl_(26) = DyrosMath::minmax_cut(rl_action_(17), -1., 1.) *torque_bound_(26);
        // // torque_rl_(27) = DyrosMath::minmax_cut(rl_action_(18), -1., 1.) *torque_bound_(27);
        // torque_rl_(29) = DyrosMath::minmax_cut(rl_action_(19), -1., 1.) *torque_bound_(29);

        torque_rl_(15) = DyrosMath::minmax_cut(rl_action_(12), -1., 1.) *torque_bound_(15);
        torque_rl_(16) = DyrosMath::minmax_cut(rl_action_(13), -1., 1.) *torque_bound_(16);
        torque_rl_(19) = DyrosMath::minmax_cut(rl_action_(14), -1., 1.) *torque_bound_(19);

        torque_rl_(25) = DyrosMath::minmax_cut(rl_action_(15), -1., 1.) *torque_bound_(25);
        torque_rl_(26) = DyrosMath::minmax_cut(rl_action_(16), -1., 1.) *torque_bound_(26);
        torque_rl_(29) = DyrosMath::minmax_cut(rl_action_(17), -1., 1.) *torque_bound_(29);

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

            rd_.torque_desired = kp_ * (q_stop_ - q_noise_) - kv_*q_vel_noise_;

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


std::string CustomController::loadPathFromConfig(const std::string &config_file)
{
    std::cout << "LOAD WEIGHT!!" << std::endl;
    std::ifstream file(config_file);
    if (!file.is_open())
    {
        throw std::runtime_error("Cannot open configuration file: " + config_file);
    }

    std::string line, key, value;
    while (std::getline(file, line))
    {
        std::istringstream line_stream(line);
        if (std::getline(line_stream, key, '=') && std::getline(line_stream, value))
        {
            if (key == "weights_path")
            {
                file.close();
                return value; // Return the weights path
            }
        }
    }

    file.close();
    throw std::runtime_error("weights_path not found in configuration file.");
}

void CustomController::loadCommand(const std::string &command_file)
    {
        std::ifstream file(command_file);
        if (!file)
        {
        throw std::runtime_error("Cannot open command file: " + command_file);
        }

        std::string line;
        while (std::getline(file, line))
        {
        if (line.empty())
        continue;

        std::istringstream iss(line);
        std::string key;
        iss >> key;

        VectorXd vec;
        vec.setZero(planned_step_number);
        for (int i = 0; i < planned_step_number; i++)
        {
        if(!(iss >> vec(i)))
        {
        throw std::runtime_error("Error parsing values for key: " + key);
        }
        }

        if (key == "foothold_x_planned")
        foothold_x_planned = vec;
        else if (key == "foothold_y_planned")
        foothold_y_planned = vec;
        else if (key == "foothold_yaw_planned")
        foothold_yaw_planned = vec;
        else if (key == "t_dsp_planned")
        t_dsp_planned = vec;
        else if (key == "t_ssp_planned")
        t_ssp_planned = vec;
        else if (key == "foot_height_planned")
        foot_height_planned = vec;
        else
        std::cerr << "Warning: Unknown key '" << key << "' in file " << command_file << std::endl;
        }

    file.close();
    }

void CustomController::updateInitialState()
{

    pelv_rpy_current_.setZero();
    pelv_rpy_current_ = DyrosMath::rot2Euler(rd_cc_.link_[Pelvis].rotm); //ZYX multiply

    pelv_yaw_rot_current_from_global_ = DyrosMath::rotateWithZ(pelv_rpy_current_(2));
    pelv_yaw_rot_current_from_global_.linear() = rd_cc_.link_[Pelvis].rotm;

    rfoot_rpy_current_.setZero();
    lfoot_rpy_current_.setZero();
    rfoot_rpy_current_ = DyrosMath::rot2Euler(rd_cc_.link_[Right_Foot].rotm);
    lfoot_rpy_current_ = DyrosMath::rot2Euler(rd_cc_.link_[Left_Foot].rotm);

    rfoot_roll_rot_ = DyrosMath::rotateWithX(rfoot_rpy_current_(0));
    lfoot_roll_rot_ = DyrosMath::rotateWithX(lfoot_rpy_current_(0));
    rfoot_pitch_rot_ = DyrosMath::rotateWithY(rfoot_rpy_current_(1));
    lfoot_pitch_rot_ = DyrosMath::rotateWithY(lfoot_rpy_current_(1));
    rfoot_yaw_rot_ = DyrosMath::rotateWithZ(rfoot_rpy_current_(2));
    lfoot_yaw_rot_ = DyrosMath::rotateWithZ(lfoot_rpy_current_(2));


    if (foot_step_(0, 6) == 0) //right foot support
    {
        supportfoot_global_init_.translation() = rd_cc_.link_[Right_Foot].xpos;
        supportfoot_global_init_.linear() = rd_cc_.link_[Right_Foot].rotm;
    }
    else if (foot_step_(0, 6) == 1)
    {
        supportfoot_global_init_.translation() = rd_cc_.link_[Left_Foot].xpos;
        supportfoot_global_init_.linear() = rd_cc_.link_[Left_Foot].rotm;
    }

    // yaw only
    supportfoot_global_init_yaw_ = supportfoot_global_init_;
    supportfoot_global_init_yaw_.linear() = DyrosMath::rotateWithZ(DyrosMath::rot2Euler(supportfoot_global_init_.linear())(2));

    Eigen::Isometry3d ref_frame;
    ref_frame = supportfoot_global_init_yaw_;
    pelv_support_init_yaw_.translation() =DyrosMath::multiplyIsometry3dVector3d(DyrosMath::inverseIsometry3d(ref_frame) , rd_cc_.link_[Pelvis].xpos);
    pelv_support_init_yaw_.linear() = ref_frame.linear().transpose() *  rd_cc_.link_[Pelvis].rotm;
    com_support_init_yaw_ = DyrosMath::multiplyIsometry3dVector3d(DyrosMath::inverseIsometry3d(ref_frame), rd_cc_.link_[COM_id].xpos);
    com_support_init_dot_yaw_ = ref_frame.linear().transpose() * rd_cc_.link_[COM_id].v;
    pelv_support_euler_init_yaw_ = DyrosMath::rot2Euler(pelv_support_init_yaw_.linear());
    lfoot_global_init_.translation() = rd_cc_.link_[Left_Foot].xpos;
    lfoot_global_init_.linear() = rd_cc_.link_[Left_Foot].rotm;
    rfoot_global_init_.translation() = rd_cc_.link_[Right_Foot].xpos;
    rfoot_global_init_.linear() = rd_cc_.link_[Right_Foot].rotm;

    lfoot_support_init_yaw_ = DyrosMath::multiplyIsometry3d(DyrosMath::inverseIsometry3d(ref_frame), lfoot_global_init_);
    rfoot_support_init_yaw_ = DyrosMath::multiplyIsometry3d(DyrosMath::inverseIsometry3d(ref_frame), rfoot_global_init_);
    rfoot_support_euler_init_yaw_ = DyrosMath::rot2Euler(rfoot_support_init_yaw_.linear());
    lfoot_support_euler_init_yaw_ = DyrosMath::rot2Euler(lfoot_support_init_yaw_.linear());
}



// // Stepping stones version
void CustomController::updateFootstepCommand(){

    if (walking_tick == 0){

        // Compute Global Foot States, estimatesy
        lfoot_global_state.segment(0,2) = rd_cc_.link_[Left_Foot].xpos.segment(0,2);
        rfoot_global_state.segment(0,2) = rd_cc_.link_[Right_Foot].xpos.segment(0,2);
        lfoot_global_state(2) = DyrosMath::rot2Euler(rd_cc_.link_[Left_Foot].rotm)(2);
        rfoot_global_state(2) = DyrosMath::rot2Euler(rd_cc_.link_[Right_Foot].rotm)(2);
        lfoot_global_state(0) += 0.03;
        rfoot_global_state(0) += 0.03;
        for (int step = 0; step < number_of_foot_step; step++){
            if (step == 0) phase_indicator_(step) = first_stance_foot_;
            else phase_indicator_(step) = 1-phase_indicator_(step-1);

            step_length_x_(step) = phase_indicator_(step)*Lcommand_step_length_x_ + (1-phase_indicator_(step))*Rcommand_step_length_x_;
            step_length_y_(step) = (2*phase_indicator_(step) - 1) *(phase_indicator_(step)*Lcommand_step_length_y_ + (1-phase_indicator_(step))*Rcommand_step_length_y_);
            step_yaw_(step) = (2*phase_indicator_(step) - 1) * (phase_indicator_(step)*Lcommand_step_yaw_ + (1-phase_indicator_(step))*Rcommand_step_yaw_);
            foot_height_(step) = phase_indicator_(step)*Lcommand_foot_height_ + (1-phase_indicator_(step))*Rcommand_foot_height_;
            t_dsp_(step) = std::floor((phase_indicator_(step)*Lcommand_t_dsp_ + (1-phase_indicator_(step))*Rcommand_t_dsp_) * hz_);
            t_dsp_seconds(step) = phase_indicator_(step)*Lcommand_t_dsp_ + (1-phase_indicator_(step))*Rcommand_t_dsp_;
            t_ssp_(step) = std::floor((phase_indicator_(step)*Lcommand_t_ssp_ + (1-phase_indicator_(step))*Rcommand_t_ssp_ )* hz_);
            t_ssp_seconds(step) = phase_indicator_(step)*Lcommand_t_ssp_ + (1-phase_indicator_(step))*Rcommand_t_ssp_;
            t_total_(step) = 2*t_dsp_(step) + t_ssp_(step);

            
        }

        switch (ctrl_mode)
        {
        case 0:
            break;
        case 1: // Stepping stone
        {
            for (int step = 0; step < number_of_foot_step; step++){

                if (step == 0) phase_indicator_(step) = first_stance_foot_;
                else phase_indicator_(step) = 1-phase_indicator_(step-1);

                foot_height_(step) = foot_height_planned(step);
                t_dsp_(step) = std::floor(t_dsp_planned(step) * hz_);
                t_dsp_seconds(step) = t_dsp_planned(step);
                t_ssp_(step) = std::floor(t_ssp_planned(step) * hz_);
                t_ssp_seconds(step) = t_ssp_planned(step);
                t_total_(step) = 2*t_dsp_(step) + t_ssp_(step);
            }        
            Eigen::Vector3d &stance = (first_stance_foot_) ? rfoot_global_state : lfoot_global_state;
            Eigen::Vector3d &swing = (first_stance_foot_) ? lfoot_global_state : rfoot_global_state;
            double x_len = foothold_x_planned(0) - stance(0);
            double y_len = foothold_y_planned(0) - stance(1);
            step_length_x_(0) = cos(-stance(2))*x_len - sin(-stance(2))*y_len;
            step_length_y_(0) = sin(-stance(2))*x_len + cos(-stance(2))*y_len;
            step_yaw_(0) = foothold_yaw_planned(0) - stance(2);

            for (int step = 1; step < number_of_foot_step; step++){
                x_len = foothold_x_planned(step) - foothold_x_planned(step-1);
                y_len = foothold_y_planned(step) - foothold_y_planned(step-1);
                step_length_x_(step) = cos(-foothold_yaw_planned(step-1))*x_len - sin(-foothold_yaw_planned(step-1))*y_len;
                step_length_y_(step) = sin(-foothold_yaw_planned(step-1))*x_len + cos(-foothold_yaw_planned(step-1))*y_len;
                step_yaw_(step) = foothold_yaw_planned(step) - foothold_yaw_planned(step-1);

            }
            break;
        }

        case 2:
        {
            for (int step = 0; step < number_of_foot_step; step++){

                if (step == 0) phase_indicator_(step) = first_stance_foot_;
                else phase_indicator_(step) = 1-phase_indicator_(step-1);

                step_length_x_(step) = foothold_x_planned(step);
                step_length_y_(step) = (2*phase_indicator_(step) - 1) * foothold_y_planned(step);
                step_yaw_(step) = (2*phase_indicator_(step) - 1) * foothold_yaw_planned(step);
                foot_height_(step) = foot_height_planned(step);
                t_dsp_(step) = std::floor(t_dsp_planned(step) * hz_);
                t_dsp_seconds(step) = t_dsp_planned(step);
                t_ssp_(step) = std::floor(t_ssp_planned(step) * hz_);
                t_ssp_seconds(step) = t_ssp_planned(step);
                t_total_(step) = 2*t_dsp_(step) + t_ssp_(step);
            }
            break;
        }
        case 3:
        {
            for (int step = 0; step < number_of_foot_step; step++){

                if (step == 0) phase_indicator_(step) = first_stance_foot_;
                else phase_indicator_(step) = 1-phase_indicator_(step-1);

                step_length_x_(step) = step_length_x_fixed;
                step_length_y_(step) = (2*phase_indicator_(step) - 1) * step_length_y_fixed;
                step_yaw_(step) = (heading_error > 0) ? (phase_indicator_(step)) * 2 * heading_error : (1-phase_indicator_(step)) * 2 * heading_error;
                foot_height_(step) = foot_height_fixed;
                t_dsp_(step) = std::floor(t_dsp_fixed * hz_);
                t_dsp_seconds(step) = t_dsp_fixed;
                t_ssp_(step) = std::floor(t_ssp_fixed * hz_);
                t_ssp_seconds(step) = t_ssp_fixed;
                t_total_(step) = 2*t_dsp_(step) + t_ssp_(step);
            }
            break;
        }

        default:
            break;
        }

        calculateFootStepTotal();
        getRobotState();
        updateInitialState();
        getZmpTrajectory();
        if (ideal_preview)
            resetPreviewState();
    }

    else if (walking_tick > t_total_(0)){

        std::cout << "Foot Position error : " << sqrt(pow(swing_state_stance_frame_(0) - step_length_x_(0), 2) + pow(swing_state_stance_frame_(1) - step_length_y_(0), 2)) << " [m]" << std::endl;
        std::cout << ">> X error : " << sqrt(pow(swing_state_stance_frame_(0) - step_length_x_(0), 2)) << " [m]" << std::endl;
        std::cout << ">> Y error : " << sqrt(pow(swing_state_stance_frame_(1) - step_length_y_(0), 2)) << " [m]" << std::endl;
        Eigen::Quaterniond q1(swing_state_stance_frame_(6), swing_state_stance_frame_(3), swing_state_stance_frame_(4), swing_state_stance_frame_(5));
        double swing_yaw = std::atan2(2.0 * (q1.w() * q1.z() + q1.x() * q1.y()),
                             1.0 - 2.0 * (q1.y() * q1.y() + q1.z() * q1.z()));
        std::cout << "Foot Yaw error : " << sqrt(pow(CustomController::wrap_to_pi(swing_yaw - step_yaw_(0)), 2)) << " [rad]" << std::endl;
        

        x_error = (step_length_x_(0) - swing_state_stance_frame_(0)) ;
        y_error = (step_length_y_(0) - swing_state_stance_frame_(1)) ;
        yaw_error = (CustomController::wrap_to_pi(step_yaw_(0) - swing_yaw)) ;


        if (ctrl_mode){
            if (current_step_number < planned_step_number){
                evalFile << sqrt(pow(swing_state_stance_frame_(0) - step_length_x_(0), 2) + pow(swing_state_stance_frame_(1) - step_length_y_(0), 2)) << "\t";
                evalFile << sqrt(pow(swing_state_stance_frame_(0) - step_length_x_(0), 2)) << "\t";
                evalFile << sqrt(pow(swing_state_stance_frame_(1) - step_length_y_(0), 2)) << "\t";
                evalFile << sqrt(pow(CustomController::wrap_to_pi(swing_yaw - step_yaw_(0)), 2)) << "\t";
                evalFile << std::endl;

            }
        }

        step_length_x_.segment(0,number_of_foot_step-1) = step_length_x_.segment(1,number_of_foot_step-1);

        step_length_y_.segment(0,number_of_foot_step-1) = step_length_y_.segment(1,number_of_foot_step-1);

        step_yaw_.segment(0,number_of_foot_step-1) = step_yaw_.segment(1,number_of_foot_step-1);

        t_dsp_.segment(0,number_of_foot_step-1) = t_dsp_.segment(1,number_of_foot_step-1);
        t_dsp_seconds.segment(0,number_of_foot_step-1) = t_dsp_seconds.segment(1,number_of_foot_step-1);

        t_ssp_.segment(0,number_of_foot_step-1) = t_ssp_.segment(1,number_of_foot_step-1);
        t_ssp_seconds.segment(0,number_of_foot_step-1) = t_ssp_seconds.segment(1,number_of_foot_step-1);

        t_total_.segment(0, number_of_foot_step-1) = t_total_.segment(1, number_of_foot_step-1);

        foot_height_.segment(0,number_of_foot_step-1) = foot_height_.segment(1,number_of_foot_step-1);

        phase_indicator_.segment(0,number_of_foot_step-1) = phase_indicator_.segment(1,number_of_foot_step-1);

        phase_indicator_(number_of_foot_step-1) = 1-phase_indicator_(number_of_foot_step-2);

        int step = number_of_foot_step-1;
        step_length_x_(step) = phase_indicator_(step)*(Lcommand_step_length_x_) + (1-phase_indicator_(step))*(Rcommand_step_length_x_);
        step_length_y_(step) = (2*phase_indicator_(step) - 1) *(phase_indicator_(step)*(Lcommand_step_length_y_) + (1-phase_indicator_(step))*(Rcommand_step_length_y_));
        step_yaw_(step) = (2*phase_indicator_(step) - 1) * (phase_indicator_(step)*(Lcommand_step_yaw_) + (1-phase_indicator_(step))*(Rcommand_step_yaw_));
        foot_height_(step) = phase_indicator_(step)*Lcommand_foot_height_ + (1-phase_indicator_(step))*Rcommand_foot_height_;
        t_dsp_(step) = std::floor((phase_indicator_(step)*Lcommand_t_dsp_ + (1-phase_indicator_(step))*Rcommand_t_dsp_) * hz_);
        t_dsp_seconds(step) = phase_indicator_(step)*Lcommand_t_dsp_ + (1-phase_indicator_(step))*Rcommand_t_dsp_;
        t_ssp_(step) = std::floor((phase_indicator_(step)*Lcommand_t_ssp_ + (1-phase_indicator_(step))*Rcommand_t_ssp_ )* hz_);
        t_ssp_seconds(step) = phase_indicator_(step)*Lcommand_t_ssp_ + (1-phase_indicator_(step))*Rcommand_t_ssp_;
        t_total_(step) = 2*t_dsp_(step) + t_ssp_(step);


        walking_tick = 0;
        current_step_number++;

        switch (ctrl_mode)
        {
        case 0:
            break;
        case 1:
        {
            if (current_step_number < planned_step_number){

                Eigen::Vector3d &stance = (phase_indicator_(0)) ? rfoot_global_state : lfoot_global_state;
                Eigen::Vector3d &swing = (phase_indicator_(0)) ? lfoot_global_state : rfoot_global_state;
                double x_len = foothold_x_planned(current_step_number) - stance(0);
                double y_len = foothold_y_planned(current_step_number) - stance(1);
                
                step_length_x_(0) = cos(-stance(2))*x_len - sin(-stance(2))*y_len;
                step_length_y_(0) = sin(-stance(2))*x_len + cos(-stance(2))*y_len;
                step_yaw_(0) = foothold_yaw_planned(current_step_number) - stance(2);
                for (int step = 1; step < number_of_foot_step; step++){
                    if (step + current_step_number < planned_step_number){
                        x_len = foothold_x_planned(step + current_step_number) - foothold_x_planned(step + current_step_number-1);
                        y_len = foothold_y_planned(step + current_step_number) - foothold_y_planned(step + current_step_number-1);
                        step_length_x_(step) = cos(-foothold_yaw_planned(step + current_step_number-1))*x_len - sin(-foothold_yaw_planned(step + current_step_number-1))*y_len;
                        step_length_y_(step) = sin(-foothold_yaw_planned(step + current_step_number-1))*x_len + cos(-foothold_yaw_planned(step + current_step_number-1))*y_len;
                        step_yaw_(step) = foothold_yaw_planned(step + current_step_number) - foothold_yaw_planned(step + current_step_number-1);

                        foot_height_(step) = foot_height_planned(step + current_step_number);
                        t_dsp_(step) = std::floor(t_dsp_planned(step + current_step_number) * hz_);
                        t_dsp_seconds(step) = t_dsp_planned(step + current_step_number);
                        t_ssp_(step) = std::floor(t_ssp_planned(step + current_step_number) * hz_);
                        t_ssp_seconds(step) = t_ssp_planned(step + current_step_number);
                        t_total_(step) = 2*t_dsp_(step) + t_ssp_(step);
                    }
                    else{
                        step_length_x_(step) = 0.;
                        step_length_y_(step) = (2*phase_indicator_(step) - 1) *0.21;
                        step_yaw_(step) = 0.;
                        foot_height_(step) = 0.1;
                        t_dsp_(step) = std::floor(0.1* hz_);
                        t_dsp_seconds(step) = 0.1;
                        t_ssp_(step) = std::floor(0.7 * hz_);
                        t_ssp_seconds(step) =0.7;
                        t_total_(step) = 2*t_dsp_(step) + t_ssp_(step);

                    }
                }
            }
            else {
                step_length_x_(step) = 0.;
                step_length_y_(step) = (2*phase_indicator_(step) - 1) *0.21;
                step_yaw_(step) = 0.;
                foot_height_(step) = 0.1;
                t_dsp_(step) = std::floor(0.1* hz_);
                t_dsp_seconds(step) = 0.1;
                t_ssp_(step) = std::floor(0.7 * hz_);
                t_ssp_seconds(step) =0.7;
                t_total_(step) = 2*t_dsp_(step) + t_ssp_(step);

            }
            break;            
        }

        case 2:
        {
            if (step + current_step_number < planned_step_number){
                step_length_x_(step) = foothold_x_planned(step + current_step_number);
                step_length_y_(step) = (2*phase_indicator_(step) - 1) * foothold_y_planned(step + current_step_number);
                step_yaw_(step) = (2*phase_indicator_(step) - 1) * (foothold_yaw_planned(step + current_step_number) );
                foot_height_(step) = foot_height_planned(step + current_step_number);
                t_dsp_(step) = std::floor(t_dsp_planned(step + current_step_number) * hz_);
                t_dsp_seconds(step) = t_dsp_planned(step + current_step_number);
                t_ssp_(step) = std::floor(t_ssp_planned(step + current_step_number) * hz_);
                t_ssp_seconds(step) = t_ssp_planned(step + current_step_number);
                t_total_(step) = 2*t_dsp_(step) + t_ssp_(step);

            }
            else{
                step_length_x_(step) = 0.;
                step_length_y_(step) = (2*phase_indicator_(step) - 1) *0.21;
                step_yaw_(step) = 0.;
                foot_height_(step) = 0.1;
                t_dsp_(step) = std::floor(0.1* hz_);
                t_dsp_seconds(step) = 0.1;
                t_ssp_(step) = std::floor(0.7 * hz_);
                t_ssp_seconds(step) =0.7;
                t_total_(step) = 2*t_dsp_(step) + t_ssp_(step);

            }
            break;
        }

        case 3:
        {
            step_length_x_(step) = step_length_x_fixed;
            step_length_y_(step) = (2*phase_indicator_(step) - 1) * step_length_y_fixed;
            step_yaw_(step) = (heading_error > 0) ? (phase_indicator_(step)) * 2 * heading_error : (1-phase_indicator_(step)) * 2 * heading_error;
            foot_height_(step) = foot_height_fixed;
            t_dsp_(step) = std::floor(t_dsp_fixed * hz_);
            t_dsp_seconds(step) = t_dsp_fixed;
            t_ssp_(step) = std::floor(t_ssp_fixed * hz_);
            t_ssp_seconds(step) = t_ssp_fixed;
            t_total_(step) = 2*t_dsp_(step) + t_ssp_(step);
            break;
        }

        default:
            break;
        }
        calculateFootStepTotal();
        getRobotState();
        updateInitialState();
        getZmpTrajectory();
        if (ideal_preview)
            resetPreviewState();
    }

}


void CustomController::getRobotState()
{

    pelv_rpy_current_.setZero();
    pelv_rpy_current_ = DyrosMath::rot2Euler(rd_cc_.link_[Pelvis].rotm); //ZYX multiply

    R_angle = pelv_rpy_current_(0);
    P_angle = pelv_rpy_current_(1);
    pelv_yaw_rot_current_from_global_ = DyrosMath::rotateWithZ(pelv_rpy_current_(2));
    pelv_yaw_rot_current_from_global_.linear() = rd_cc_.link_[Pelvis].rotm;

    rfoot_rpy_current_.setZero();
    lfoot_rpy_current_.setZero();
    rfoot_rpy_current_ = DyrosMath::rot2Euler(rd_cc_.link_[Right_Foot].rotm);
    lfoot_rpy_current_ = DyrosMath::rot2Euler(rd_cc_.link_[Left_Foot].rotm);

    rfoot_roll_rot_ = DyrosMath::rotateWithX(rfoot_rpy_current_(0));
    lfoot_roll_rot_ = DyrosMath::rotateWithX(lfoot_rpy_current_(0));
    rfoot_pitch_rot_ = DyrosMath::rotateWithY(rfoot_rpy_current_(1));
    lfoot_pitch_rot_ = DyrosMath::rotateWithY(lfoot_rpy_current_(1));

    pelv_global_current_.translation() = rd_cc_.link_[Pelvis].xpos;
    pelv_global_current_.linear() = rd_cc_.link_[Pelvis].rotm;
    lfoot_global_current_.translation() = rd_cc_.link_[Left_Foot].xpos;
    lfoot_global_current_.linear() = rd_cc_.link_[Left_Foot].rotm;
    rfoot_global_current_.translation() = rd_cc_.link_[Right_Foot].xpos;
    rfoot_global_current_.linear() = rd_cc_.link_[Right_Foot].rotm;
    com_global_current_ = rd_cc_.link_[COM_id].xpos;
    com_global_current_dot_prev_ = com_global_current_dot_;
    com_global_current_dot_ = rd_cc_.link_[COM_id].v;

    double support_foot_flag = foot_step_(0, 6);
    if (support_foot_flag == 0)
    {
        supportfoot_global_current_.translation() = rd_cc_.link_[Right_Foot].xpos;
        supportfoot_global_current_.linear() = DyrosMath::rotateWithZ(DyrosMath::rot2Euler(rd_cc_.link_[Right_Foot].rotm)(2));
    }
    else if (support_foot_flag == 1)
    {
        supportfoot_global_current_.translation() = rd_cc_.link_[Left_Foot].xpos;
        supportfoot_global_current_.linear() = DyrosMath::rotateWithZ(DyrosMath::rot2Euler(rd_cc_.link_[Left_Foot].rotm)(2));
    }

    pelv_support_current_ = DyrosMath::inverseIsometry3d(supportfoot_global_current_) * pelv_global_current_;
    lfoot_support_current_ = DyrosMath::inverseIsometry3d(supportfoot_global_current_) * lfoot_global_current_;
    rfoot_support_current_ = DyrosMath::inverseIsometry3d(supportfoot_global_current_) * rfoot_global_current_;

    com_support_current_ = DyrosMath::multiplyIsometry3dVector3d(DyrosMath::inverseIsometry3d(supportfoot_global_current_), com_global_current_);
    com_support_current_dot_ = DyrosMath::multiplyIsometry3dVector3d(DyrosMath::inverseIsometry3d(supportfoot_global_current_), com_global_current_dot_);
    // std::cout << "Support foot is : " << ((phase_indicator_(0)) ? "right" : "left") << std::endl;
    // std::cout << "support foot global pos : " << supportfoot_global_init_.translation().transpose() << std::endl;
    // std::cout << "Left foot global pos : " << rd_.link_[Left_Foot].xpos.transpose() << std::endl;
    // std::cout << "Right foot global pos : " << rd_.link_[Right_Foot].xpos.transpose() << std::endl;
    swing_state_stance_frame_.segment(0,3) = DyrosMath::multiplyIsometry3dVector3d(DyrosMath::inverseIsometry3d(supportfoot_global_current_), phase_indicator_(0)*rd_cc_.link_[Left_Foot].xpos + (1-phase_indicator_(0))*rd_cc_.link_[Right_Foot].xpos); // Compute swing foot state from init support foot
    Eigen::Quaterniond swing_quat(phase_indicator_(0)*(supportfoot_global_current_.linear().transpose() * rd_cc_.link_[Left_Foot].rotm) + (1-phase_indicator_(0))*(supportfoot_global_current_.linear().transpose() * rd_cc_.link_[Right_Foot].rotm));
    swing_state_stance_frame_(3) = swing_quat.x();
    swing_state_stance_frame_(4) = swing_quat.y();
    swing_state_stance_frame_(5) = swing_quat.z();
    swing_state_stance_frame_(6) = swing_quat.w();

    com_state_stance_frame_.segment(0,3) = DyrosMath::multiplyIsometry3dVector3d(DyrosMath::inverseIsometry3d(supportfoot_global_current_), rd_cc_.link_[COM_id].xpos); // Compute swing foot state from init support foot
    Eigen::Quaterniond com_quat(supportfoot_global_current_.linear().transpose() * rd_cc_.link_[COM_id].rotm);
    com_state_stance_frame_(3) = com_quat.x();
    com_state_stance_frame_(4) = com_quat.y();
    com_state_stance_frame_(5) = com_quat.z();
    com_state_stance_frame_(6) = com_quat.w();

    if (!ideal_preview){
        x_preview_(0) = com_support_current_(0);
        y_preview_(0) = com_support_current_(1);
        
        x_preview_(1) = com_support_current_dot_(0);
        y_preview_(1) = com_support_current_dot_(1);
    
        x_preview_(2) = DyrosMath::multiplyIsometry3dVector3d(DyrosMath::inverseIsometry3d(supportfoot_global_current_), com_global_current_dot_ - com_global_current_dot_prev_)(0) * hz_;
        y_preview_(2) = DyrosMath::multiplyIsometry3dVector3d(DyrosMath::inverseIsometry3d(supportfoot_global_current_), com_global_current_dot_ - com_global_current_dot_prev_)(1) * hz_;
    }

    // Compute Global Foot States, estimates

    Eigen::Vector3d &stance = (phase_indicator_(0)) ? rfoot_global_state : lfoot_global_state;
    Eigen::Vector3d &swing = (phase_indicator_(0)) ? lfoot_global_state : rfoot_global_state;
    Eigen::Isometry3d &swing_stance = (phase_indicator_(0)) ? lfoot_support_current_ : rfoot_support_current_;
    
    double swing_yaw_stance = DyrosMath::rot2Euler(swing_stance.linear())(2);
    swing(0) = stance(0) + cos(stance(2))*swing_stance.translation()(0) - sin(stance(2))*swing_stance.translation()(1);
    swing(1) = stance(1) + sin(stance(2))*swing_stance.translation()(0) + cos(stance(2))*swing_stance.translation()(1);
    swing(2) = stance(2) + swing_yaw_stance;

  
}

void CustomController::calculateFootStepTotal()
{

    foot_step_.resize(number_of_foot_step, 7);
    foot_step_.setZero();
    foot_step_support_frame_.resize(number_of_foot_step, 7);
    foot_step_support_frame_.setZero();
    
    // foot_step_ is foothold command in that step's stance foot
    // foot_step_support_frame_ is foothold command in the first stance foot frame
    foot_step_(0,0) = step_length_x_(0);
    foot_step_(0,1) = step_length_y_(0);
    foot_step_(0,5) = step_yaw_(0);
    foot_step_(0,6) = 1-phase_indicator_(0);
    foot_step_support_frame_(0, 0) = step_length_x_(0);
    foot_step_support_frame_(0, 1) = step_length_y_(0);
    foot_step_support_frame_(0, 5) = step_yaw_(0);
    foot_step_support_frame_(0, 6) = 1-phase_indicator_(0);

    for (int i = 1; i < number_of_foot_step; i++){
        foot_step_(i,0) = step_length_x_(i);
        foot_step_(i,1) = step_length_y_(i);
        foot_step_(i,5) = step_yaw_(i);
        foot_step_(i,6) = 1-phase_indicator_(i);
        
        foot_step_support_frame_(i, 0) = foot_step_support_frame_(i-1, 0) + cos(foot_step_support_frame_(i-1, 5)) * step_length_x_(i) - sin(foot_step_support_frame_(i-1, 5)) * step_length_y_(i);
        foot_step_support_frame_(i, 1) = foot_step_support_frame_(i-1, 1) + sin(foot_step_support_frame_(i-1, 5)) * step_length_x_(i) + cos(foot_step_support_frame_(i-1, 5)) * step_length_y_(i);
        foot_step_support_frame_(i, 5) = foot_step_support_frame_(i-1, 5) + step_yaw_(i);
        foot_step_support_frame_(i, 6) = 1-phase_indicator_(i);
    }
}


void CustomController::addZmpOffset()
{

    foot_step_support_frame_offset_ = foot_step_support_frame_;
    foot_step_offset_ = foot_step_;
    // ZMP OFFSET
    for (int i = 0; i < number_of_foot_step; i++){
        foot_step_offset_(i, 0) += phase_indicator_(i)*zmp_offset*sin(foot_step_(i, 5)) - (1 - phase_indicator_(i))*zmp_offset*sin(foot_step_(i, 5));
        foot_step_offset_(i, 1) += -phase_indicator_(i)*zmp_offset*cos(foot_step_(i, 5)) + (1 - phase_indicator_(i))*zmp_offset*cos(foot_step_(i, 5));
        foot_step_offset_(i, 0) += zmp_offset_x*cos(foot_step_(i, 5));
        foot_step_offset_(i, 1) += zmp_offset_x*sin(foot_step_(i, 5));

        foot_step_support_frame_offset_(i, 0) += phase_indicator_(i)*zmp_offset*sin(foot_step_support_frame_(i, 5)) - (1 - phase_indicator_(i))*zmp_offset*sin(foot_step_support_frame_(i, 5));
        foot_step_support_frame_offset_(i, 1) += -phase_indicator_(i)*zmp_offset*cos(foot_step_support_frame_(i, 5)) + (1 - phase_indicator_(i))*zmp_offset*cos(foot_step_support_frame_(i, 5));
        foot_step_support_frame_offset_(i, 0) += zmp_offset_x*cos(foot_step_support_frame_(i, 5));
        foot_step_support_frame_offset_(i, 1) += zmp_offset_x*sin(foot_step_support_frame_(i, 5));
    }
}


void CustomController::getZmpTrajectory()
{
    unsigned int norm_size = 0;

    norm_size = 4.0*hz_ ; // compute zmp over the three planned steps
    addZmpOffset(); 

    zmpGenerator(norm_size);

}


void CustomController::zmpGenerator(const unsigned int norm_size)
{
    /*
    Goal
    ----------
    -> To position the ZMP at the center of the foot sole according to the footstep planning to prevent the robot from falling during walking.

    Parameters
    ----------
    -> norm_size : The size of the previewed vector for the ZMP reference.

    -> planning_step_num : The number of footsteps to be predicted.

    Returns
    -------
    -> ref_zmp_ : The ZMP reference vector calculated based on the foot sole plan (size: 2 x norm_size).
    */

    ref_zmp_.setZero(norm_size, 2);
    ref_zmp_thread3.setZero(norm_size, 2);
    ref_com_yaw_.setZero(norm_size);
    ref_com_yawvel_.setZero(norm_size);

    Eigen::VectorXd temp_px;
    Eigen::VectorXd temp_py;
    Eigen::VectorXd temp_yaw;
    Eigen::VectorXd temp_yawvel;

    unsigned int index = 0;

  
    for (unsigned int i = 0; i < number_of_foot_step; i++)
    {   
        onestepZmp(i, temp_px, temp_py, temp_yaw, temp_yawvel); // save 1-step zmp into temp px, py
        ref_zmp_.block(index, 0, t_total_(i), 1) = temp_px; 
        ref_zmp_.block(index, 1, t_total_(i), 1) = temp_py;
        ref_com_yaw_.segment(index, t_total_(i)) = temp_yaw;
        ref_com_yawvel_.segment(index, t_total_(i)) = temp_yawvel;

        index = index + t_total_(i);                                                          
    }
    if (t_total_.sum() < norm_size){
        for (int i = t_total_.sum(); i < norm_size; i++){
            ref_zmp_(i, 0) = ref_zmp_(i-1, 0);
            ref_zmp_(i, 1) = ref_zmp_(i-1, 1);
            ref_com_yaw_(i) = ref_com_yaw_(i-1);
            ref_com_yawvel_(i) = ref_com_yawvel_(i-1);
        }
    }

}

void CustomController::onestepZmp(unsigned int current_step_number, Eigen::VectorXd &temp_px, Eigen::VectorXd &temp_py, Eigen::VectorXd &temp_yaw, Eigen::VectorXd &temp_yawvel) // CoM Yaw as well.
{
    temp_px.setZero(t_total_(current_step_number));  
    temp_py.setZero(t_total_(current_step_number));
    temp_yaw.setZero(t_total_(current_step_number));
    temp_yawvel.setZero(t_total_(current_step_number));

    double v0_x_dsp1 = 0.0; double v0_y_dsp1 = 0.0;
    double vT_x_dsp1 = 0.0; double vT_y_dsp1 = 0.0;
    double v0_yaw_dsp1 = 0.0; double vT_yaw_dsp1 = 0.0;
    double v0_x_ssp  = 0.0; double v0_y_ssp = 0.0;
    double vT_x_ssp  = 0.0; double vT_y_ssp = 0.0;
    double v0_yaw_ssp  = 0.0; double vT_yaw_ssp = 0.0;
    double v0_x_dsp2 = 0.0; double v0_y_dsp2 = 0.0;
    double vT_x_dsp2 = 0.0; double vT_y_dsp2 = 0.0;
    double v0_yaw_dsp2 = 0.0; double vT_yaw_dsp2 = 0.0;

    double t_dsp1_ = t_dsp_(current_step_number);
    double t_dsp2_ = t_dsp_(current_step_number);
    double t_ssp = t_ssp_(current_step_number);
    double t_total = t_total_(current_step_number);    
    
    //TODO CoM Yaw implement
    if (current_step_number == 0)
    {
        // v0_x_dsp1 = com_support_init_yaw_(0);
        v0_x_dsp1 = (phase_indicator_(0)*lfoot_support_init_yaw_.translation()(0) + (1-phase_indicator_(0))*rfoot_support_init_yaw_.translation()(0))/2;
        vT_x_dsp1 = zmp_offset_x;
        // v0_y_dsp1 = com_support_init_yaw_(1);
        v0_y_dsp1 = (phase_indicator_(0)*lfoot_support_init_yaw_.translation()(1) + (1-phase_indicator_(0))*rfoot_support_init_yaw_.translation()(1))/2;
        vT_y_dsp1 = (2*phase_indicator_(current_step_number)-1)*zmp_offset;
        // v0_yaw_dsp1 =  pelv_support_euler_init_yaw_(2);
        v0_yaw_dsp1 = DyrosMath::rot2Euler(phase_indicator_(0)*lfoot_support_init_yaw_.linear() + (1-phase_indicator_(0))*rfoot_support_init_yaw_.linear())(2)/2;
        vT_yaw_dsp1 = v0_yaw_dsp1;

        v0_x_ssp = vT_x_dsp1;
        vT_x_ssp = v0_x_ssp;
        v0_y_ssp = vT_y_dsp1;
        vT_y_ssp = v0_y_ssp;
        v0_yaw_ssp =  vT_yaw_dsp1;
        vT_yaw_ssp = foot_step_support_frame_offset_(current_step_number - 0, 5) / 2.0;

        v0_x_dsp2 = vT_x_ssp;
        vT_x_dsp2 = zmp_offset_x + (foot_step_support_frame_offset_(current_step_number - 0, 0)) / 2.0;
        v0_y_dsp2 = vT_y_ssp;
        vT_y_dsp2 = ((2*phase_indicator_(current_step_number - 0)-1)*zmp_offset + foot_step_support_frame_offset_(current_step_number - 0, 1)) / 2.0;
        v0_yaw_dsp2 = vT_yaw_ssp;
        vT_yaw_dsp2 = v0_yaw_dsp2;


    }
    else if (current_step_number == 1)
    { 
        v0_x_dsp1 = zmp_offset_x + (foot_step_support_frame_offset_(current_step_number - 1, 0) ) / 2.0;
        vT_x_dsp1 =  foot_step_support_frame_offset_(current_step_number - 1, 0);
        v0_y_dsp1 = ((2*phase_indicator_(current_step_number - 1)-1)*zmp_offset + foot_step_support_frame_offset_(current_step_number - 1, 1)) / 2.0;
        vT_y_dsp1 =  foot_step_support_frame_offset_(current_step_number - 1, 1);
        v0_yaw_dsp1 = foot_step_support_frame_offset_(current_step_number - 1, 5) / 2.0;
        vT_yaw_dsp1 = v0_yaw_dsp1;

        v0_x_ssp = vT_x_dsp1;
        vT_x_ssp = v0_x_ssp;
        v0_y_ssp = vT_y_dsp1;
        vT_y_ssp = v0_y_ssp;
        v0_yaw_ssp = (foot_step_support_frame_offset_(current_step_number - 1, 5) )/ 2.0;
        vT_yaw_ssp = (foot_step_support_frame_offset_(current_step_number - 0, 5) + foot_step_support_frame_offset_(current_step_number - 1, 5))/ 2;

        v0_x_dsp2 =  vT_x_ssp;
        vT_x_dsp2 = (foot_step_support_frame_offset_(current_step_number, 0) + foot_step_support_frame_offset_(current_step_number - 1, 0)) / 2.0;
        v0_y_dsp2 = vT_y_ssp;
        vT_y_dsp2 = (foot_step_support_frame_offset_(current_step_number, 1) + foot_step_support_frame_offset_(current_step_number - 1, 1)) / 2.0;
        v0_yaw_dsp2 = (foot_step_support_frame_offset_(current_step_number - 0, 5) + foot_step_support_frame_offset_(current_step_number - 1, 5))/ 2;
        vT_yaw_dsp2 = (foot_step_support_frame_offset_(current_step_number - 0, 5) + foot_step_support_frame_offset_(current_step_number - 1, 5))/ 2;

    }
    else
    {   
        v0_x_dsp1 = (foot_step_support_frame_offset_(current_step_number - 2, 0) + foot_step_support_frame_offset_(current_step_number - 1, 0)) / 2.0;
        vT_x_dsp1 =  foot_step_support_frame_offset_(current_step_number - 1, 0);
        v0_y_dsp1 = (foot_step_support_frame_offset_(current_step_number - 2, 1) + foot_step_support_frame_offset_(current_step_number - 1, 1)) / 2.0;
        vT_y_dsp1 =  foot_step_support_frame_offset_(current_step_number - 1, 1);
        v0_yaw_dsp1 = (foot_step_support_frame_offset_(current_step_number - 2, 5) + foot_step_support_frame_offset_(current_step_number - 1, 5))/ 2;
        vT_yaw_dsp1 = (foot_step_support_frame_offset_(current_step_number - 2, 5) + foot_step_support_frame_offset_(current_step_number - 1, 5))/ 2;

        v0_x_ssp = foot_step_support_frame_offset_(current_step_number - 1, 0);
        vT_x_ssp = foot_step_support_frame_offset_(current_step_number - 1, 0);
        v0_y_ssp = foot_step_support_frame_offset_(current_step_number - 1, 1);
        vT_y_ssp = foot_step_support_frame_offset_(current_step_number - 1, 1);
        v0_yaw_ssp =  (foot_step_support_frame_offset_(current_step_number - 2, 5) + foot_step_support_frame_offset_(current_step_number - 1, 5))/ 2;
        vT_yaw_ssp = (foot_step_support_frame_offset_(current_step_number - 0, 5) + foot_step_support_frame_offset_(current_step_number - 1, 5))/ 2;

        v0_x_dsp2 =  foot_step_support_frame_offset_(current_step_number - 1, 0);
        vT_x_dsp2 = (foot_step_support_frame_offset_(current_step_number - 1, 0) + foot_step_support_frame_offset_(current_step_number - 0, 0)) / 2.0;
        v0_y_dsp2 =  foot_step_support_frame_offset_(current_step_number - 1, 1);
        vT_y_dsp2 = (foot_step_support_frame_offset_(current_step_number - 1, 1) + foot_step_support_frame_offset_(current_step_number - 0, 1)) / 2.0;
        v0_yaw_dsp2 = (foot_step_support_frame_offset_(current_step_number - 0, 5) + foot_step_support_frame_offset_(current_step_number - 1, 5))/ 2;
        vT_yaw_dsp2 = (foot_step_support_frame_offset_(current_step_number - 0, 5) + foot_step_support_frame_offset_(current_step_number - 1, 5))/ 2;
    }

    double lin_interpol = 0.0;
    for (int i = 0; i < t_total; i++)
    {
        if (i < t_dsp1_) 
        { 
            // lin_interpol = i / (t_dsp1_);
            // temp_px(i) = (1.0 - lin_interpol) * v0_x_dsp1 + lin_interpol * vT_x_dsp1;
            // temp_py(i) = (1.0 - lin_interpol) * v0_y_dsp1 + lin_interpol * vT_y_dsp1;
            
            temp_px(i) = DyrosMath::minmax_cut(DyrosMath::cubic(i, 0.0, t_dsp1_, v0_x_dsp1, vT_x_dsp1, 0.0, 0.0), min(v0_x_dsp1, vT_x_dsp1), max(v0_x_dsp1, vT_x_dsp1));
            temp_py(i) = DyrosMath::minmax_cut(DyrosMath::cubic(i, 0.0, t_dsp1_, v0_y_dsp1, vT_y_dsp1, 0.0, 0.0), min(v0_y_dsp1, vT_y_dsp1), max(v0_y_dsp1, vT_y_dsp1));
            temp_yaw(i) = DyrosMath::minmax_cut(DyrosMath::cubic(i, 0.0, t_dsp1_, v0_yaw_dsp1, vT_yaw_dsp1, 0.0, 0.0), min(v0_yaw_dsp1, vT_yaw_dsp1), max(v0_yaw_dsp1, vT_yaw_dsp1));
            temp_yawvel(i) = DyrosMath::cubicDot(i, 0.0, t_dsp1_, v0_yaw_dsp1, vT_yaw_dsp1, 0., 0.);
        }
        else if (i >= t_dsp1_ && i < t_dsp1_ + t_ssp)
        {
            // lin_interpol = (i - t_dsp1_) / t_ssp_;
            // temp_px(i) = (1.0 - lin_interpol) * v0_x_ssp + lin_interpol * vT_x_ssp;
            // temp_py(i) = (1.0 - lin_interpol) * v0_y_ssp + lin_interpol * vT_y_ssp;

            temp_px(i) = DyrosMath::minmax_cut(DyrosMath::cubic(i, t_dsp1_, t_dsp1_ + t_ssp, v0_x_ssp, vT_x_ssp, 0.0, 0.0), min(v0_x_ssp, vT_x_ssp), max(v0_x_ssp, vT_x_ssp));
            temp_py(i) = DyrosMath::minmax_cut(DyrosMath::cubic(i, t_dsp1_, t_dsp1_ + t_ssp, v0_y_ssp, vT_y_ssp, 0.0, 0.0), min(v0_y_ssp, vT_y_ssp), max(v0_y_ssp, vT_y_ssp));
            temp_yaw(i) = DyrosMath::minmax_cut(DyrosMath::cubic(i, t_dsp1_, t_dsp1_ + t_ssp, v0_yaw_ssp, vT_yaw_ssp, 0.0, 0.0), min(v0_yaw_ssp, vT_yaw_ssp), max(v0_yaw_ssp, vT_yaw_ssp));
            temp_yawvel(i) = DyrosMath::cubicDot(i, t_dsp1_, t_dsp1_+t_ssp, v0_yaw_ssp, vT_yaw_ssp, 0., 0.);
        }
        else
        {
            // lin_interpol = (i - t_dsp1_ - t_ssp_) / t_dsp2_;
            // temp_px(i) = (1.0 - lin_interpol) * v0_x_dsp2 + lin_interpol * vT_x_dsp2;
            // temp_py(i) = (1.0 - lin_interpol) * v0_y_dsp2 + lin_interpol * vT_y_dsp2;
            temp_px(i) = DyrosMath::minmax_cut(DyrosMath::cubic(i, t_dsp1_ + t_ssp, t_total, v0_x_dsp2, vT_x_dsp2, 0.0, 0.0), min(v0_x_dsp2, vT_x_dsp2), max(v0_x_dsp2, vT_x_dsp2));
            temp_py(i) = DyrosMath::minmax_cut(DyrosMath::cubic(i, t_dsp1_ + t_ssp, t_total, v0_y_dsp2, vT_y_dsp2, 0.0, 0.0), min(v0_y_dsp2, vT_y_dsp2), max(v0_y_dsp2, vT_y_dsp2));
            temp_yaw(i) = DyrosMath::minmax_cut(DyrosMath::cubic(i, t_dsp1_ + t_ssp, t_total, v0_yaw_dsp2, vT_yaw_dsp2, 0.0, 0.0), min(v0_yaw_dsp2, vT_yaw_dsp2), max(v0_yaw_dsp2, vT_yaw_dsp2));
            temp_yawvel(i) = DyrosMath::cubicDot(i, t_dsp1_ + t_ssp, t_total, v0_yaw_dsp2, vT_yaw_dsp2, 0., 0.);
        }
        // std::cout << current_step_number << " step's temp_py " << i << " : " << temp_py(i) << std::endl;
    }

}

void CustomController::resetPreviewState(){
    x_preview_.setZero(); y_preview_.setZero(); 
    x_preview_(0) = com_support_current_(0);
    y_preview_(0) = com_support_current_(1);
    
    x_preview_(1) = com_support_current_dot_(0);
    y_preview_(1) = com_support_current_dot_(1);

    x_preview_(2) = DyrosMath::multiplyIsometry3dVector3d(DyrosMath::inverseIsometry3d(supportfoot_global_current_), com_global_current_dot_ - com_global_current_dot_prev_)(0) * hz_;
    y_preview_(2) = DyrosMath::multiplyIsometry3dVector3d(DyrosMath::inverseIsometry3d(supportfoot_global_current_), com_global_current_dot_ - com_global_current_dot_prev_)(1) * hz_;
    UX_preview_ = 0;
    UY_preview_ = 0;
    windupPreview();
}

void CustomController::getComTrajectory()
{
    double dt_preview_ = 1.0 / hz_; // : sampling time of preview [s]
    double NL_preview  = 1.6 * hz_;      // : number of preview horizons

    if (is_preview_ctrl_init == true)
    {
        Gi_preview_.setZero();
        Gd_preview_.setZero();
        Gx_preview_.setZero();

        preview_Parameter(dt_preview_, NL_preview, Gi_preview_, Gd_preview_, Gx_preview_, A_preview_, B_preview_, C_preview_);
        
        resetPreviewState();

        is_preview_ctrl_init = false;

        std::cout << "PREVIEW PARAMETERS ARE SUCCESSFULLY INITIALIZED" << std::endl;
    }

    previewcontroller(dt_preview_, NL_preview, walking_tick, 
                      x_preview_, y_preview_, UX_preview_, UY_preview_,
                      Gi_preview_, Gd_preview_, Gx_preview_, 
                      A_preview_, B_preview_, C_preview_);

    com_desired_(0) = x_preview_(0);
    com_desired_(1) = y_preview_(0);

    com_desired_(2) = com_height_;
    com_desired_dot_(0) = x_preview_(1);
    com_desired_dot_(1) = y_preview_(1);
    com_desired_dot_(2) = 0.;

    if (!ideal_preview) windupPreview();

    
}

void CustomController::preview_Parameter(double dt, int NL, Eigen::MatrixXd &Gi, Eigen::VectorXd &Gd, Eigen::MatrixXd &Gx, Eigen::MatrixXd &A, Eigen::VectorXd &B, Eigen::MatrixXd &C)
{
    A.resize(3, 3);
    A(0, 0) = 1.0;
    A(0, 1) = dt;
    A(0, 2) = dt * dt * 0.5;
    A(1, 0) = 0;
    A(1, 1) = 1.0;
    A(1, 2) = dt;
    A(2, 0) = 0;
    A(2, 1) = 0;
    A(2, 2) = 1;

    B.resize(3);
    B(0) = dt * dt * dt / 6;
    B(1) = dt * dt / 2;
    B(2) = dt;

    C.resize(1, 3);
    C(0, 0) = 1;
    C(0, 1) = 0;
    C(0, 2) = -com_height_ / GRAVITY;

    Eigen::MatrixXd A_bar;
    Eigen::VectorXd B_bar;

    B_bar.setZero(4);
    B_bar.segment(0, 1) = C * B;
    B_bar.segment(1, 3) = B;

    Eigen::Matrix1x4d B_bar_tran;
    B_bar_tran = B_bar.transpose();

    Eigen::MatrixXd I_bar;
    Eigen::MatrixXd F_bar;
    A_bar.setZero(4, 4);
    I_bar.setZero(4, 1);
    F_bar.setZero(4, 3);

    F_bar.block<1, 3>(0, 0) = C * A;
    F_bar.block<3, 3>(1, 0) = A;

    I_bar.setZero();
    I_bar(0, 0) = 1.0;

    A_bar.block<4, 1>(0, 0) = I_bar;
    A_bar.block<4, 3>(0, 1) = F_bar;

    Eigen::MatrixXd Qe;
    Qe.setZero(1, 1);
    Qe(0, 0) = 1.0;

    Eigen::MatrixXd R;
    R.setZero(1, 1);
    R(0, 0) = 0.000001;

    Eigen::MatrixXd Qx;
    Qx.setZero(3, 3);

    Eigen::MatrixXd Q_bar;
    Q_bar.setZero(3, 3);
    Q_bar(0, 0) = Qe(0, 0);

    Eigen::Matrix4d K; K.setZero();
    if (com_height_ == 0.7){
        K(0, 0) = 68.7921510770868;
        K(0, 1) = 2331.78394937073;
        K(0, 2) = 632.495145628673;
        K(0, 3) = 2.57540412848618;
        K(1, 0) = 2331.78394937073;
        K(1, 1) = 81346.5405208585;
        K(1, 2) = 22074.2517082842;
        K(1, 3) = 92.2651171309122;
        K(2, 0) = 632.495145628673;
        K(2, 1) = 22074.2517082842;
        K(2, 2) = 5990.22394297225;
        K(2, 3) = 25.0754425025208;
        K(3, 0) = 2.57540412848618;
        K(3, 1) = 92.2651171309122;
        K(3, 2) = 25.0754425025208;
        K(3, 3) = 0.115349533788953;
    }

    else if (com_height_ == 0.728){
        K(0, 0) = 70.0813009391222;
        K(0, 1) = 2420.65372019101;
        K(0, 2) = 669.387885399937;
        K(0, 3) = 2.72127899669319;
        K(1, 0) = 2420.65372019101;
        K(1, 1) = 85969.0761591658;
        K(1, 2) = 23782.1340369413;
        K(1, 3) = 99.0855302680049;
        K(2, 0) = 669.387885399937;
        K(2, 1) = 23782.1340369413;
        K(2, 2) = 6579.12894080612;
        K(2, 3) = 27.4486400273521;
        K(3, 0) = 2.72127899669319;
        K(3, 1) = 99.0855302680049;
        K(3, 2) = 27.4486400273521;
        K(3, 3) = 0.125016968799296;
    }

    // // 0.65m
    else if (com_height_ == 0.65){
        K(0, 0) = 66.4281134896190;
        K(0, 1) = 2173.13307414793;
        K(0, 2) = 568.379709117415;
        K(0, 3) = 2.32202159217899;
        K(1, 0) = 2173.13307414793;
        K(1, 1) = 73309.6668377235;
        K(1, 2) = 19183.2680991541;
        K(1, 3) = 80.7128585695845;
        K(2, 0) = 568.379709117415;
        K(2, 1) = 19183.2680991541;
        K(2, 2) = 5019.91885890998;
        K(2, 3) = 21.1593587923669;
        K(3, 0) = 2.32202159217899;
        K(3, 1) = 80.7128585695845;
        K(3, 2) = 21.1593587923669;
        K(3, 3) = 0.0993494956670487;
    }


    // 0.68m
    else if (com_height_ == 0.68){
        K(0, 0) = 67.8563860285047;
        K(0, 1) = 2268.31636940822;
        K(0, 2) = 606.574465948491;
        K(0, 3) = 2.47294507780013;
        K(1, 0) = 2268.31636940822;
        K(1, 1) = 78097.9429536234;
        K(1, 2) = 20893.4283786941;
        K(1, 3) = 87.5477908301230;
        K(2, 0) = 606.574465948491;
        K(2, 1) = 20893.4283786941;
        K(2, 2) = 5589.73130292243;
        K(2, 3) = 23.4600661791894;
        K(3, 0) = 2.47294507780013;
        K(3, 1) = 87.5477908301230;
        K(3, 2) = 23.4600661791894;
        K(3, 3) = 0.108757490098902;
    }


    Eigen::MatrixXd Temp_mat;
    Eigen::MatrixXd Temp_mat_inv;
    Eigen::MatrixXd Ac_bar;
    Temp_mat.setZero(1, 1);
    Temp_mat_inv.setZero(1, 1);
    Ac_bar.setZero(4, 4);

    Temp_mat = R + B_bar_tran * K * B_bar;
    Temp_mat_inv = Temp_mat.inverse();

    Ac_bar = A_bar - B_bar * Temp_mat_inv * B_bar_tran * K * A_bar;

    Eigen::MatrixXd Ac_bar_tran(4, 4);
    Ac_bar_tran = Ac_bar.transpose();

    Gi.setZero(1, 1);
    Gx.setZero(1, 3);
    
    // 0.7m    
    if (com_height_ == 0.7){
        Gi(0, 0) = 562.367264784428; //Temp_mat_inv * B_bar_tran * K * I_bar ;
        //Gx = Temp_mat_inv * B_bar_tran * K * F_bar ;
        Gx(0, 0) = 38686.4538399671;
        Gx(0, 1) = 10800.0433241572;
        Gx(0, 2) = 128.255400226113;
    }

    // 0.728m
    else if (com_height_ == 0.728){
        Gi(0, 0) = 556.382091536873; //Temp_mat_inv * B_bar_tran * K * I_bar ;
        //Gx = Temp_mat_inv * B_bar_tran * K * F_bar ;
        Gx(0, 0) = 38991.9807941560;
        Gx(0, 1) = 11086.4028841705;
        Gx(0, 2) = 130.234568101964;
    }

    // 0.65m
    else if (com_height_ == 0.65){
        Gi(0, 0) = 573.462734668883; //Temp_mat_inv * B_bar_tran * K * I_bar ;
        //Gx = Temp_mat_inv * B_bar_tran * K * F_bar ;
        Gx(0, 0) = 38094.0476205746;
        Gx(0, 1) = 10274.4390649459;
        Gx(0, 2) = 124.583981245267;
    }

    // 0.68m
    else if (com_height_ == 0.68){
        Gi(0, 0) = 566.740708905623; //Temp_mat_inv * B_bar_tran * K * I_bar ;
        //Gx = Temp_mat_inv * B_bar_tran * K * F_bar ;
        Gx(0, 0) = 38456.9763214859;
        Gx(0, 1) = 10592.0336283141;
        Gx(0, 2) = 126.808547874660;
    }


    Eigen::MatrixXd X_bar;
    Eigen::Vector4d X_bar_col;
    X_bar.setZero(4, NL);
    X_bar_col.setZero();
    X_bar_col = -Ac_bar_tran * K * I_bar;

    for (int i = 0; i < NL; i++)
    {
        X_bar.block<4, 1>(0, i) = X_bar_col;
        X_bar_col = Ac_bar_tran * X_bar_col;
    }

    Gd.setZero(NL);
    Eigen::VectorXd Gd_col(1);
    Gd_col(0) = -Gi(0, 0);

    for (int i = 0; i < NL; i++)
    {
        Gd.segment(i, 1) = Gd_col;
        Gd_col = Temp_mat_inv * B_bar_tran * X_bar.col(i);
    }
}

void CustomController::previewcontroller(double dt, int NL, int tick, 
                                         Eigen::Vector3d &x_k, Eigen::Vector3d &y_k, double &UX, double &UY,
                                         const Eigen::MatrixXd &Gi, const Eigen::VectorXd &Gd, const Eigen::MatrixXd &Gx, 
                                         const Eigen::MatrixXd &A,  const Eigen::VectorXd &B,  const Eigen::MatrixXd &C)
{

    Eigen::VectorXd px, py;
    px.setZero(1); px = C * x_k;
    py.setZero(1); py = C * y_k;
    EX_preview_ -= (px(0) - ref_zmp_(tick,0)) * Gi(0, 0);
    EY_preview_ -= (py(0) - ref_zmp_(tick,1)) * Gi(0, 0);
    double sum_Gd_px_ref = 0, sum_Gd_py_ref = 0;
    for (int i = 0; i < NL; i++)
    {
        sum_Gd_px_ref = sum_Gd_px_ref - Gd(i) * (ref_zmp_(tick + 1 + i,0));
        sum_Gd_py_ref = sum_Gd_py_ref - Gd(i) * (ref_zmp_(tick + 1 + i,1));
    }
    Eigen::VectorXd GX_X; GX_X.setZero(1);
    GX_X = -Gx * (x_k);
    Eigen::VectorXd GX_Y; GX_Y.setZero(1);
    GX_Y = -Gx * (y_k);

    UX = EX_preview_ + sum_Gd_px_ref + GX_X(0);
    UY = EY_preview_ + sum_Gd_py_ref + GX_Y(0);    

    x_k = A * x_k + B * UX;
    y_k = A * y_k + B * UY;    
}


void CustomController::getFootTrajectory() 
{
    Eigen::Vector6d target_swing_foot; target_swing_foot.setZero();
    target_swing_foot = foot_step_support_frame_.row(0).transpose().segment(0,6);
    Eigen::Isometry3d &support_foot_traj           = (is_lfoot_support == true && is_rfoot_support == false) ? lfoot_trajectory_support_ : rfoot_trajectory_support_;
    Eigen::Vector3d &support_foot_traj_euler       = (is_lfoot_support == true && is_rfoot_support == false) ? lfoot_trajectory_euler_support_ : rfoot_trajectory_euler_support_;

    Eigen::Isometry3d &swing_foot_traj             = (is_lfoot_support == true && is_rfoot_support == false) ? rfoot_trajectory_support_ : lfoot_trajectory_support_;
    Eigen::Vector3d &swing_foot_traj_euler         = (is_lfoot_support == true && is_rfoot_support == false) ? rfoot_trajectory_euler_support_ : lfoot_trajectory_euler_support_;
    const Eigen::Isometry3d &swing_foot_init       = (is_lfoot_support == true && is_rfoot_support == false) ? rfoot_support_init_yaw_ : lfoot_support_init_yaw_;
    const Eigen::Vector3d &swing_foot_euler_init   = (is_lfoot_support == true && is_rfoot_support == false) ? rfoot_support_euler_init_yaw_ : lfoot_support_euler_init_yaw_;

    if (is_dsp1 == true)
    {
        support_foot_traj.translation().setZero();
        support_foot_traj_euler.setZero();

        swing_foot_traj.translation() = swing_foot_init.translation();
        // swing_foot_traj.translation()(2) = 0.0;
        swing_foot_traj_euler = swing_foot_euler_init;

        target_swing_state_stance_frame_.segment(7, 6) << 0., 0, 0, 0, 0, 0;
    }
    else if (is_ssp == true)
    {
        support_foot_traj.translation().setZero();
        support_foot_traj_euler.setZero();

        if (walking_tick < t_dsp_(0) + t_ssp_(0) / 2.0)
        {
            swing_foot_traj.translation()(2) = DyrosMath::cubic(walking_tick, 
                                                                t_dsp_(0),
                                                                t_dsp_(0) + t_ssp_(0) / 2.0, 
                                                                0.0, foot_height_(0), 
                                                                0.0, 0.0);
            target_swing_state_stance_frame_(9) = DyrosMath::cubicDot(walking_tick, 
                                                                t_dsp_(0),
                                                                t_dsp_(0) + t_ssp_(0) / 2.0, 
                                                                0.0, foot_height_(0), 
                                                                0.0, 0.0);

        }
        else
        {
            swing_foot_traj.translation()(2) = DyrosMath::cubic(walking_tick, 
                                                                t_dsp_(0) + t_ssp_(0) / 2.0,
                                                                t_dsp_(0) + t_ssp_(0), 
                                                                foot_height_(0), target_swing_foot(2), 
                                                                0.0, 0.0);
            target_swing_state_stance_frame_(9) = DyrosMath::cubicDot(walking_tick, 
                                                                t_dsp_(0) + t_ssp_(0) / 2.0,
                                                                t_dsp_(0) + t_ssp_(0), 
                                                                foot_height_(0), target_swing_foot(2), 
                                                                0.0, 0.0);                                                                
        }

        swing_foot_traj.translation().segment(0,2) = DyrosMath::cubicVector<2>(walking_tick, 
                                                                               t_dsp_(0),
                                                                               t_dsp_(0) + t_ssp_(0), 
                                                                               swing_foot_init.translation().segment(0,2), target_swing_foot.segment(0,2),
                                                                               Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero());
        target_swing_state_stance_frame_(7) = DyrosMath::cubicDot(walking_tick, 
                                                                               t_dsp_(0),
                                                                               t_dsp_(0) + t_ssp_(0), 
                                                                               swing_foot_init.translation()(0), target_swing_foot(0),
                                                                               0., 0.);
        target_swing_state_stance_frame_(8) = DyrosMath::cubicDot(walking_tick, 
                                                                               t_dsp_(0),
                                                                               t_dsp_(0) + t_ssp_(0), 
                                                                               swing_foot_init.translation()(1), target_swing_foot(1),
                                                                               0., 0.);

        swing_foot_traj_euler.setZero();
        swing_foot_traj_euler(2) = DyrosMath::cubic(walking_tick, 
                                                    t_dsp_(0), 
                                                    t_dsp_(0) + t_ssp_(0), 
                                                    swing_foot_euler_init(2), target_swing_foot(5), 
                                                    0.0, 0.0);
        target_swing_state_stance_frame_(12) =  DyrosMath::cubicDot(walking_tick, 
                                                    t_dsp_(0), 
                                                    t_dsp_(0) + t_ssp_(0), 
                                                    swing_foot_euler_init(2), target_swing_foot(5), 
                                                    0.0, 0.0);
    }
    else if (is_dsp2 == true)
    {
        support_foot_traj_euler.setZero();
        
        swing_foot_traj.translation() = target_swing_foot.segment(0,3);
        swing_foot_traj_euler = target_swing_foot.segment(3,3);
        target_swing_state_stance_frame_.segment(7, 6) << 0., 0, 0, 0, 0, 0;
    }

    swing_foot_traj.linear() = DyrosMath::Euler2rot(0., 0., swing_foot_traj_euler(2));
    support_foot_traj.linear() = DyrosMath::Euler2rot(support_foot_traj_euler(0), support_foot_traj_euler(1), support_foot_traj_euler(2));
    target_swing_state_stance_frame_.segment(0,3) = swing_foot_traj.translation();
    Eigen::Quaterniond swing_quat(swing_foot_traj.linear());
    target_swing_state_stance_frame_(3) = swing_quat.x();
    target_swing_state_stance_frame_(4) = swing_quat.y();
    target_swing_state_stance_frame_(5) = swing_quat.z();
    target_swing_state_stance_frame_(6) = swing_quat.w();

}

void CustomController::computeIkControl(const Eigen::Isometry3d &float_trunk_transform, const Eigen::Isometry3d &float_lleg_transform, const Eigen::Isometry3d &float_rleg_transform, Eigen::Vector12d &q_des)
{
    Eigen::Vector3d R_r, R_D, L_r, L_D;

    L_D << 0.11, +0.1025, -0.1025;
    R_D << 0.11, -0.1025, -0.1025;

    L_r = float_lleg_transform.rotation().transpose() * (float_trunk_transform.translation() + float_trunk_transform.rotation() * L_D - float_lleg_transform.translation());
    R_r = float_rleg_transform.rotation().transpose() * (float_trunk_transform.translation() + float_trunk_transform.rotation() * R_D - float_rleg_transform.translation());

    double R_C = 0, L_C = 0, L_upper = 0.351, L_lower = 0.351, R_alpha = 0, L_alpha = 0;

    L_C = sqrt(pow(L_r(0), 2) + pow(L_r(1), 2) + pow(L_r(2), 2));
    R_C = sqrt(pow(R_r(0), 2) + pow(R_r(1), 2) + pow(R_r(2), 2));
     
    double knee_acos_var_L = 0;
    double knee_acos_var_R = 0;

    knee_acos_var_L = (pow(L_upper, 2) + pow(L_lower, 2) - pow(L_C, 2))/ (2 * L_upper * L_lower);
    knee_acos_var_R = (pow(L_upper, 2) + pow(L_lower, 2) - pow(R_C, 2))/ (2 * L_upper * L_lower);

    knee_acos_var_L = DyrosMath::minmax_cut(knee_acos_var_L, -0.99, + 0.99);
    knee_acos_var_R = DyrosMath::minmax_cut(knee_acos_var_R, -0.99, + 0.99);

    q_des(3) = (-acos(knee_acos_var_L) + M_PI);  
    q_des(9) = (-acos(knee_acos_var_R) + M_PI);
    q_des(5) = atan2(L_r(1), L_r(2));                                                                                  // Ankle roll
    q_des(11) = atan2(R_r(1), R_r(2));


    // L_alpha = asin(DyrosMath::minmax_cut(L_upper / L_C * sin(M_PI - q_des(3)), -0.99, 0.99) );

    // R_alpha = asin(DyrosMath::minmax_cut(L_upper / R_C * sin(M_PI - q_des(9)), -0.99, 0.99));

    L_alpha =asin( L_upper / L_C * sin(M_PI - q_des(3)));

    R_alpha = asin(L_upper / R_C * sin(M_PI - q_des(9)));
    
    q_des(4) = -atan2(L_r(0), sqrt(pow(L_r(1), 2) + pow(L_r(2), 2))) - L_alpha;
    q_des(10) = -atan2(R_r(0), sqrt(pow(R_r(1), 2) + pow(R_r(2), 2))) - R_alpha;

    Eigen::Matrix3d R_Knee_Ankle_Y_rot_mat, L_Knee_Ankle_Y_rot_mat;
    Eigen::Matrix3d R_Ankle_X_rot_mat, L_Ankle_X_rot_mat;
    Eigen::Matrix3d R_Hip_rot_mat, L_Hip_rot_mat;

    L_Knee_Ankle_Y_rot_mat = DyrosMath::rotateWithY(-q_des(3) - q_des(4));
    L_Ankle_X_rot_mat = DyrosMath::rotateWithX(-q_des(5));
    R_Knee_Ankle_Y_rot_mat = DyrosMath::rotateWithY(-q_des(9) - q_des(10));
    R_Ankle_X_rot_mat = DyrosMath::rotateWithX(-q_des(11));

    L_Hip_rot_mat.setZero();
    R_Hip_rot_mat.setZero();

    L_Hip_rot_mat = float_trunk_transform.rotation().transpose() * float_lleg_transform.rotation() * L_Ankle_X_rot_mat * L_Knee_Ankle_Y_rot_mat;
    R_Hip_rot_mat = float_trunk_transform.rotation().transpose() * float_rleg_transform.rotation() * R_Ankle_X_rot_mat * R_Knee_Ankle_Y_rot_mat;

    q_des(0) = atan2(-L_Hip_rot_mat(0, 1), L_Hip_rot_mat(1, 1));                                                       // Hip yaw
    q_des(1) = atan2(L_Hip_rot_mat(2, 1), -L_Hip_rot_mat(0, 1) * sin(q_des(0)) + L_Hip_rot_mat(1, 1) * cos(q_des(0))); // Hip roll
    q_des(2) = atan2(-L_Hip_rot_mat(2, 0), L_Hip_rot_mat(2, 2));                                                       // Hip pitch
    q_des(3) = q_des(3);                                                                                               // Knee pitch
    q_des(4) = q_des(4);                                                                                               // Ankle pitch

    q_des(6) = atan2(-R_Hip_rot_mat(0, 1), R_Hip_rot_mat(1, 1));
    q_des(7) = atan2(R_Hip_rot_mat(2, 1), -R_Hip_rot_mat(0, 1) * sin(q_des(6)) + R_Hip_rot_mat(1, 1) * cos(q_des(6)));
    q_des(8) = atan2(-R_Hip_rot_mat(2, 0), R_Hip_rot_mat(2, 2));
    q_des(9) = q_des(9);
    q_des(10) = q_des(10);
}



void CustomController::getTargetState(){

    target_com_state_stance_frame_.segment(0, 3) = com_desired_;
    Eigen::Quaterniond com_quat(DyrosMath::rotateWithZ(ref_com_yaw_(walking_tick+1)));
    target_com_state_stance_frame_(3) = com_quat.x();
    target_com_state_stance_frame_(4) = com_quat.y();
    target_com_state_stance_frame_(5) = com_quat.z();
    target_com_state_stance_frame_(6) = com_quat.w();
    target_com_state_stance_frame_.segment(7, 3) = com_desired_dot_;
    target_com_state_stance_frame_(10) = 0.;
    target_com_state_stance_frame_(11) = 0.;
    target_com_state_stance_frame_(12) = ref_com_yawvel_(walking_tick+1);

    if (com_height_ == 0.68){
        target_com_state_float_frame_.translation() << (rd_cc_.link_[Pelvis].rotm.transpose() * (rd_cc_.link_[Pelvis].xpos-rd_cc_.link_[COM_id].xpos))(0),
        (rd_cc_.link_[Pelvis].rotm.transpose() * (rd_cc_.link_[Pelvis].xpos-rd_cc_.link_[COM_id].xpos))(1), 
        (rd_cc_.link_[Pelvis].rotm.transpose() * (rd_cc_.link_[Pelvis].xpos-rd_cc_.link_[COM_id].xpos))(2);

    }
    else if (com_height_ == 0.728){
        target_com_state_float_frame_.translation() << DyrosMath::minmax_cut((rd_cc_.link_[Pelvis].rotm.transpose() * (rd_cc_.link_[Pelvis].xpos-rd_cc_.link_[COM_id].xpos))(0), -0.15, 0.),
        DyrosMath::minmax_cut((rd_cc_.link_[Pelvis].rotm.transpose() * (rd_cc_.link_[Pelvis].xpos-rd_cc_.link_[COM_id].xpos))(1), -0.04, 0.04), 
        DyrosMath::minmax_cut((rd_cc_.link_[Pelvis].rotm.transpose() * (rd_cc_.link_[Pelvis].xpos-rd_cc_.link_[COM_id].xpos))(2), 0., 0.04);
    }
    else{
        std::cout << "WRONG COM HEIGHT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
        exit(0);
    }
    target_com_state_float_frame_.linear() = Eigen::Matrix3d::Identity();

    Eigen::Vector4d swingq = target_swing_state_stance_frame_.segment(3,4);
    Eigen::Vector4d comq = target_com_state_stance_frame_.segment(3,4);
    Eigen::Quaterniond target_swing_state_stance_frame_quat_(swingq(3), swingq(0), swingq(1), swingq(2));
    Eigen::Quaterniond target_com_state_stance_frame_quat_(comq(3), comq(0), comq(1), comq(2));

    target_lfoot_state_float_frame_.translation() = phase_indicator_(0)*target_com_state_stance_frame_quat_.toRotationMatrix().transpose() * (target_swing_state_stance_frame_.segment(0,3) - target_com_state_stance_frame_.segment(0,3)) + (1-phase_indicator_(0)) * target_com_state_stance_frame_quat_.toRotationMatrix().transpose()*(-target_com_state_stance_frame_.segment(0,3));
    target_lfoot_state_float_frame_.linear() = phase_indicator_(0)*target_com_state_stance_frame_quat_.toRotationMatrix().transpose()*target_swing_state_stance_frame_quat_.toRotationMatrix() + (1-phase_indicator_(0))*target_com_state_stance_frame_quat_.toRotationMatrix().transpose();
    target_rfoot_state_float_frame_.translation() = (1-phase_indicator_(0))*target_com_state_stance_frame_quat_.toRotationMatrix().transpose() * (target_swing_state_stance_frame_.segment(0,3) - target_com_state_stance_frame_.segment(0,3)) + phase_indicator_(0) * target_com_state_stance_frame_quat_.toRotationMatrix().transpose()*(-target_com_state_stance_frame_.segment(0,3));
    target_rfoot_state_float_frame_.linear() = (1-phase_indicator_(0))*target_com_state_stance_frame_quat_.toRotationMatrix().transpose()*target_swing_state_stance_frame_quat_.toRotationMatrix() + phase_indicator_(0)*target_com_state_stance_frame_quat_.toRotationMatrix().transpose();

    computeIkControl(target_com_state_float_frame_, target_lfoot_state_float_frame_, target_rfoot_state_float_frame_, q_leg_desired_);
}


void CustomController::windupPreview(){
    EX_preview_ = 0.;
    EY_preview_ = 0.;
}

void CustomController::updateNextStepTime()
{       
    // std::cout << "walking time : " << (walking_tick) / hz_ <<  ", Value : " << value_ << std::endl;
    walking_tick++;
}

void CustomController::walkingStateMachine()
{
    if (foot_step_(0, 6) == 1) 
    {
        is_lfoot_support = true;
        is_rfoot_support = false;
    }
    else if (foot_step_(0, 6) == 0) 
    {
        is_lfoot_support = false;
        is_rfoot_support = true;
    }

    if (walking_tick < t_dsp_(0))
    {
        is_dsp1 = true;
        is_ssp  = false;
        is_dsp2 = false;
    }
    else if (walking_tick >= t_dsp_(0) && walking_tick < t_total_(0) - t_dsp_(0))
    {
        is_dsp1 = false;
        is_ssp  = true;
        is_dsp2 = false;
    }
    else
    {
        is_dsp1 = false;
        is_ssp  = false;
        is_dsp2 = true;
    } 
}