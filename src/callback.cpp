#include "callback.hpp"
#include "initial.hpp"
#include <ros/ros.h>
#include <tf2/LinearMath/Quaternion.h>
#include <tf2/LinearMath/Matrix3x3.h>

// グローバル変数の定義
geometry_msgs::PoseStamped body_link_pose;
bool body_link_pose_received = false;

geometry_msgs::PoseStamped v1_rear_left_steering_pose;
bool v1_rear_left_steering_pose_received = false;

geometry_msgs::PoseStamped v1_front_left_steering_pose;
bool v1_front_left_steering_pose_received = false;


void jointStateCallback(const sensor_msgs::JointState::ConstPtr& msg)
{
    for (size_t i = 0; i < msg->name.size(); ++i) {
        const auto& joint = msg->name[i];
        double pos = msg->position[i];  // [m] for prismatic
        double vel = msg->velocity[i];
        double wheel_angle_FL,wheel_angle_FR,wheel_angle_RL,wheel_angle_RR
        double wheel_angle_vel_FL, wheel_angle_vel_FR, wheel_angle_vel_RL, wheel_angle_vel_RR
        // // — プリズマティックジョイントの変位を delta_pos に格納 —
        // if (joint == "linkage_point_to_v1") {
        //     delta_pos[0] = pos;
        //     linkage_point1_pose_received = true;
        //     // （PoseStamped の更新は必要に応じて）
        // }
        // else if (joint == "linkage_point_to_v2") {
        //     delta_pos[1] = pos;
        //     linkage_point2_pose_received = true;
        // }
        // else if (joint == "linkage_point_to_v3") {
        //     delta_pos[2] = pos;
        //     linkage_point3_pose_received = true;
        // }

        // …既存の revolute ジョイント処理…
        if (msg->name[i] == "front_right_steering") {
            phidot = vel;
            qdot_twist[3] = phidot
            x_d[4] = phidot;
        }
        if (joint == "front_left_wheel") {
            wheel_angle_FL     = pos;
            wheel_angle_vel_FL  = vel;
        }
        else if (joint == "front_right_wheel") {
            wheel_angle_FR     = pos;
            wheel_angle_vel_FR   = vel;
        }
        else if (joint == "rear_left_wheel") {
            wheel_angle_RL     = pos;
            wheel_angle_ve_RL   = vel;
        }
        else if (joint == "rear_right_wheel") {
            wheel_angle_RR     = pos;
            wheel_angle_vel_RR   = vel;
        }
        q_twist[4] = (wheel_angle_FL + wheel_angle_FR)/2
        qdot_twist[4] = (wheel_angle_vel_FL + wheel_angle_vel_FR)/2
        q_twist[5] = (wheel_angle_RL + wheel_angle_RR)/2
        qdot_twist[5] = (wheel_angle_ve_RL + wheel_angle_vel_RR)/2

    }
}

// body_link 用のコールバック
void trueBodyLinkCallback(const nav_msgs::Odometry::ConstPtr& msg)
{
    // — 既存の姿勢・オイラー角取得処理 —
    body_link_pose.header = msg->header;
    body_link_pose.pose   = msg->pose.pose;
    body_link_pose_received = true;

    tf2::Quaternion q(
      body_link_pose.pose.orientation.x,
      body_link_pose.pose.orientation.y,
      body_link_pose.pose.orientation.z,
      body_link_pose.pose.orientation.w);
    double roll, pitch, yaw;
    tf2::Matrix3x3(q).getRPY(roll, pitch, yaw);

    true_body_pos[0] = body_link_pose.pose.position.x;
    true_body_pos[1] = body_link_pose.pose.position.y;
    true_body_yaw    = yaw;
    x_old[1] = body_link_pose.pose.position.x - cos(yaw) * (lv/2);
    x_old[2] = body_link_pose.pose.position.y - sin(yaw) * (lv/2);
    x_old[3] = yaw;
    q_twist[0] = body_link_pose.pose.position.x - cos(yaw) * (lv/2);
    q_twist[1] = body_link_pose.pose.position.y - sin(yaw) * (lv/2);
    q_twist[2] = yaw;  
    got_body_pose = true;

    // — u1_act の取得ロジック追加 —
    // body frame 前方速度 vx, 横方向速度 vy を取得
    qdot_twist[0] = msg->twist.twist.linear.x;
    qdot_twist[1] = msg->twist.twist.linear.y;
    qdot_twist[2] = msg->twist.twist.angular.z;
    x_d[1] = qdot_twist[0];
    x_d[2] = qdot_twist[1];
    x_d[3] = qdot_twist[2];
    // twist が body_frame ならそのまま vx を前進速度とする
    // もし world_frame であれば下記のように投影
    // u1_act = vx * cos(yaw) + vy * sin(yaw);
}

void trueV1FrontLeftSteeringCallback(const nav_msgs::Odometry::ConstPtr& msg)
{
    v1_front_left_steering_pose.header = msg->header;
    v1_front_left_steering_pose.pose = msg->pose.pose;
    v1_front_left_steering_pose_received = true;

    tf2::Quaternion q;
    q.setX(v1_front_left_steering_pose.pose.orientation.x);
    q.setY(v1_front_left_steering_pose.pose.orientation.y);
    q.setZ(v1_front_left_steering_pose.pose.orientation.z);
    q.setW(v1_front_left_steering_pose.pose.orientation.w);

    double roll, pitch, yaw;
    tf2::Matrix3x3(q).getRPY(roll, pitch, yaw);

    true_steering = yaw;

    q_twist[3] = true_steering - true_body_yaw;
    x_old[4] = q_twist[3];
    got_steering_pose = true;

    
}
