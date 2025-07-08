#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <ctime>
#include <sstream>
#include <iostream>
#include "initial.hpp"
#include "mathFunc.h"		// 数学関数のヘッダーファイル
#include "Bezier.h"			// Bezier 曲線の関数
#include "vehicle.hpp"       // Vehicle クラスの宣言
#include "callback.hpp"     // コールバック関数の宣言
#include "initial.hpp"		// 初期値設定のヘッダーファイル
#include "getInputValue_dynamics.hpp"
#include "dynamics_integrator.hpp"
#include "pidTau.hpp"
#include "csvLogger.hpp"
#include <ros/package.h>

using namespace std;


Search searchP(std::vector<double>& x);
Search searchPP(std::vector<double>& x);



int main(int argc, char** argv) {
    ros::init(argc, argv, "steering_desired_controller");
    ros::NodeHandle nh;
	ros::Rate loop_rate(100);
	// jointをsubscribe
    ros::Subscriber joint_state_sub    = nh.subscribe("/vehicle_4ws/joint_states", 10, jointStateCallback);
    ros::Subscriber true_bodylink_sub = nh.subscribe("/vehicle_4ws/true_body_link", 10, trueBodyLinkCallback);
	ros::Subscriber front_left_steering_sub = nh.subscribe("/vehicle_4ws/true_front_left_steering_link", 10, trueV1FrontLeftSteeringCallback);

	//データファイル作成
  	std::string pkg = ros::package::getPath("4ws_backstepping");
  	std::string data_dir = pkg + "/data";

  	CSVLogger logger(data_dir, 100000);



	PIDGains driveGains{400.0, 40.0, 0.01};  // (Kp, Ki, Kd)
    PIDGains steerGains{2.0, 0.2, 0.005};
    // Vehicle クラスのインスタンス生成
    Vehicle vehicle1(nh, "v1");
	//制御入力のクラスをインスタンス化
	getInputValue getInputValue(0.01);
	//動力学計算のクラスをインスタンス化
	DynamicsIntegrator integrator(
        M_mass, I_theta, lv, GRAV, rho,
        driveGains, steerGains,
        0.01,
		getInputValue
    );

	//gazebo上の初期値をx_oldに代入
	while (ros::ok() && !(got_body_pose && got_steering_pose)) {
		ros::spinOnce();
		loop_rate.sleep();
	}
	//デバッグ用ログ出力
	ROS_INFO("Locking initial pose and calling initial(): xdot=%.3f, ydot=%.3f, thetadot=%.3f ,phiRdot=%.3f, varphiRdot=%.3f, phiFdot=%.3f, varphiFdot=%.3f",
           x_d[1], x_d[2], x_d[3], x_d[4], x_d[5], x_d[6], x_d[7]);
		   
	//初期値を設定
	initial(t_max, h);


	//曲率計算
	qs[0] = 0.0;
	double dt = 1.0 / (Q - 1);
	for (int i = 1; i < Q; ++i) {
    	qs[i] = i * dt;
	}

	for (int i = 0; i < Q; i++) {

		R[i][0] = Rx(Bx, qs, i);
		R[i][1] = Ry(By, qs, i);
		dRdq[i][0] = d1Rxdq1(Bx, qs, i);
		dRdq[i][1] = d1Rydq1(By, qs, i);

		cs[i][0] = (-(d1Rydq1(By, qs, i) * d2Rxdq2(Bx, qs, i)) + d1Rxdq1(Bx, qs, i) * d2Rydq2(By, qs, i)) /
			Power(Power(d1Rxdq1(Bx, qs, i), 2) + Power(d1Rydq1(By, qs, i), 2), 1.5);


		cs[i][1] = (Power(d1Rydq1(By, qs, i), 2) * (3 * d2Rxdq2(Bx, qs, i) * d2Rydq2(By, qs, i) - d1Rydq1(By, qs, i) * d3Rxdq3(Bx, qs, i))
			- Power(d1Rxdq1(Bx, qs, i), 2) * (3 * d2Rxdq2(Bx, qs, i) * d2Rydq2(By, qs, i) + d1Rydq1(By, qs, i) * d3Rxdq3(Bx, qs, i))
			+ Power(d1Rxdq1(Bx, qs, i), 3) * d3Rydq3(By, qs, i) + d1Rxdq1(Bx, qs, i) * d1Rydq1(By, qs, i) * (3 * Power(d2Rxdq2(Bx, qs, i), 2)
				- 3 * Power(d2Rydq2(By, qs, i), 2) + d1Rydq1(By, qs, i) * d3Rydq3(By, qs, i))) / Power(Power(d1Rxdq1(Bx, qs, i), 2) + Power(d1Rydq1(By, qs, i), 2), 3);


		cs[i][2] = (-(Power(d1Rxdq1(Bx, qs, i), 4) * (4 * d2Rydq2(By, qs, i) * d3Rxdq3(Bx, qs, i)
			+ 6 * d2Rxdq2(Bx, qs, i) * d3Rydq3(By, qs, i) + d1Rydq1(By, qs, i) * d4Rxdq4(Bx, qs, i)))
			+ Power(d1Rxdq1(Bx, qs, i), 2) * d1Rydq1(By, qs, i) * (-15 * Power(d2Rxdq2(Bx, qs, i), 3)
				+ d2Rxdq2(Bx, qs, i) * (39 * Power(d2Rydq2(By, qs, i), 2) - 2 * d1Rydq1(By, qs, i) * d3Rydq3(By, qs, i))
				+ 2 * d1Rydq1(By, qs, i) * (d2Rydq2(By, qs, i) * d3Rxdq3(Bx, qs, i) - d1Rydq1(By, qs, i) * d4Rxdq4(Bx, qs, i)))
			+ Power(d1Rydq1(By, qs, i), 3) * (3 * Power(d2Rxdq2(Bx, qs, i), 3) + d2Rxdq2(Bx, qs, i) * (-15 * Power(d2Rydq2(By, qs, i), 2)
				+ 4 * d1Rydq1(By, qs, i) * d3Rydq3(By, qs, i)) + d1Rydq1(By, qs, i) * (6 * d2Rydq2(By, qs, i) * d3Rxdq3(Bx, qs, i)
					- d1Rydq1(By, qs, i) * d4Rxdq4(Bx, qs, i))) + Power(d1Rxdq1(Bx, qs, i), 5) * d4Rydq4(By, qs, i)
			+ d1Rxdq1(Bx, qs, i) * Power(d1Rydq1(By, qs, i), 2) * (-39 * Power(d2Rxdq2(Bx, qs, i), 2) * d2Rydq2(By, qs, i)
				+ 15 * Power(d2Rydq2(By, qs, i), 3) + 10 * d1Rydq1(By, qs, i) * d2Rxdq2(Bx, qs, i) * d3Rxdq3(Bx, qs, i)
				- 10 * d1Rydq1(By, qs, i) * d2Rydq2(By, qs, i) * d3Rydq3(By, qs, i) + Power(d1Rydq1(By, qs, i), 2) * d4Rydq4(By, qs, i))
			+ Power(d1Rxdq1(Bx, qs, i), 3) * (15 * Power(d2Rxdq2(Bx, qs, i), 2) * d2Rydq2(By, qs, i) - 3 * Power(d2Rydq2(By, qs, i), 3)
				+ 10 * d1Rydq1(By, qs, i) * d2Rxdq2(Bx, qs, i) * d3Rxdq3(Bx, qs, i) - 10 * d1Rydq1(By, qs, i) * d2Rydq2(By, qs, i) * d3Rydq3(By, qs, i)
				+ 2 * Power(d1Rydq1(By, qs, i), 2) * d4Rydq4(By, qs, i))) / Power(Power(d1Rxdq1(Bx, qs, i), 2) + Power(d1Rydq1(By, qs, i), 2), 4.5);
		
		cs[i][3] = (3*Power(d1Rydq1(By, qs, i)*d2Rxdq2(Bx, qs, i) - d1Rxdq1(Bx, qs, i)*d2Rydq2(By, qs, i),2)*(Power(d1Rydq1(By, qs, i),2)*(3*d2Rxdq2(Bx, qs, i)*d2Rydq2(By, qs, i) - d1Rydq1(By, qs, i)*d3Rxdq3(Bx, qs, i)) - Power(d1Rxdq1(Bx, qs, i),2)*(3*d2Rxdq2(Bx, qs, i)*d2Rydq2(By, qs, i) +    d1Rydq1(By, qs, i)*d3Rxdq3(Bx, qs, i)) + Power(d1Rxdq1(Bx, qs, i),3)*d3Rydq3(By, qs, i) + d1Rxdq1(Bx, qs, i)*d1Rydq1(By, qs, i)* (3*Power(d2Rxdq2(Bx, qs, i),2) - 3*Power(d2Rydq2(By, qs, i),2) +    d1Rydq1(By, qs, i)*d3Rydq3(By, qs, i))))/  Power(Power(d1Rxdq1(Bx, qs, i),2) + Power(d1Rydq1(By, qs, i),2),6) +  (3*Power(d1Rydq1(By, qs, i)*d2Rxdq2(Bx, qs, i) - d1Rxdq1(Bx, qs, i)*d2Rydq2(By, qs, i),2)*    (-(Power(d1Rxdq1(Bx, qs, i),4)*(4*d2Rydq2(By, qs, i)*d3Rxdq3(Bx, qs, i) +            6*d2Rxdq2(Bx, qs, i)*d3Rydq3(By, qs, i) + d1Rydq1(By, qs, i)*d4Rxdq4(Bx, qs, i))) +       Power(d1Rxdq1(Bx, qs, i),2)*d1Rydq1(By, qs, i)*       (-15*Power(d2Rxdq2(Bx, qs, i),3) +          d2Rxdq2(Bx, qs, i)*(39*Power(d2Rydq2(By, qs, i),2) -             2*d1Rydq1(By, qs, i)*d3Rydq3(By, qs, i)) +          2*d1Rydq1(By, qs, i)*(d2Rydq2(By, qs, i)*d3Rxdq3(Bx, qs, i) -             d1Rydq1(By, qs, i)*d4Rxdq4(Bx, qs, i))) +       Power(d1Rydq1(By, qs, i),3)*(3*Power(d2Rxdq2(Bx, qs, i),3) +          d2Rxdq2(Bx, qs, i)*(-15*Power(d2Rydq2(By, qs, i),2) +             4*d1Rydq1(By, qs, i)*d3Rydq3(By, qs, i)) +          d1Rydq1(By, qs, i)*(6*d2Rydq2(By, qs, i)*d3Rxdq3(Bx, qs, i) -             d1Rydq1(By, qs, i)*d4Rxdq4(Bx, qs, i))) +       Power(d1Rxdq1(Bx, qs, i),5)*d4Rydq4(By, qs, i) + d1Rxdq1(Bx, qs, i)*Power(d1Rydq1(By, qs, i),2)*(-39*Power(d2Rxdq2(Bx, qs, i),2)*d2Rydq2(By, qs, i) + 15*Power(d2Rydq2(By, qs, i),3) + 10*d1Rydq1(By, qs, i)*d2Rxdq2(Bx, qs, i)*d3Rxdq3(Bx, qs, i) - 10*d1Rydq1(By, qs, i)*d2Rydq2(By, qs, i)*d3Rydq3(By, qs, i) + Power(d1Rydq1(By, qs, i),2)*d4Rydq4(By, qs, i)) + Power(d1Rxdq1(Bx, qs, i),3)*(15*Power(d2Rxdq2(Bx, qs, i),2)*d2Rydq2(By, qs, i) - 3*Power(d2Rydq2(By, qs, i),3) + 10*d1Rydq1(By, qs, i)*d2Rxdq2(Bx, qs, i)*d3Rxdq3(Bx, qs, i) - 10*d1Rydq1(By, qs, i)*d2Rydq2(By, qs, i)*d3Rydq3(By, qs, i) + 2*Power(d1Rydq1(By, qs, i),2)*d4Rydq4(By, qs, i))))/Power(Power(d1Rxdq1(Bx, qs, i),2) + Power(d1Rydq1(By, qs, i),2),7.5) + (d1Rxdq1(Bx, qs, i)*(Power(d1Rydq1(By, qs, i),5)*(13*Power(d2Rxdq2(Bx, qs, i),4) - 3*Power(d1Rydq1(By, qs, i),2)*Power(d3Rxdq3(Bx, qs, i),2) + Power(d2Rxdq2(Bx, qs, i),2)*(-87*Power(d2Rydq2(By, qs, i),2) + 16*d1Rydq1(By, qs, i)*d3Rydq3(By, qs, i)) + d2Rxdq2(Bx, qs, i)*(42*d1Rydq1(By, qs, i)*d2Rydq2(By, qs, i)*d3Rxdq3(Bx, qs, i) - 4*Power(d1Rydq1(By, qs, i),2)*d4Rxdq4(Bx, qs, i))) - 
  					    Power(d1Rxdq1(Bx, qs, i),7)*(10*d3Rxdq3(Bx, qs, i)*d3Rydq3(By, qs, i) + 
  					    5*d2Rydq2(By, qs, i)*d4Rxdq4(Bx, qs, i) + 10*d2Rxdq2(Bx, qs, i)*d4Rydq4(By, qs, i) + 
  					    d1Rydq1(By, qs, i)*d5Rxdq5(Bx, qs, i)) + 
  					      Power(d1Rxdq1(Bx, qs, i),3)*Power(d1Rydq1(By, qs, i),2)*
  					       (714*Power(d2Rxdq2(Bx, qs, i),3)*d2Rydq2(By, qs, i) - 
  					         22*d1Rydq1(By, qs, i)*Power(d2Rxdq2(Bx, qs, i),2)*d3Rxdq3(Bx, qs, i) + 
  					         d2Rxdq2(Bx, qs, i)*(-766*Power(d2Rydq2(By, qs, i),3) + 108*d1Rydq1(By, qs, i)*d2Rydq2(By, qs, i)*d3Rydq3(By, qs, i) + 8*Power(d1Rydq1(By, qs, i),2)*d4Rydq4(By, qs, i)) + d1Rydq1(By, qs, i)*(14*Power(d2Rydq2(By, qs, i),2)*d3Rxdq3(Bx, qs, i) + 23*d1Rydq1(By, qs, i)*d2Rydq2(By, qs, i)*d4Rxdq4(Bx, qs, i) + d1Rydq1(By, qs, i)*(22*d3Rxdq3(Bx, qs, i)*d3Rydq3(By, qs, i) - 3*d1Rydq1(By, qs, i)*d5Rxdq5(Bx, qs, i)))) + d1Rxdq1(Bx, qs, i)*Power(d1Rydq1(By, qs, i),4)*(-301*Power(d2Rxdq2(Bx, qs, i),3)*d2Rydq2(By, qs, i) + 83*d1Rydq1(By, qs, i)*Power(d2Rxdq2(Bx, qs, i),2)*d3Rxdq3(Bx, qs, i) + d2Rxdq2(Bx, qs, i)*(279*Power(d2Rydq2(By, qs, i),3) - 134*d1Rydq1(By, qs, i)*d2Rydq2(By, qs, i)*d3Rydq3(By, qs, i) + 9*Power(d1Rydq1(By, qs, i),2)*d4Rydq4(By, qs, i)) + d1Rydq1(By, qs, i)*(-87*Power(d2Rydq2(By, qs, i),2)*d3Rxdq3(Bx, qs, i) + 14*d1Rydq1(By, qs, i)*d2Rydq2(By, qs, i)*d4Rxdq4(Bx, qs, i) + d1Rydq1(By, qs, i)*(16*d3Rxdq3(Bx, qs, i)*d3Rydq3(By, qs, i) - d1Rydq1(By, qs, i)*d5Rxdq5(Bx, qs, i)))) + Power(d1Rxdq1(Bx, qs, i),5)*(-105*Power(d2Rxdq2(Bx, qs, i),3)*d2Rydq2(By, qs, i) - 105*d1Rydq1(By, qs, i)*Power(d2Rxdq2(Bx, qs, i),2)*d3Rxdq3(Bx, qs, i) + d2Rxdq2(Bx, qs, i)*(75*Power(d2Rydq2(By, qs, i),3) + 242*d1Rydq1(By, qs, i)*d2Rydq2(By, qs, i)*d3Rydq3(By, qs, i) - 11*Power(d1Rydq1(By, qs, i),2)*d4Rydq4(By, qs, i)) + 
  					         d1Rydq1(By, qs, i)*(101*Power(d2Rydq2(By, qs, i),2)*d3Rxdq3(Bx, qs, i) + 
  					            4*d1Rydq1(By, qs, i)*d2Rydq2(By, qs, i)*d4Rxdq4(Bx, qs, i) - 
  					            d1Rydq1(By, qs, i)*(4*d3Rxdq3(Bx, qs, i)*d3Rydq3(By, qs, i) + 
  					               3*d1Rydq1(By, qs, i)*d5Rxdq5(Bx, qs, i)))) + 
  					      Power(d1Rxdq1(Bx, qs, i),8)*d5Rydq5(By, qs, i) + 
  					      Power(d1Rxdq1(Bx, qs, i),6)*(45*Power(d2Rxdq2(Bx, qs, i),2)*d3Rydq3(By, qs, i) - 
  					         25*Power(d2Rydq2(By, qs, i),2)*d3Rydq3(By, qs, i) + 
  					         15*d2Rxdq2(Bx, qs, i)*(4*d2Rydq2(By, qs, i)*d3Rxdq3(Bx, qs, i) + 
  					            d1Rydq1(By, qs, i)*d4Rxdq4(Bx, qs, i)) + 
  					         d1Rydq1(By, qs, i)*(10*Power(d3Rxdq3(Bx, qs, i),2) - 13*Power(d3Rydq3(By, qs, i),2) - 
  					            19*d2Rydq2(By, qs, i)*d4Rydq4(By, qs, i)) + 
  					         3*Power(d1Rydq1(By, qs, i),2)*d5Rydq5(By, qs, i)) + 
  					      Power(d1Rxdq1(Bx, qs, i),2)*Power(d1Rydq1(By, qs, i),3)*
  					       (-162*Power(d2Rxdq2(Bx, qs, i),4) - 192*Power(d2Rydq2(By, qs, i),4) + 
  					         163*d1Rydq1(By, qs, i)*Power(d2Rydq2(By, qs, i),2)*d3Rydq3(By, qs, i) + 
  					         3*Power(d2Rxdq2(Bx, qs, i),2)*
  					          (322*Power(d2Rydq2(By, qs, i),2) - 37*d1Rydq1(By, qs, i)*d3Rydq3(By, qs, i)) + 
  					         d1Rydq1(By, qs, i)*d2Rxdq2(Bx, qs, i)*
  					          (-232*d2Rydq2(By, qs, i)*d3Rxdq3(Bx, qs, i) + 
  					            7*d1Rydq1(By, qs, i)*d4Rxdq4(Bx, qs, i)) - 
  					         19*Power(d1Rydq1(By, qs, i),2)*d2Rydq2(By, qs, i)*d4Rydq4(By, qs, i) + 
  					         Power(d1Rydq1(By, qs, i),2)*(4*Power(d3Rxdq3(Bx, qs, i),2) - 
  					            13*Power(d3Rydq3(By, qs, i),2) + d1Rydq1(By, qs, i)*d5Rydq5(By, qs, i))) + 
  					      Power(d1Rxdq1(Bx, qs, i),4)*d1Rydq1(By, qs, i)*
  					       (105*Power(d2Rxdq2(Bx, qs, i),4) + 88*Power(d2Rydq2(By, qs, i),4) + 
  					         138*d1Rydq1(By, qs, i)*Power(d2Rydq2(By, qs, i),2)*d3Rydq3(By, qs, i) - 
  					         Power(d2Rxdq2(Bx, qs, i),2)*(627*Power(d2Rydq2(By, qs, i),2) + 
  					            82*d1Rydq1(By, qs, i)*d3Rydq3(By, qs, i)) + 
  					         2*d1Rydq1(By, qs, i)*d2Rxdq2(Bx, qs, i)*
  					          (-107*d2Rydq2(By, qs, i)*d3Rxdq3(Bx, qs, i) + 
  					            13*d1Rydq1(By, qs, i)*d4Rxdq4(Bx, qs, i)) - 
  					         38*Power(d1Rydq1(By, qs, i),2)*d2Rydq2(By, qs, i)*d4Rydq4(By, qs, i) + 
  					         Power(d1Rydq1(By, qs, i),2)*(17*Power(d3Rxdq3(Bx, qs, i),2) - 
  					            26*Power(d3Rydq3(By, qs, i),2) + 3*d1Rydq1(By, qs, i)*d5Rydq5(By, qs, i)))))/
  					  Power(Power(d1Rxdq1(Bx, qs, i),2) + Power(d1Rydq1(By, qs, i),2),7) - 
  					 (d1Rydq1(By, qs, i)*(Power(d1Rxdq1(Bx, qs, i),6)*
  					       (16*Power(d2Rydq2(By, qs, i),2)*d3Rxdq3(Bx, qs, i) + 
  					         d2Rydq2(By, qs, i)*(42*d2Rxdq2(Bx, qs, i)*d3Rydq3(By, qs, i) + 
  					            9*d1Rydq1(By, qs, i)*d4Rxdq4(Bx, qs, i)) + 
  					         d1Rydq1(By, qs, i)*(16*d3Rxdq3(Bx, qs, i)*d3Rydq3(By, qs, i) + 
  					            14*d2Rxdq2(Bx, qs, i)*d4Rydq4(By, qs, i) + d1Rydq1(By, qs, i)*d5Rxdq5(Bx, qs, i)))\
  					       + Power(d1Rydq1(By, qs, i),5)*(75*Power(d2Rxdq2(Bx, qs, i),3)*d2Rydq2(By, qs, i) - 
  					         25*d1Rydq1(By, qs, i)*Power(d2Rxdq2(Bx, qs, i),2)*d3Rxdq3(Bx, qs, i) - 
  					         5*d2Rxdq2(Bx, qs, i)*(21*Power(d2Rydq2(By, qs, i),3) - 
  					            12*d1Rydq1(By, qs, i)*d2Rydq2(By, qs, i)*d3Rydq3(By, qs, i) + 
  					            Power(d1Rydq1(By, qs, i),2)*d4Rydq4(By, qs, i)) + 
  					         d1Rydq1(By, qs, i)*(45*Power(d2Rydq2(By, qs, i),2)*d3Rxdq3(Bx, qs, i) - 
  					            10*d1Rydq1(By, qs, i)*d2Rydq2(By, qs, i)*d4Rxdq4(Bx, qs, i) + 
  					            d1Rydq1(By, qs, i)*(-10*d3Rxdq3(Bx, qs, i)*d3Rydq3(By, qs, i) + 
  					               d1Rydq1(By, qs, i)*d5Rxdq5(Bx, qs, i)))) + 
  					      Power(d1Rxdq1(Bx, qs, i),2)*Power(d1Rydq1(By, qs, i),3)*
  					       (-766*Power(d2Rxdq2(Bx, qs, i),3)*d2Rydq2(By, qs, i) + 
  					         138*d1Rydq1(By, qs, i)*Power(d2Rxdq2(Bx, qs, i),2)*d3Rxdq3(Bx, qs, i) + 
  					         d2Rxdq2(Bx, qs, i)*(714*Power(d2Rydq2(By, qs, i),3) - 
  					            214*d1Rydq1(By, qs, i)*d2Rydq2(By, qs, i)*d3Rydq3(By, qs, i) + 
  					            4*Power(d1Rydq1(By, qs, i),2)*d4Rydq4(By, qs, i)) + 
  					         d1Rydq1(By, qs, i)*(-82*Power(d2Rydq2(By, qs, i),2)*d3Rxdq3(Bx, qs, i) - 
  					            11*d1Rydq1(By, qs, i)*d2Rydq2(By, qs, i)*d4Rxdq4(Bx, qs, i) + 
  					            d1Rydq1(By, qs, i)*(-4*d3Rxdq3(Bx, qs, i)*d3Rydq3(By, qs, i) + 
  					               3*d1Rydq1(By, qs, i)*d5Rxdq5(Bx, qs, i)))) + 
  					      Power(d1Rxdq1(Bx, qs, i),4)*d1Rydq1(By, qs, i)*
  					       (279*Power(d2Rxdq2(Bx, qs, i),3)*d2Rydq2(By, qs, i) + 
  					         163*d1Rydq1(By, qs, i)*Power(d2Rxdq2(Bx, qs, i),2)*d3Rxdq3(Bx, qs, i) + 
  					         d2Rxdq2(Bx, qs, i)*(-301*Power(d2Rydq2(By, qs, i),3) - 
  					            232*d1Rydq1(By, qs, i)*d2Rydq2(By, qs, i)*d3Rydq3(By, qs, i) + 
  					            23*Power(d1Rydq1(By, qs, i),2)*d4Rydq4(By, qs, i)) + 
  					         d1Rydq1(By, qs, i)*(-111*Power(d2Rydq2(By, qs, i),2)*d3Rxdq3(Bx, qs, i) + 
  					            8*d1Rydq1(By, qs, i)*d2Rydq2(By, qs, i)*d4Rxdq4(Bx, qs, i) + 
  					            d1Rydq1(By, qs, i)*(22*d3Rxdq3(Bx, qs, i)*d3Rydq3(By, qs, i) + 
  					               3*d1Rydq1(By, qs, i)*d5Rxdq5(Bx, qs, i)))) - 
  					      Power(d1Rxdq1(Bx, qs, i),7)*(3*Power(d3Rydq3(By, qs, i),2) + 
  					         4*d2Rydq2(By, qs, i)*d4Rydq4(By, qs, i) + d1Rydq1(By, qs, i)*d5Rydq5(By, qs, i)) + 
  					      Power(d1Rxdq1(Bx, qs, i),3)*Power(d1Rydq1(By, qs, i),2)*
  					       (-192*Power(d2Rxdq2(Bx, qs, i),4) - 162*Power(d2Rydq2(By, qs, i),4) - 
  					         22*d1Rydq1(By, qs, i)*Power(d2Rydq2(By, qs, i),2)*d3Rydq3(By, qs, i) + 
  					         14*Power(d2Rxdq2(Bx, qs, i),2)*
  					          (69*Power(d2Rydq2(By, qs, i),2) + d1Rydq1(By, qs, i)*d3Rydq3(By, qs, i)) + 
  					         2*d1Rydq1(By, qs, i)*d2Rxdq2(Bx, qs, i)*
  					          (54*d2Rydq2(By, qs, i)*d3Rxdq3(Bx, qs, i) - 19*d1Rydq1(By, qs, i)*d4Rxdq4(Bx, qs, i))
  					           + 26*Power(d1Rydq1(By, qs, i),2)*d2Rydq2(By, qs, i)*d4Rydq4(By, qs, i) + 
  					         Power(d1Rydq1(By, qs, i),2)*(-26*Power(d3Rxdq3(Bx, qs, i),2) + 
  					            17*Power(d3Rydq3(By, qs, i),2) - 3*d1Rydq1(By, qs, i)*d5Rydq5(By, qs, i))) + 
  					      d1Rxdq1(Bx, qs, i)*Power(d1Rydq1(By, qs, i),4)*
  					       (88*Power(d2Rxdq2(Bx, qs, i),4) + 105*Power(d2Rydq2(By, qs, i),4) - 
  					         105*d1Rydq1(By, qs, i)*Power(d2Rydq2(By, qs, i),2)*d3Rydq3(By, qs, i) + 
  					         Power(d2Rxdq2(Bx, qs, i),2)*(-627*Power(d2Rydq2(By, qs, i),2) + 
  					            101*d1Rydq1(By, qs, i)*d3Rydq3(By, qs, i)) + 
  					         d1Rydq1(By, qs, i)*d2Rxdq2(Bx, qs, i)*
  					          (242*d2Rydq2(By, qs, i)*d3Rxdq3(Bx, qs, i) - 
  					            19*d1Rydq1(By, qs, i)*d4Rxdq4(Bx, qs, i)) + 
  					         15*Power(d1Rydq1(By, qs, i),2)*d2Rydq2(By, qs, i)*d4Rydq4(By, qs, i) - 
  					         Power(d1Rydq1(By, qs, i),2)*(13*Power(d3Rxdq3(Bx, qs, i),2) - 
  					            10*Power(d3Rydq3(By, qs, i),2) + d1Rydq1(By, qs, i)*d5Rydq5(By, qs, i))) - 
  					      Power(d1Rxdq1(Bx, qs, i),5)*(-13*Power(d2Rydq2(By, qs, i),4) - 
  					         83*d1Rydq1(By, qs, i)*Power(d2Rydq2(By, qs, i),2)*d3Rydq3(By, qs, i) + 
  					         87*Power(d2Rxdq2(Bx, qs, i),2)*
  					          (Power(d2Rydq2(By, qs, i),2) + d1Rydq1(By, qs, i)*d3Rydq3(By, qs, i)) + 
  					         d1Rydq1(By, qs, i)*d2Rxdq2(Bx, qs, i)*
  					          (134*d2Rydq2(By, qs, i)*d3Rxdq3(Bx, qs, i) + 
  					            19*d1Rydq1(By, qs, i)*d4Rxdq4(Bx, qs, i)) - 
  					         7*Power(d1Rydq1(By, qs, i),2)*d2Rydq2(By, qs, i)*d4Rydq4(By, qs, i) + 
  					         Power(d1Rydq1(By, qs, i),2)*(13*Power(d3Rxdq3(Bx, qs, i),2) - 
  					            4*Power(d3Rydq3(By, qs, i),2) + 3*d1Rydq1(By, qs, i)*d5Rydq5(By, qs, i)))))/
  					  Power(Power(d1Rxdq1(Bx, qs, i),2) + Power(d1Rydq1(By, qs, i),2),7);
	}

	
	ROS_INFO("Bezier[0] = (%.3f, %.3f), qs[0]=%.3f, cs[0][0]=%.3f",
         R[0][0], R[0][1], qs[0], cs[0][0]);


	//全探索
	searchP(x_old);

	//制御入力を計算し、それらをルンゲクッタ法で更新
	getInputValue.getU(x_old, sr.j);
	integrator.step(q_map, qdot_map, u1, u2);
	getInputValue.ddrungeKutta(x_d, x_dd);
	getInputValue.rungeKutta(x_old, x_d);


	// 各車両へ steering コマンドと車輪のトルクコマンドを送信
    vehicle1.publishSteeringCommand(Q_phiF/2.0, Q_phiF/2.0, Q_phiR/2.0, Q_phiR/2.0);
    vehicle1.publishWheelCommand(torque_front[0], torque_front[1], torque_rear[0], torque_rear[1]);

	//誤差平均
	roop_sum ++;
	d_sum += std::abs(sr.d);


	//Gazeboのフィードバックをもとに計算
	while(ros::ok()) {
		//ROSのコールバックを処理
		double current_time = ros::Time::now().toSec();
		ros::spinOnce();
	
		//部分探索
		searchPP(x_old);

		roop_sum ++;
		d_sum += std::abs(sr.d);
		d_ave = d_sum / roop_sum;
	
		//制御入力を計算し、それらをルンゲクッタ法で更新
		getInputValue.getU(x_old, sr.j);
		integrator.step(q_map, qdot_map, u1, u2);
		getInputValue.ddrungeKutta(x_d, x_dd);
		getInputValue.rungeKutta(x_old, x_d);

		
	

		//デバッグ用ログ出力
    	// ROS_INFO_THROTTLE(1.0, "x_old = [%.3f, %.3f, %.3f, %.3f, %.3f]",
    	//     x_old[0], x_old[1], x_old[2], x_old[3], x_old[4]);

    	ROS_INFO_THROTTLE(0.01, "sr: j=%d, Psx=%.3f, Psy=%.3f, d=%.3f, Cs=%.6f, dCs1=%.6f, dCs2=%.6f, dCs3=%.6f, d_ave=%.6f, thetap=%.6f",
    	    sr.j, sr.Psx, sr.Psy, sr.d, sr.Cs, sr.Cs1, sr.Cs2, sr.Cs3, d_ave, Thetap);
		
		
		// ROS_INFO_THROTTLE(1.0, "x_d: x_d=%.3f, y_d=%.3f, theta_d=%.3f",
		//     x_d[1], x_d[2],x_d[3]);

		ROS_INFO_THROTTLE(0.01, "Q: Q_phir=%.3f, Q_phif=%.3f, Q_varphir=%.3f, Q_varphif=%.3f",
		    Q_phiR, Q_phiF, Q_varphiR, Q_varphiF);

		ROS_INFO_THROTTLE(0.01, "nu: nu1=%.3f, nu2=%.3f, nu3=%.3f",
		    nu1, nu2, nu3);

		// ROS_INFO_THROTTLE(0.01, "Torque: front_left=%.3f, front_right=%.3f, rear_left=%.3f, rear_right=%.3f\n",
		//     torque_front[0], torque_front[1], torque_rear[0], torque_rear[1]);
		


		// ROS_INFO_THROTTLE(1.0, "x_dd: x_dd=%.3f, y_dd=%.3f, theta_dd=%.3f",
		//     x_dd[1], x_dd[2],x_dd[3]);
	
		// 各車両へ steering コマンドと車輪のトルクコマンドを送信
		vehicle1.publishSteeringCommand(Q_phiF/2.0, Q_phiF/2.0, Q_phiR/2.0, Q_phiR/2.0);
		vehicle1.publishWheelCommand(torque_front[0], torque_front[1], torque_rear[0], torque_rear[1]);

		// ループレートを維持
		logger.logData();
		loop_rate.sleep();
	}
	return 0;
}



//Ps探索
Search searchP(std::vector<double>& x) {
	int i;
	double dot = 0.0;
	double dist = 0.0;
	double dist0 = DBL_MAX;

	for (i = 0; i < Q; i++) {
		dot = (x[1] - R[i][0]) * (dRdq[i][0] / norm(Bx, By, qs, i)) + (x[2] - R[i][1]) * (dRdq[i][1] / norm(Bx, By, qs, i));

		dist = sqrt(pow((x[1] - R[i][0]), 2) + pow((x[2] - R[i][1]), 2));

		if (-0.0001 < dot && dot < 0.0001) {
			if (dist < dist0) {
				dist0 = dist;
				sr.Psx = R[i][0];
				sr.Psy = R[i][1];
				sr.d = (x[1] - R[i][0]) * (-(dRdq[i][1] / norm(Bx, By, qs, i))) + (x[2] - R[i][1]) * (dRdq[i][0] / norm(Bx, By, qs, i));
				sr.Cs = cs[i][0];
				sr.Cs1 = cs[i][1];
				sr.Cs2 = cs[i][2];
				sr.Cs3 = cs[i][3];
				sr.j = i;
			}
		}
	}
	return(sr);
}

Search searchPP(std::vector<double>& x) {
	int i;
	double dot = 0.0;
	double dist = 0.0;
	double dist0= DBL_MAX;

	if (sr.j < PSdist) {
		for (i = 0; i < sr.j + PSdist; i++) {
			dot = (x[1] - R[i][0]) * (dRdq[i][0] / norm(Bx, By, qs, i)) + (x[2] - R[i][1]) * (dRdq[i][1] / norm(Bx, By, qs, i));

			dist = sqrt(pow((x[1] - R[i][0]), 2) + pow((x[2] - R[i][1]), 2));

			if (-0.001 < dot && dot < 0.001) {
				if (dist < dist0) {
					dist0 = dist;
					sr.Psx = R[i][0];
					sr.Psy = R[i][1];
					sr.d = (x[1] - R[i][0]) * (-(dRdq[i][1] / norm(Bx, By, qs, i))) + (x[2] - R[i][1]) * (dRdq[i][0] / norm(Bx, By, qs, i));
					sr.Cs = cs[i][0];
					sr.Cs1 = cs[i][1];
					sr.Cs2 = cs[i][2];
					sr.Cs3 = cs[i][3];
					sr.j = i;
				}
		}
	}
	}else if (sr.j > Q - PSdist) {
		for (i = sr.j - PSdist; i < Q; i++) {
			dot = (x[1] - R[i][0]) * (dRdq[i][0] / norm(Bx, By, qs, i)) + (x[2] - R[i][1]) * (dRdq[i][1] / norm(Bx, By, qs, i));

			dist = sqrt(pow((x[1] - R[i][0]), 2) + pow((x[2] - R[i][1]), 2));

			if (-0.001 < dot && dot < 0.001) {
				if (dist < dist0) {
					dist0 = dist;
					sr.Psx = R[i][0];
					sr.Psy = R[i][1];
					sr.d = (x[1] - R[i][0]) * (-(dRdq[i][1] / norm(Bx, By, qs, i))) + (x[2] - R[i][1]) * (dRdq[i][0] / norm(Bx, By, qs, i));
					sr.Cs = cs[i][0];
					sr.Cs1 = cs[i][1];
					sr.Cs2 = cs[i][2];
					sr.Cs3 = cs[i][3];
					sr.j = i;
				}
			}
		}
	}
	else {
		for (i = sr.j - PSdist; i < sr.j+PSdist; i++) {
			dot = (x[1] - R[i][0]) * (dRdq[i][0] / norm(Bx, By, qs, i)) + (x[2] - R[i][1]) * (dRdq[i][1] / norm(Bx, By, qs, i));

			dist = sqrt(pow((x[1] - R[i][0]), 2) + pow((x[2] - R[i][1]), 2));

			if (-0.001 < dot && dot < 0.001) {
				if (dist < dist0) {
					dist0 = dist;
					sr.Psx = R[i][0];
					sr.Psy = R[i][1];
					sr.d = (x[1] - R[i][0]) * (-(dRdq[i][1] / norm(Bx, By, qs, i))) + (x[2] - R[i][1]) * (dRdq[i][0] / norm(Bx, By, qs, i));
					sr.Cs = cs[i][0];
					sr.Cs1 = cs[i][1];
					sr.Cs2 = cs[i][2];
					sr.Cs3 = cs[i][3];
					sr.j = i;
				}
			}
		}
	}
	return(sr);
}