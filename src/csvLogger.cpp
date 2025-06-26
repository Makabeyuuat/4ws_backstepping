#include "csvLogger.hpp"
#include "initial.hpp"


CSVLogger::CSVLogger(const std::string& filename, int close_threshold)
  : closed_(false)
  , close_threshold_(close_threshold){
  csv_.open(filename, std::ios::out);
  if (!csv_.is_open()) {
    ROS_ERROR("CSVLogger: failed to open '%s'", filename.c_str());
    return;
  }
  // ヘッダー行
  csv_ << "t,x,y,theta,phi,"
       << "sr_j,Psx,Psy,d,Cs,Cs1,Cs2,Cs3,d_ave,"
       << "x_d1,x_d2,x_d3,"
       << "nu1,nu2,"
       << "torque_fl,torque_fr,torque_rl,torque_rr,Q_phi"
       << "\n";
}

CSVLogger::~CSVLogger() {
  if (csv_.is_open()) csv_.close();
}

void CSVLogger::logData() {
  if (!csv_.is_open()) return;

  double t = ros::Time::now().toSec();
  csv_ << x_old[0] << ',' << x_old[1] << ',' << x_old[2] << ','
       << x_old[3] << ',' << x_old[4] << ','
       // sr
       << sr.j       << ',' << sr.Psx  << ',' << sr.Psy  << ','
       << sr.d       << ',' << sr.Cs   << ',' << sr.Cs1   << ','
       << sr.Cs2     << ',' << sr.Cs3  << ',' << d_ave    << ','
       // x_d
       << x_d[1]     << ',' << x_d[2] << ',' << x_d[3]    << ','
       // nu
       << nu1        << ',' << nu2     << ','
       // torque
       << torque_front[0] << ',' << torque_front[1] << ','
       << torque_rear[0]  << ',' << torque_rear[1]  << ','
       << Q_phi
       << "\n";
  csv_.flush();
  // 閾値に到達したらクローズ
  if (sr.j == close_threshold_) {
    csv_.close();
    closed_ = true;
    ROS_INFO("CSVLogger: closed CSV file at sr.j = %d", sr.j);
  }
}
