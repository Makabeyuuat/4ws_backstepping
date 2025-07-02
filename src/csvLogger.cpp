#include "csvLogger.hpp"
#include "initial.hpp"


CSVLogger::CSVLogger(const std::string& dir, int close_threshold)
  : closed_(false)
  , close_threshold_(close_threshold){
    // タイムスタンプを生成してファイル名を作成
    std::string ts = makeTimeStamp();
    std::string filename = dir + "/data_log_" + ts + ".csv";

    csv_.open(filename, std::ios::out);
    if (!csv_.is_open()) {
      ROS_ERROR("CSVLogger: failed to open '%s'", filename.c_str());
      return;
  }
  // ヘッダー行
  csv_ << "t,x,y,theta,phiR,phiF,thetap"
       << "sr_j,Psx,Psy,d,Cs,Cs1,Cs2,Cs3,d_ave,"
       << "x_d1,x_d2,x_d3,"
       << "nu1,nu2,nu3"
       << "torque_fl,torque_fr,torque_rl,torque_rr,"
       << "lamda1,lambda2,lambda3,lambda4,"
        << "Q_phiR, Q_varphiR, Q_phiF, Q_varphiF"
       << "\n";
}

CSVLogger::~CSVLogger() {
  if (csv_.is_open()) csv_.close();
}

void CSVLogger::logData() {
  if (!csv_.is_open()) return;

  double t = ros::Time::now().toSec();
  csv_ << x_old[0] << ',' << x_old[1] << ',' << x_old[2] << ','
       << x_old[3] << ',' << x_old[4] << ','<< x_old[5] << ',' << Thetap << ','
       // sr
       << sr.j       << ',' << sr.Psx  << ',' << sr.Psy  << ','
       << sr.d       << ',' << sr.Cs   << ',' << sr.Cs1   << ','
       << sr.Cs2     << ',' << sr.Cs3  << ',' << d_ave    << ','
       // x_d
       << x_d[1]     << ',' << x_d[2] << ',' << x_d[3]    << ','
       // nu
       << nu1        << ',' << nu2     << ','<< nu3     << ','
       // torque
       << torque_front[0] << ',' << torque_front[1] << ','
       << torque_rear[0]  << ',' << torque_rear[1]  << ','
       // lamda
       << lamda_data(0) << ',' << lamda_data(1) << ','
       << lamda_data(2)  << ',' << lamda_data(3)  << ','
       << Q_phiR << ',' << Q_varphiR << ',' << Q_phiF<< ',' << Q_varphiF
       << "\n";
  csv_.flush();
  // 閾値に到達したらクローズ
  if (sr.j == close_threshold_) {
    csv_.close();
    closed_ = true;
    ROS_INFO("CSVLogger: closed CSV file at sr.j = %d", sr.j);
  }
}

std::string CSVLogger::makeTimeStamp() {
  std::time_t t = std::time(nullptr);
  std::tm    tm = *std::localtime(&t);
  std::ostringstream oss;
  oss << std::put_time(&tm, "%Y%m%d_%H%M%S");
  return oss.str();
}
