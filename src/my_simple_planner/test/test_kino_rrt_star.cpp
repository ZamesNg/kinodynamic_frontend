#include <iostream>
#include <memory>

#include <ros/ros.h>
#include <nav_msgs/Path.h>
#include <visualization_msgs/Marker.h>

#include "local_perception/sdf_map.h"
#include "kinodyn_rrt_star/kinodyn_rrt_star.h"

using namespace std;

shared_ptr<SDFMap> sdf_map;
KinodynRRTStarPlanner planner;
nav_msgs::Odometry odom;
ros::Publisher vis_traj_pub;

void recvGoalCallback(const geometry_msgs::PoseStamped &wp);
void recvOdomCallback(const nav_msgs::OdometryConstPtr &odom_);
void visualizeESDFGrid(const ros::TimerEvent &event);

int main(int argc, char **argv)
{
  ros::init(argc, argv, "test_esdf");
  ros::NodeHandle nh("~");

  ros::Subscriber pts_sub = nh.subscribe("/test_planner/goal", 1, recvGoalCallback);
  ros::Subscriber odom_sub = nh.subscribe("/test_planner/odom", 1, recvOdomCallback);
  vis_traj_pub = nh.advertise<visualization_msgs::Marker>("RRTStar_path_vis", 1);

  sdf_map.reset(new SDFMap);
  sdf_map->initMap(nh);
  planner.initPlanner(nh, sdf_map);

  ros::Timer vis_esdf_timer = nh.createTimer(ros::Duration(1.0), visualizeESDFGrid);
  ros::spin();
  return 0;
}

void visualizeESDFGrid(const ros::TimerEvent &event)
{
  if (!sdf_map->hasDepthObservation())
    return;

  sdf_map->publishMap();
  sdf_map->publishMapInflate();
  sdf_map->publishUpdateRange();
  sdf_map->publishESDF();
}

void recvGoalCallback(const geometry_msgs::PoseStamped &wp)
{
  auto pos = odom.pose.pose.position;
  Eigen::Vector3d start = {pos.x, pos.y, pos.z};
  pos = wp.pose.position;
  Eigen::Vector3d goal = {pos.x, pos.y, pos.z};

  ROS_WARN("recv goal start planning!");
  planner.searchTraj(start, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), goal);
  ROS_WARN("finish planning!");

  std::vector<Eigen::Matrix<double, 6, 3>> coeffs = planner.getTrajCoeff();
  std::vector<double> intervals = planner.getTrajInterval();

  double T;
  Eigen::Matrix<double, 6, 3> coeff;
  Eigen::Matrix<double, 6, 1> nature_bais;
  Eigen::Vector3d position;

  visualization_msgs::Marker Line;

  Line.header.frame_id = "world";
  Line.header.stamp = ros::Time::now();
  Line.ns = "test_rrt_star_node/traj";
  Line.action = visualization_msgs::Marker::ADD;
  Line.pose.orientation.w = 1.0;
  Line.type = visualization_msgs::Marker::LINE_STRIP;
  Line.scale.x = 0.2;

  Line.color.r = 1.0;
  Line.color.g = 0.0;
  Line.color.b = 0.0;
  Line.color.a = 1.0;

  ROS_WARN("%ld picese trajectory", coeffs.size());
  for (size_t i = 0; i < coeffs.size(); ++i)
  {
    coeff = coeffs[i];
    T = intervals[i];

    for (double t = 0.0; t < T; t += 0.05)
    {
      nature_bais << 1.0, t, t * t, t * t * t,
          t * t * t * t, t * t * t * t * t;

      position = coeff.transpose() * nature_bais;
      geometry_msgs::Point pt;
      pt.x = position(0);
      pt.y = position(1);
      pt.z = position(2);
      Line.points.push_back(pt);
    }
    nature_bais << 1.0, T, T * T, T * T * T,
        T * T * T * T, T * T * T * T * T;

    position = coeff.transpose() * nature_bais;
    geometry_msgs::Point pt;
    pt.x = position(0);
    pt.y = position(1);
    pt.z = position(2);
    Line.points.push_back(pt);
  }
  vis_traj_pub.publish(Line);
}

void recvOdomCallback(const nav_msgs::OdometryConstPtr &odom_)
{
  odom = *odom_;
}