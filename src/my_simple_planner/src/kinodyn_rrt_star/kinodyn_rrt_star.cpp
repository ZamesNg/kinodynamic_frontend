#include <memory>
#include <math.h>

#include <ros/ros.h>

#include "kinodyn_rrt_star/kinodyn_rrt_star.h"

double StateNode::calcOptimalTrajWithFullState_(DroneState_t start, DroneState_t end, double &optimal_T, Eigen::Matrix<double, 6, 3> &coeff)
{
  Eigen::Vector3d p0 = start.pos, v0 = start.vel, a0 = start.acc;
  Eigen::Vector3d pf = end.pos, vf = end.vel, af = end.acc;
  Eigen::Matrix<double, 6, 6> m;

  double c4 = -9.0 * a0.squaredNorm() - 9.0 * af.squaredNorm() + 6.0 * a0.transpose() * af;
  double c3 = (af.transpose() * (96.0 * v0 + 144.0 * vf) - a0.transpose() * (96.0 * vf + 144.0 * v0))(0);
  double c2 = -576.0 * v0.squaredNorm() - 1008.0 * v0.transpose() * vf - 576.0 * vf.squaredNorm() +
              -360.0 * (p0 - pf).transpose() * (a0 - af);
  double c1 = -2880.0 * (v0 + vf).transpose() * (p0 - pf);
  double c0 = -3600.0 * (p0 - pf).squaredNorm();

  m << 0.0, 0.0, 0.0, 0.0, 0.0, -c0,
      1.0, 0.0, 0.0, 0.0, 0.0, -c1,
      0.0, 1.0, 0.0, 0.0, 0.0, -c2,
      0.0, 0.0, 1.0, 0.0, 0.0, -c3,
      0.0, 0.0, 0.0, 1.0, 0.0, -c4,
      0.0, 0.0, 0.0, 0.0, 1.0, 0.0;

  double cost = StateNode::kErrorCost, optimal_cost = StateNode::kErrorCost;
  double T = -1.0;
  optimal_T = -1.0;

  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> eigen_values;
  eigen_values = m.eigenvalues();

  for (int i = 0; i < 6; ++i)
  {
    T = std::real(eigen_values(i));
    double img = std::imag(eigen_values(i));
    double freq = 1.0 / T;

    if (T <= 0 || std::abs(img) >= 1e-16)
      continue;

    Eigen::Vector3d k_alpha, k_beta, k_gamma;
    k_gamma = freq * (3.0 * af - 9.0 * a0 + freq * (-36.0 * v0 - 24.0 * vf + freq * 60.0 * (pf - p0)));
    k_beta = freq * freq * (-24.0 * af + 36.0 * a0 + freq * (192.0 * v0 + 168.0 * vf + freq * 360.0 * (p0 - pf)));
    k_alpha = freq * freq * freq * (60.0 * (af - a0) + freq * (-360.0 * (v0 + vf) + freq * 720.0 * (pf - p0)));

    cost = T * (1.0 + k_gamma.squaredNorm() +
                T * ((k_beta.transpose() * k_gamma)(0) +
                     T * ((k_beta.squaredNorm() + (k_alpha.transpose() * k_gamma)(0)) / 3.0 +
                          T * ((k_alpha.transpose() * k_beta)(0) / 4.0 +
                               T * k_alpha.squaredNorm() / 20.0))));
    if (cost <= optimal_cost)
    {
      optimal_cost = cost;
      optimal_T = T;
      coeff.row(3) = k_gamma / 6.0;
      coeff.row(4) = k_beta / 24.0;
      coeff.row(5) = k_alpha / 120.0;
    }
  }

  coeff.row(0) = p0;
  coeff.row(1) = v0;
  coeff.row(2) = a0 / 2.0;

  if (optimal_T == -1.0)
    ROS_ERROR("from point(%.4f, %.4f, %.4f) to point(%.4f, %.4f, %.4f) optimal_T cannot be found!", p0(0), p0(1), p0(2), pf(0), pf(1), pf(2));

  return optimal_cost;
}

double StateNode::calcOptimalTrajWithPartialState_(DroneState_t start, Eigen::Vector3d pf, double &optimal_T, Eigen::Matrix<double, 6, 3> &coeff)
{
  Eigen::Vector3d p0 = start.pos, v0 = start.vel, a0 = start.acc;
  Eigen::Matrix<double, 6, 6> m;

  double c4 = -5.0 * a0.squaredNorm();
  double c3 = -40.0 * a0.transpose() * v0;
  double c2 = -60.0 * (v0.squaredNorm() + a0.transpose() * (p0 - pf));
  double c1 = -160.0 * v0.transpose() * (p0 - pf);
  double c0 = -100.0 * (p0 - pf).squaredNorm();

  m << 0.0, 0.0, 0.0, 0.0, 0.0, -c0,
      1.0, 0.0, 0.0, 0.0, 0.0, -c1,
      0.0, 1.0, 0.0, 0.0, 0.0, -c2,
      0.0, 0.0, 1.0, 0.0, 0.0, -c3,
      0.0, 0.0, 0.0, 1.0, 0.0, -c4,
      0.0, 0.0, 0.0, 0.0, 1.0, 0.0;

  double cost = StateNode::kErrorCost, optimal_cost = StateNode::kErrorCost;
  double T = -1.0;
  optimal_T = -1.0;

  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> eigen_values;
  eigen_values = m.eigenvalues();

  for (int i = 0; i < 6; ++i)
  {
    T = std::real(eigen_values(i));
    double img = std::imag(eigen_values(i));
    double freq = 1.0 / T;

    if (T <= 0 || std::abs(img) >= 1e-16)
      continue;

    Eigen::Vector3d k_alpha, k_beta, k_gamma;
    k_gamma = (-5.0 * a0 + (-10.0 * v0 - 10.0 * (p0 - pf) * freq) * freq) * freq;
    k_beta = -2.0 * freq * k_gamma;
    k_alpha = -freq * k_beta;

    // cost = (-k_gamma / 10.0).squaredNorm() * 20.0 * T + T;
    cost = T * (1.0 + k_gamma.squaredNorm() +
                T * ((k_beta.transpose() * k_gamma)(0) +
                     T * ((k_beta.squaredNorm() + (k_alpha.transpose() * k_gamma)(0)) / 3.0 +
                          T * ((k_alpha.transpose() * k_beta)(0) / 4.0 +
                               T * k_alpha.squaredNorm() / 20.0))));
    if (cost <= optimal_cost)
    {
      optimal_cost = cost;
      optimal_T = T;
      coeff.row(3) = k_gamma / 6.0;
      coeff.row(4) = k_beta / 24.0;
      coeff.row(5) = k_alpha / 120.0;
    }
  }

  coeff.row(0) = p0;
  coeff.row(1) = v0;
  coeff.row(2) = a0 / 2.0;

  if (optimal_T == -1.0)
    ROS_ERROR("from point(%.4f, %.4f, %.4f) to point(%.4f, %.4f, %.4f) optimal_T cannot be found!", p0(0), p0(1), p0(2), pf(0), pf(1), pf(2));

  return optimal_cost;
}

void StateNode::updateStateFromCoeff()
{
  Eigen::Matrix<double, 6, 1> nature_bias;
  nature_bias << 1.0, T_, T_ * T_, T_ * T_ * T_,
      T_ * T_ * T_ * T_, T_ * T_ * T_ * T_ * T_;

  state_.pos = polynomial_coeff.transpose() * nature_bias;

  nature_bias
      << 0.0,
      1.0, 2.0 * T_, 3.0 * T_ * T_,
      4.0 * T_ * T_ * T_, 5.0 * T_ * T_ * T_ * T_;

  state_.vel = polynomial_coeff.transpose() * nature_bias;

  nature_bias << 0.0, 0.0, 2.0, 6.0 * T_,
      12.0 * T_ * T_, 20.0 * T_ * T_ * T_;

  state_.acc = polynomial_coeff.transpose() * nature_bias;
}

bool StateNode::setParent(Ptr parent, double cost_to_parent, double T, Eigen::Matrix<double, 6, 3> coeff)
{
  parent_ = parent;
  cost_to_parent_ = cost_to_parent;
  cost_ = parent->getCost() + cost_to_parent_;
  T_ = T;
  polynomial_coeff = coeff;

  updateStateFromCoeff();
  return true;
}

bool StateNode::setGoal(DroneState_t goal_state)
{
  double T;
  Eigen::Matrix<double, 6, 3> coeff;
  heuristic_ = calcOptimalTrajWithFullState_(state_, goal_state, T, coeff);
  return true;
}

void KinodynRRTStarPlanner::initPlanner(ros::NodeHandle &nh, shared_ptr<SDFMap> map)
{
  // get parameter
  nh.param("planner/goal_tolerance", goal_tolerance, 0.1);
  nh.param("planner/resolution", resolution, 0.05);

  Eigen::Vector3d map_size;
  nh.param("planner/map_size/x", map_size(0), 10.0);
  nh.param("planner/map_size/y", map_size(1), 10.0);
  nh.param("planner/map_size/z", map_size(2), 10.0);
  map_size_(0, 0) = map_size(0) / 2.0;
  map_size_(1, 1) = map_size(1) / 2.0;
  map_size_(2, 2) = map_size(2) / 2.0;

  nh.param("planner/origin/x", map_origin_(0), 0.0);
  nh.param("planner/origin/y", map_origin_(1), 0.0);
  nh.param("planner/origin/z", map_origin_(2), 0.0);

  nh.param("planner/origin/z", step_size_, 1.0);

  ROS_INFO("------------Kinodynamic RRT Star Parameter List------------");
  ROS_INFO("goal_tolerance:%.4fm", goal_tolerance);
  ROS_INFO("resolution:%.4fs", resolution);
  ROS_INFO("map_size:(%.4f, %.4f, %.4f)", map_size(0), map_size(1), map_size(2));
  ROS_INFO("map_origin:(%.4f, %.4f, %.4f)", map_origin_(0), map_origin_(1), map_origin_(2));

  // vis_timer = nh.createTimer(ros::Duration(0.5), &KinodynRRTStarPlanner::visualizeTraj, this);
  vis_traj_library_pub_ = nh.advertise<visualization_msgs::MarkerArray>("/planner/trajectory_lib", 1);

  vis_sample_trajectory_.header.frame_id = "world";
  vis_sample_trajectory_.header.stamp = ros::Time::now();
  vis_sample_trajectory_.ns = "planner/trajectory_library";
  vis_sample_trajectory_.action = visualization_msgs::Marker::ADD;
  vis_sample_trajectory_.pose.orientation.w = 1.0;
  vis_sample_trajectory_.type = visualization_msgs::Marker::LINE_STRIP;
  vis_sample_trajectory_.scale.x = 0.02;

  vis_sample_trajectory_.color.r = 0.0;
  vis_sample_trajectory_.color.g = 1.0;
  vis_sample_trajectory_.color.b = 0.0;
  vis_sample_trajectory_.color.a = 0.2;

  vis_final_trajectory_.header.frame_id = "world";
  vis_final_trajectory_.header.stamp = ros::Time::now();
  vis_final_trajectory_.ns = "planner/trajectory_library";
  vis_final_trajectory_.action = visualization_msgs::Marker::ADD;
  vis_final_trajectory_.pose.orientation.w = 1.0;
  vis_final_trajectory_.type = visualization_msgs::Marker::LINE_STRIP;
  vis_final_trajectory_.scale.x = 0.08;

  vis_final_trajectory_.color.r = 1.0;
  vis_final_trajectory_.color.g = 0.0;
  vis_final_trajectory_.color.b = 0.0;
  vis_final_trajectory_.color.a = 0.7;

  vis_sample_points_.header.frame_id = "world";
  vis_sample_points_.header.stamp = ros::Time::now();
  vis_sample_points_.ns = "planner/trajectory_library";
  vis_sample_points_.action = visualization_msgs::Marker::ADD;
  vis_sample_points_.pose.orientation.w = 1.0;
  vis_sample_points_.type = visualization_msgs::Marker::POINTS;
  vis_sample_points_.scale.x = 0.05;
  vis_sample_points_.scale.y = 0.05;
  vis_sample_points_.scale.z = 0.05;

  vis_sample_points_.color.r = 1.0;
  vis_sample_points_.color.g = 0.0;
  vis_sample_points_.color.b = 0.0;
  vis_sample_points_.color.a = 0.8;

  // init planner
  reach_goal_ = false;
  sample_node_ = nullptr;
  total_cost_ = -1.0;
  kd_tree = kd_create(3);
  map_ = map;
  goal_ = new StateNode(Eigen::Vector3d::Zero());
  goal_->setCost(StateNode::kErrorCost);
}
// always find a trajectory but may not reach the goal
// when time limit is reached, will return a trajectory
// that closest to the goal.
bool KinodynRRTStarPlanner::searchTraj(Eigen::Vector3d start_pos, Eigen::Vector3d start_vel, Eigen::Vector3d start_acc, Eigen::Vector3d end_pos)
{
  goal_->setState(end_pos, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
  // nothing update
  if (!map_->hasDepthObservation())
  {
    ROS_WARN("The ESDF map does not have any observation!");
    return false;
  }

  if (start_)
    delete start_;

  start_ = new StateNode(start_pos);
  start_->setState(start_pos, start_vel, start_acc);
  start_->setGoal(goal_->getState());
  updateNode(start_);

  reach_goal_ = false;

  // TODO: add time limits
  while (!reach_goal_)
  {
    samplePos();
    // Eigen::Vector3d end = sample_node_->getState().pos;
    // ROS_INFO("---------------- new pos(%.4f, %.4f, %.4f) is sampled ----------------", end(0), end(1), end(2));
    // steer()
    double min_cost = StateNode::kErrorCost;
    double optimal_T = -1.0;
    StateNode::Ptr parent = nullptr;
    double T, cost = StateNode::kErrorCost;

    Eigen::Matrix<double, 6, 3> coeff, optimal_coeff;

    nearestBackwardNeighbors();
    for (auto near = neighbors_.cbegin(); near != neighbors_.cend(); ++near)
    {
      // from near to sample_node
      cost = sample_node_->calcOptimalTrajWithPartialState((*near)->getState(), T, coeff);
      if (checkCollision(coeff, T) && cost < min_cost)
      {
        min_cost = cost;
        parent = (*near);
        optimal_T = T;
        optimal_coeff = coeff;
      }
    }

    if (parent && min_cost < StateNode::kErrorCost)
    {
      // Eigen::Vector3d start = parent->getState().pos;
      // ROS_WARN("Find a collision free path from (%.4f, %.4f, %.4f) to (%.4f, %.4f, %.4f) with time:%.4fs and cost:%.4f!", start(0), start(1), start(2), end(0), end(1), end(2), optimal_T, min_cost);
      sample_node_->setParent(parent, min_cost, optimal_T, optimal_coeff);
      sample_node_->setGoal(goal_->getState());
      updateNode(sample_node_);
    }
    else
    {
      // Eigen::Vector3d pos = sample_node_->getState().pos;
      // ROS_WARN("Can not find a feasiable free path from pos:(%.4f, %.4f, %.4f), this state will be dropped!", pos(0), pos(1), pos(2));
      continue;
    }

    // rewrite()
    nearestForwardNeighbors();
    for (auto near = neighbors_.begin(); near != neighbors_.end(); ++near)
    {
      // from sample_node to near
      cost = (*near)->calcOptimalTrajWithFullState(sample_node_->getState(), T, coeff);
      if (checkCollision(coeff, T) && cost + sample_node_->getCost() < (*near)->getCost())
      {
        (*near)->setParent(sample_node_, cost, T, coeff);
        // Eigen::Vector3d start = sample_node_->getState().pos;
        // end = (*near)->getState().pos;
        // ROS_WARN("rewrite a path from (%.4f, %.4f, %.4f) to (%.4f, %.4f, %.4f) with time:%.4fs and cost:%.4f!", start(0), start(1), start(2), end(0), end(1), end(2), T, cost);
      }
    }

    cost = goal_->calcOptimalTrajWithPartialState(sample_node_->getState(), T, coeff);
    if (checkCollision(coeff, T) && cost + sample_node_->getCost() < goal_->getCost())
    {
      ROS_INFO("A feasiable trajectort to goal was found with cost = %.4f!", cost + sample_node_->getCost());
      goal_->setParent(sample_node_, cost, T, coeff);
      retrieveTraj(goal_);
    }

    // reach_goal_ = checkGoal(sample_node_->getState().pos, goal_->getState().pos);
    visTrajLib();
  }

  // TODO: reach goal and merge informed RRT* to it?
  // retrieveTraj(sample_node_);
  // state_nodes_.sort();

  return true;
}

void KinodynRRTStarPlanner::samplePos()
{
  int collision_flag = 1;
  Eigen::Vector3d sample_pos;
  Eigen::Vector3d near_pos;
  kdres *nearest_node;

  while (collision_flag == 1)
  {
    // sample randomly
    sample_pos = Eigen::Vector3d::Random();
    sample_pos = map_size_ * sample_pos + map_origin_;

    double pos[3] = {sample_pos(0), sample_pos(1), sample_pos(2)};
    nearest_node = kd_nearest(kd_tree, pos);
    near_pos = ((StateNode::Ptr)kd_res_item_data(nearest_node))->getState().pos;

    // ROS_INFO("nearest_pos(%.4f, %.4f, %.4f)\t random_sample(%.4f, %.4f, %.4f)", near_pos(0), near_pos(0), near_pos(0), sample_pos(0), sample_pos(0), sample_pos(0));
    sample_pos = (sample_pos - near_pos).normalized() * step_size_ + near_pos;
    // ROS_INFO("final sampled pos(%.4f, %.4f, %.4f)", sample_pos(0), sample_pos(0), sample_pos(0));

    // check if sample in obscale
    collision_flag = map_->getInflateOccupancy(sample_pos);
  }

  sample_node_ = new StateNode(sample_pos);
}

bool KinodynRRTStarPlanner::checkCollision(Eigen::Matrix<double, 6, 3> coeff, double T)
{
  Eigen::Matrix<double, 6, 1> nature_bais;
  Eigen::Vector3d pos, relative_pos;
  geometry_msgs::Point pt;

  for (double t = 0.0; t < T; t += resolution)
  {
    nature_bais << 1.0, t, t * t, t * t * t,
        t * t * t * t, t * t * t * t * t;
    pos = coeff.transpose() * nature_bais;

    pt.x = pos(0);
    pt.y = pos(1);
    pt.z = pos(2);

    relative_pos = pos - map_origin_;
    if (fabs(relative_pos(0)) > map_size_(0, 0) || fabs(relative_pos(1)) > map_size_(1, 1) || fabs(relative_pos(2)) > map_size_(2, 2) || map_->getInflateOccupancy(pos) == 1)
    {
      return false;
    }
  }
  return true;
}

void KinodynRRTStarPlanner::nearestBackwardNeighbors()
{
  // calc the backward bouding box of state
  double radius = calcSampleRadius();
  double tau = 0.75 * radius;

  DroneState_t state = sample_node_->getState();
  Eigen::Vector3d p = state.pos, v = state.vel, a = state.acc;

  double state_tau[3] = {0.5 * a(0) * tau * tau - v(0) * tau + p(0),
                         0.5 * a(1) * tau * tau - v(1) * tau + p(1),
                         0.5 * a(2) * tau * tau - v(2) * tau + p(2)};
  // find the nodes using kd tree

  kdres *set = kd_nearest_range(kd_tree, state_tau, 4.0 * sqrt(5.0) * tau * tau * tau);

  neighbors_.clear();
  StateNode::Ptr neighbor;
  while (!kd_res_end(set))
  {
    neighbor = (StateNode::Ptr)kd_res_item_data(set);
    if (neighbor != sample_node_)
    {
      neighbors_.push_back(neighbor);
    }

    kd_res_next(set);
  }
  kd_res_free(set);
}

void KinodynRRTStarPlanner::nearestForwardNeighbors()
{
  // calc the forward bounding box of state
  double radius = calcSampleRadius();
  double tau = 0.75 * radius;

  DroneState_t state = sample_node_->getState();
  Eigen::Vector3d p = state.pos, v = state.vel, a = state.acc;

  double state_tau[3] = {0.5 * a(0) * tau * tau + v(0) * tau + p(0),
                         0.5 * a(1) * tau * tau + v(1) * tau + p(1),
                         0.5 * a(2) * tau * tau + v(2) * tau + p(2)};

  // find the nodes using kd tree
  kdres *set = kd_nearest_range(kd_tree, state_tau, 4.0 * sqrt(5.0) * tau * tau * tau);

  neighbors_.clear();
  StateNode::Ptr neighbor;
  while (!kd_res_end(set))
  {
    neighbor = (StateNode::Ptr)kd_res_item_data(set);

    if (neighbor != sample_node_)
    {
      neighbors_.push_back(neighbor);
    }

    kd_res_next(set);
  }
  kd_res_free(set);
}

void KinodynRRTStarPlanner::visTrajLib()
{
  Eigen::Matrix<double, 6, 1> nature_bais;
  Eigen::Matrix<double, 6, 3> coeff;
  double T;
  Eigen::Vector3d pos;
  geometry_msgs::Point pt;
  int marker_id = 0;

  for (auto node = state_nodes_.cbegin(); node != state_nodes_.cend(); ++node, ++marker_id)
  {
    coeff = (*node)->getPolyTraj();
    T = (*node)->getInterval();

    vis_sample_trajectory_.points.clear();
    vis_sample_trajectory_.id = marker_id;
    for (double t = 0.0; t < T; t += 0.1)
    {
      nature_bais << 1.0, t, t * t, t * t * t,
          t * t * t * t, t * t * t * t * t;
      pos = coeff.transpose() * nature_bais;

      pt.x = pos(0);
      pt.y = pos(1);
      pt.z = pos(2);
      vis_sample_trajectory_.points.push_back(pt);
    }

    nature_bais << 1.0, T, T * T, T * T * T,
        T * T * T * T, T * T * T * T * T;

    pos = coeff.transpose() * nature_bais;
    pt.x = pos(0);
    pt.y = pos(1);
    pt.z = pos(2);
    vis_sample_trajectory_.points.push_back(pt);
    vis_sample_points_.points.push_back(pt);
    vis_trajectorylib_.markers.push_back(vis_sample_trajectory_);
  }

  if (goal_->getParent())
  {
    for (auto node = goal_; node->getParent(); node = node->getParent(), ++marker_id)
    {
      coeff = node->getPolyTraj();
      T = node->getInterval();

      vis_final_trajectory_.points.clear();
      vis_final_trajectory_.id = marker_id;
      for (double t = 0.0; t < T; t += 0.1)
      {
        nature_bais << 1.0, t, t * t, t * t * t,
            t * t * t * t, t * t * t * t * t;
        pos = coeff.transpose() * nature_bais;

        pt.x = pos(0);
        pt.y = pos(1);
        pt.z = pos(2);
        vis_final_trajectory_.points.push_back(pt);
      }
      nature_bais << 1.0, T, T * T, T * T * T,
          T * T * T * T, T * T * T * T * T;

      pos = coeff.transpose() * nature_bais;
      pt.x = pos(0);
      pt.y = pos(1);
      pt.z = pos(2);
      vis_final_trajectory_.points.push_back(pt);
      vis_trajectorylib_.markers.push_back(vis_final_trajectory_);
    }
  }

  vis_sample_points_.id = marker_id + 1;
  vis_trajectorylib_.markers.push_back(vis_sample_points_);
  vis_traj_library_pub_.publish(vis_trajectorylib_);
  vis_trajectorylib_.markers.clear();
  vis_sample_points_.points.clear();
};