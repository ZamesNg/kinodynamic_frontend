#ifndef _KINODYNRRTSTAR_H
#define _KINODYNRRTSTAR_H

#include <memory>
#include <vector>
#include <algorithm>
#include <limits>

#include <ros/ros.h>
#include <visualization_msgs/MarkerArray.h>
#include <visualization_msgs/Marker.h>
#include <Eigen/Eigen>

#include "local_perception/sdf_map.h"
#include "kdtree/kdtree.h"

typedef struct DroneState
{
  Eigen::Vector3d pos;
  Eigen::Vector3d vel;
  Eigen::Vector3d acc;
} DroneState_t;

class StateNode
{
public:
  // if there is collision or the trajectory could not be solve, return this value
  static constexpr double kErrorCost = numeric_limits<double>::max();

  StateNode(Eigen::Vector3d pos) : cost_(0.0), heuristic_(0.0), parent_(nullptr), children_({nullptr})
  {
    state_.pos = pos;
    state_.vel = Eigen::Vector3d(0.0, 0.0, 0.0);
    state_.acc = Eigen::Vector3d(0.0, 0.0, 0.0);
  };
  typedef StateNode *Ptr;

  // only calc, but not update
  double calcOptimalTrajWithFullState(DroneState_t parent, double &optimal_T, Eigen::Matrix<double, 6, 3> &coeff) { return calcOptimalTrajWithFullState_(parent, state_, optimal_T, coeff); };
  double calcOptimalTrajWithPartialState(DroneState_t parent, double &optimal_T, Eigen::Matrix<double, 6, 3> &coeff) { return calcOptimalTrajWithPartialState_(parent, state_.pos, optimal_T, coeff); };

  // this will update related data
  bool setParent(Ptr parent, double cost_to_parent, double T, Eigen::Matrix<double, 6, 3> coeff);
  bool setGoal(DroneState_t goal_state);
  void setState(Eigen::Vector3d pos, Eigen::Vector3d vel, Eigen::Vector3d acc)
  {
    state_.pos = pos;
    state_.vel = vel;
    state_.acc = acc;
  };

  Ptr getParent() { return parent_; };
  DroneState_t getState() { return state_; };
  double getCostToParent() { return cost_to_parent_; };
  double getCost() { return cost_; };
  double getHeuristic() { return heuristic_; };
  Eigen::Matrix<double, 6, 3> getPolyTraj() { return polynomial_coeff; };
  double getInterval() { return T_; };

  bool operator<(StateNode node)
  {
    return cost_ + heuristic_ < node.getCost() + node.getHeuristic();
  };

private:
  DroneState_t state_;

  // cost from this node to its parent
  double cost_to_parent_;

  // total cost from begin to this state
  // cost_ = parent->cost_ + cost_to_parent
  double cost_;

  // cost-to-go: the estimation from this state to goals
  double heuristic_;

  // the trajectory for this state to its parent
  // is a polynomial with 5 degree.
  // so we use 6 coeff to represent the trajectory
  Eigen::Matrix<double, 6, 3> polynomial_coeff;
  double T_;
  Ptr parent_;

  // may have many children
  std::vector<Ptr> children_;

  double calcOptimalTrajWithFullState_(DroneState_t start, DroneState_t end, double &optimal_T, Eigen::Matrix<double, 6, 3> &coeff);
  double calcOptimalTrajWithPartialState_(DroneState_t start, Eigen::Vector3d pf, double &optimal_T, Eigen::Matrix<double, 6, 3> &coeff);
  void updateStateFromCoeff();
};

class KinodynRRTStarPlanner
{
public:
  KinodynRRTStarPlanner(){};
  void initPlanner(ros::NodeHandle &nh, shared_ptr<SDFMap> map);
  bool searchTraj(Eigen::Vector3d start_pos, Eigen::Vector3d start_vel, Eigen::Vector3d start_acc,
                  Eigen::Vector3d end_pos);
  std::vector<Eigen::Matrix<double, 6, 3>> getTrajCoeff() { return poly_coeff_; };
  std::vector<double> getTrajInterval() { return intervals_; };
  double getCost() { return total_cost_; };
  void setMap(shared_ptr<SDFMap> map) { map_ = map; };

private:
  bool reach_goal_;
  double goal_tolerance;
  // use to check collision
  double resolution;
  StateNode::Ptr start_;
  StateNode::Ptr goal_;
  StateNode::Ptr sample_node_;
  std::vector<StateNode::Ptr> neighbors_;
  std::list<StateNode::Ptr> state_nodes_;

  std::vector<Eigen::Matrix<double, 6, 3>> poly_coeff_;
  std::vector<double> intervals_;
  double total_cost_;

  // map related param
  // the size = (upper bound - lowwer bound) / 2.0 in each dimension
  Eigen::Matrix3d map_size_;
  Eigen::Vector3d map_origin_;
  shared_ptr<SDFMap> map_;

  kdtree *kd_tree;

  // visualize
  visualization_msgs::MarkerArray vis_trajectorylib_;
  visualization_msgs::Marker vis_trajectory_;

  // ros::Timer vis_timer;
  ros::Publisher vis_traj_library_pub_;

  void samplePos();
  double calcSampleRadius()
  {
    // calc the forward bounding box of state
    int cnt = state_nodes_.size() + 1;
    double log_cnt = log2(cnt) / cnt;
    double ss_volume = map_size_.determinant() * 8.0;

    return 3.3812 * powf64(ss_volume * log_cnt, 1.0 / 6.0);
  };
  void nearestBackwardNeighbors();
  void nearestForwardNeighbors();
  bool checkCollision(Eigen::Matrix<double, 6, 3> coeff, double T);

  bool checkGoal(Eigen::Vector3d pos, Eigen::Vector3d goal)
  {
    double distance = (pos - goal).norm();
    return (distance <= goal_tolerance);
  };

  StateNode::Ptr getCloestNode()
  {
    state_nodes_.sort();
    return (state_nodes_.front());
  };

  void retrieveTraj(StateNode::Ptr node)
  {
    // find the current closest state to goal;
    // StateNode::Ptr node = getCloestNode();
    total_cost_ = node->getCost();
    while (node->getParent() != nullptr)
    {
      poly_coeff_.push_back(node->getPolyTraj());
      intervals_.push_back(node->getInterval());
      node = node->getParent();
    }
    reverse(poly_coeff_.begin(), poly_coeff_.end());
    reverse(intervals_.begin(), intervals_.end());
  };

  void updateNode(StateNode::Ptr node)
  {
    state_nodes_.push_back(node);
    DroneState_t state = node->getState();
    double *pt = new double[3]{state.pos(0), state.pos(1), state.pos(2)};
    kd_insert(kd_tree, pt, node);
  };

  void reBuildKdTree()
  {
    kd_clear(kd_tree);
    DroneState_t state;
    // remember to update all data needed before inqueue
    for (auto node = state_nodes_.cbegin(); node != state_nodes_.cend(); ++node)
    {
      state = (*node)->getState();
      double *pt = new double[3]{state.pos(0), state.pos(1), state.pos(2)};
      kd_insert(kd_tree, pt, *node);
    }
  };

  void visTrajLib()
  {
    Eigen::Matrix<double, 6, 1> nature_bais;
    Eigen::Matrix<double, 6, 3> coeff;
    double T;
    Eigen::Vector3d pos;
    geometry_msgs::Point pt;
    int marker_id = 0;

    for (auto node = state_nodes_.cbegin(); node != state_nodes_.cend(); ++node,++marker_id)
    {
      coeff = (*node)->getPolyTraj();
      T = (*node)->getInterval();

      vis_trajectory_.points.clear();
      vis_trajectory_.id = marker_id;
      for (double t = 0.0; t < T; t += 0.1)
      {
        nature_bais << 1.0, t, t * t, t * t * t,
            t * t * t * t, t * t * t * t * t;
        pos = coeff.transpose() * nature_bais;

        pt.x = pos(0);
        pt.y = pos(1);
        pt.z = pos(2);
        vis_trajectory_.points.push_back(pt);
      }
      vis_trajectorylib_.markers.push_back(vis_trajectory_);
    }
    vis_traj_library_pub_.publish(vis_trajectorylib_);
    vis_trajectorylib_.markers.clear();
  };

  // void visualizeTraj(const ros::TimerEvent &event)
  // {
  //   ROS_INFO("publish trajectory!");
  // }
};

#endif
