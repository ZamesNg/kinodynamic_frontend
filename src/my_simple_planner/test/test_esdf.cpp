#include <iostream>
#include <memory>

#include <ros/ros.h>
#include "local_perception/sdf_map.h"

using namespace std;

shared_ptr<SDFMap> sdf_map;

int main(int argc, char **argv)
{
  ros::init(argc,argv,"test_esdf");
  ros::NodeHandle nh("~");

  sdf_map.reset(new SDFMap);
  sdf_map->initMap(nh);

  ros::spin();
  return 0;
}