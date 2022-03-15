#include <iostream>
#include <memory>

#include <ros/ros.h>
#include "local_perception/sdf_map.h"

using namespace std;

shared_ptr<SDFMap> sdf_map;

void visualizeESDFGrid(const ros::TimerEvent& event);

int main(int argc, char **argv)
{
  ros::init(argc,argv,"test_esdf");
  ros::NodeHandle nh("~");

  sdf_map.reset(new SDFMap);
  sdf_map->initMap(nh);

  ros::Timer vis_esdf_timer = nh.createTimer(ros::Duration(1.0),visualizeESDFGrid);
  ros::spin();
  return 0;
}

void visualizeESDFGrid(const ros::TimerEvent& event)
{   
    if( !sdf_map->hasDepthObservation() ) 
        return;

    sdf_map->publishMap();
    sdf_map->publishMapInflate();
    sdf_map->publishUpdateRange();
    sdf_map->publishESDF();
}