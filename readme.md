# Kinodynamic Frontend

This project focus on quadrotor's `Kinodynamic` planning frontend, including `Kinodynamic-RRT*` and `Kinodynamic-A*`. The word `frontend` means that it will search a coarse(not optimal) trajectory(not path), and the word `kinodynamic` means that the trajectory obeys `kinamatic` constraint(collision free) and `dynamic` constraint(a differential equation).

If this project benefit you in any way, please don't skimp on your star :star2: .

---

## Kinodynamic-RRT*

The project is based on Paper: [Kinodynamic RRT*: Asymptotically optimal motion planning for robots with linear dynamics](https://ieeexplore.ieee.org/abstract/document/6631299/). And some minor changes:

- There are few impelement of `KD tree` can do a range search in each dimension, so this project only search neigbhors with position(3 dimension).

- The paper samples full state (in a quadrotor with minimum jerk performance index, full state means 9 dimension: pos, vel, acc), but this project only samples the position of quadrotor. Accordingly, when calculating trajectory we solve `partial state OBVP` for new sample point and solve `full state OBVP` for rewriting the tree.

With my test, Kinodynamic-RRT* sample the state space randomly, which it is inefficient. It might be horrible in real-time large scale planning.

### quick start

Build the project:
```
catkin_make clean && catkin_make -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -j8
```

Run in terminal:
```
roslaunch my_simple_planner rviz.launch
```

another terminal:
```
roslaunch my_simple_planner test_rrtstar.launch
```

press `G` in rviz and select your goal.

### Example

Search a trajectory in a maze-like `30m x 20m x 3m` environment:

![Kinodynamic-RRT*](image/kinodynamic_rrt_start_30x320x3.gif)


---

## Kinodynamic-A*

Not start yet, but in TODO list.

## Thanks

I use some source file from repo like [mockamap](https://github.com/HKUST-Aerial-Robotics/mockamap), [Teach-Repeat-Replan](https://github.com/HKUST-Aerial-Robotics/Teach-Repeat-Replan), [kdtree](https://github.com/sdeming/kdtree).

Thanks them for their high quality code.
