# Perception-Aware Perching on Powerlines with Multirotors

This repo contains the code for the NLP-based perching trajectory generation presented in our paper [Paneque et al. RA-L'22](https://rpg.ifi.uzh.ch/docs/RAL22_Paneque.pdf).

[![Perception-Aware Perching on Powerlines with Multirotors (Narrated Video)](media/youtube_thumbnail.jpg)](https://youtu.be/JsPavnsfpbk)

In addition, we also provide the following example trajectories:
 
Upside down             | 90º             |  90º (Perception Aware)
:-------------------------:|:-------------------------:|:-------------------------:
![Upside Down](media/upside_down.gif)  |  ![90º (No Perception Awareness)](media/90_deg_no_pa.gif)  |  ![90º (With Perception Awareness)](media/90_deg_pa.gif)

## Citing 
If you use this code in an academic context, please cite the following publication:

Paneque, J. L., Martinez-de-Dios, J. R., Ollero, A., Hanover, D., Sun, S., Romero, A., & Scaramuzza, D. (2022). **Perception-Aware Perching on Powerlines with Multirotors**. IEEE Robotics and Automation Letters ([PDF](https://rpg.ifi.uzh.ch/docs/RAL22_Paneque.pdf))

```
@article{paneque2022perching,
  author={Paneque, Julio L. and Dios, Jose Ramiro Martínez-de and Ollero, Anibal and Hanover, Drew and Sun, Sihao and Romero, Angel and Scaramuzza, Davide},
  journal={IEEE Robotics and Automation Letters}, 
  title={Perception-Aware Perching on Powerlines With Multirotors}, 
  year={2022},
  volume={7},
  number={2},
  pages={3077-3084},
  doi={10.1109/LRA.2022.3145514}}
```

## License 
MIT License. Copyright (C) 2022 Julio L. Paneque, Jose Ramiro Martínez de Dios and Anibal Ollero (GRVC Robotics Lab, Universidad de Sevilla) and Drew Hanover, Sihao Sun, Ángel Romero and Davide Scaramuzza (Robotics and Perception Group, University of Zurich).

This is research code, expect that it changes often and any fitness for a particular purpose is disclaimed.

This work has dependencies on other libraries which are individually cited when appearing.

## Requirements 
The code was tested with Ubuntu 18.04, Matlab R2020b, FORCESPRO 5.1.0, Casadi 3.5.1, g++ 10.3.0 and libYAML-cpp 0.5. Compatibility with other versions should be possible but is not tested.

## Installation and usage

### Generating the solver (`solver_generation`)
The NLP solver is generated and run using the FORCESPRO software (https://forces.embotech.com/). A free academic license of the software can be requested in [this link](https://my.embotech.com/auth/sign_up).

Once you have the license, install FORCESPRO following [these instructions](https://forces.embotech.com/Documentation/installation/obtaining.html), and then install the MATLAB client following [these instructions](https://forces.embotech.com/Documentation/installation/matlab.html). After that, navigate to the `solver_generation` folder and run:

```
generate_solver.m
```

This will generate all the required files and copy them to the `solver_interface` folder.

### Using the solver (`solver_interface`)

After the solver is generated, it can be used from the C++ interface to compute the required perching trajectories. To compile the interface, navigate to `solver_interface` and run:

```
g++ -o compute_perching_traj -Iinclude -std=c++17 src/perch_recovery_planner.cpp src/parameters/quad.cpp src/parameters/lines.cpp src/parameters/costs.cpp src/parameters/reference.cpp src/parameters/xinit.cpp extern/solver/PerchingSolver/lib/libPerchingSolver.so extern/solver/PerchingSolver_adtool2forces.c extern/solver/PerchingSolver_casadi.c -lm -lyaml-cpp
```
Once the code is compiled, it can be run given a path to a `YAML` file with the desired problem. This file provides the paths for the following sub-files (also in `YAML` format):

* The quadrotor parameters (`quad_config`)
* The powerlines setup (`lines_config`)
* The NLP costs for perching (`perching_costs`)
* The NLP costs for recovery (`recovery_costs`)
* The perching reference (`reference`)
* The starting position (`xinit`)

We provide the configuration files for three different perching maneuvers: **upside_down**, **90_deg_no_pa** and **90_deg_pa**. To run the first one, simply run: 

```
./compute_perching_traj config/upside_down.yaml
```
This will compute the maneuver and generate an `upside_down.csv` file in the current directory. We already provide the resulting `.csv` files of the three maneuvers so they can be seen without compiling any code. 

### Visualizing the results (`trajectory_visualization`)

We provide a optional python script to visualize the generated maneuvers in 3D animated plots. The required packages to run it are listed in `trajectory_visualization/requirements.txt`. The script takes as inputs the trajectory file, the quadrotor configuration file and the powerlines setup file. Additional arguments `--show_vel` and `--save` can be used to show the quadrotor velocity overtime and to save the resulting animation as a `gif` file (note: saving the animation will take some time before anything is shown).

You can play an example animation by running:
```
python visualize_trajectory.py -t ../solver_interface/upside_down.csv -q ../solver_interface/config/quad_parameters.yaml -l ../solver_interface/config/upside_down/line_parameters.yaml --show_vel
```

This will show an animation of the **upside_down** perching maneuver.
