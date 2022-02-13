# SPH-sandbox
Interactive Smoothed Particle Hydrodynamic sandbox simulating a simple 2D flow arround a circular obstacle.

## Project content

The repository includes two python scripts:
- GridSPH.py: This script implements the physics class that manages the SPH model. This "home-made" model defines the main flow parameters, the type of fluid particles, the interactions (particle-particle and particle-walls), the obstacles and the regions (inlet, outlet, actual scene). The GridSPH class is imported in the next script and used in the interactive interface.
- mainInteractive.py: This script executes an interactive interface to run the simulation and propose a few settings for displaying of the scene. 
 
![view](https://user-images.githubusercontent.com/57091666/153754191-1d3fea5f-755b-4d12-b685-0ab7fddd4a88.png)

## How to use it?

This project can be executed directly after pulling the project locally. We simply need to execute the script mainInteractive.py in a shell. Running the script in Spyder or PyCharm (not tested) might cause a crash of the python kernel when you will try to close the window. This bug does not occcur when you run the script in command-line.

