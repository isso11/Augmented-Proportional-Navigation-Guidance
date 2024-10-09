# Augmented-Proportional-Navigation-Guidance
Nonlinear Augmented Proportional Navigation (APNG) Engagement Simulator. This GUI allows the user to simulate a missile-target engagement using nonlinear APNG guidance law with optional noise and command limits

## Overview

![GUI Preview](https://github.com/isso11/Augmented-Proportional-Navigation-Guidance/blob/main/image.png?raw=true)

The **Nonlinear APNG Engagement Simulator** is a MATLAB-based graphical user interface (GUI) designed to simulate missile-target engagements using the Augmented Proportional Navigation Guidance (APNG) method. This tool allows users to input various parameters related to the missile and target dynamics and visualize the resulting trajectories and guidance commands.

## Features
- **User-Friendly GUI:** An intuitive interface to input initial conditions and simulation parameters.
- **Dynamic Simulations:** Run simulations based on user-defined initial positions, velocities, and other parameters.
- **Visualization:** Plots of engagement trajectories, guidance commands, and other relevant metrics.
- **Customization:** Input fields for noise parameters and guidance command limits.

## Requirements
- MATLAB R2020a or later
- Image Processing Toolbox (for displaying background images)

## Getting Started

### Installation
1. **Clone the Repository:**
   ```bash
   git clone https://github.com/yourusername/AUG_PRO_NAV.git
   ```

2. **Navigate to the Project Directory:**
   ```bash
   cd AUG_PRO_NAV
   ```

3. **Open the MATLAB GUI:**
   - Open `AUG_PRO_NAV.m` in MATLAB.

### Usage
1. **Set the Initial Conditions:**
   - Enter the target and missile initial positions, velocities, heading angles, and other parameters in the input fields.

2. **Run the Simulation:**
   - Click the **Simulate** button to start the simulation.
   - The results will be displayed in the respective plots.

3. **Reset:**
   - Click the **Reset** button to clear the input fields and plots.

4. **Hold on Plots:**
   - Use the **Hold on plots to compare!** button to keep the plots visible for subsequent runs.


## Author
**Islam Elnady**  
Email: islamelnady@yahoo.com  

## REFERENCES: 
 • Ben Dickinson's Toturial Series: Guidance Fundamentals
 • Zarchan, P., 2012. Tactical and Stratigic missile Guidance, AIAA Press,
