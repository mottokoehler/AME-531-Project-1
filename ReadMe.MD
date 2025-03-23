# AME 531 Project 1

## Overview
This project solves the Navier-Stokes equations using different numerical methods (explicit, implicit, and Runge-Kutta) and plots the results using `matplotlibcpp`.

## Dependencies
- C++17 compatible compiler
- Visual Studio 2019 or later
- `matplotlibcpp` for plotting
- `matplotlib` and `Python` for `matplotlibcpp`
- ****Python 3.10 or lower****

## Installation

### Install `matplotlib` and `Python`
1. Install Python from [python.org](https://www.python.org/).
    Ensure version 3.8 is installed (only version tested with, but 3.11+ will definitely not work)
2. Install `matplotlib` using pip:  
   pip install matplotlib
3. Install 'numpy' using pip:
   pip install numpy

### Install `matplotlibcpp`
1. Clone the `matplotlibcpp` repository:  https://github.com/lava/matplotlib-cpp
2. Copy the `matplotlibcpp.h` file to your project directory.

## Project Setup
1. Open the project in Visual Studio.
2. Ensure that the include directories and library directories are set up correctly:
   - Right-click on the project in Solution Explorer and select `Properties`.
   - Go to `VC++ Directories > Include Directories>` and add the path to `matplotlib` and 'eigen-3.4.0' subdirectories.
    - \Python38\Lib\site-packages\numpy\core\include
    - \Python38\include
    - \source\repos\mottokoehler\AME-531-Project-1\AME 531 Project 1\matplotlib
    - \source\repos\mottokoehler\AME-531-Project-1\AME 531 Project 1\eigen-3.4.0
   - Go to `VC++ Directories > Library Directories>` and add the path to the Python libraries.
    - \Python38\libs

## Build Instructions
1. Open the project in Visual Studio.
2. Build the project by selecting `Build > Build Solution` or pressing `Ctrl+Shift+B`.

## Running the Project
1. Place the `input.k` file in the current working directory.
2. Run the executable from Visual Studio or the command line.

## Sample `input.k` File


## Contact
For any questions or issues, please contact Michael Koehler at mokoehler@arizona.edu.
      