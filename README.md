{Computational Framework}:



We have developed a parallel finite element library called HiPerLiFE (High Performance Library for Finite Elements),
which serves as the core numerical engine for the simulations presented in this study.

The installation of the HiPerLiFE libraries can be carried out by following the guidelines provided at: \url{https://hiperlife.gitlab.io/hiperlife/Installation/index.html}.


Code Organization and Project Setup:

For the present work, we implemented the problem-specific code within the HiPerLiFE ecosystem by creating a top-level project folder named {wrinkling\_epithelial\_continuum}. This project folder contains a global {CMakeLists.txt} file, which controls the compilation, and one or more subprojects (applications), each corresponding to a specific simulation setup. For example, the application \texttt{inflate\_hold\_deflate} implements the model of interest here and resides in its own folder with a local \texttt{CMakeLists.txt}.
The code used in this study is openly available at \url{https://github.com/pradeep927/mechanics_of_epithelial_domes}.



The build system relies on CMake to configure and manage the compilation. Two additional files are kept in the project root: {cmake.project.ubuntu.20.04.sh}, a shell script to automate the configuration process on Ubuntu systems; 
and \texttt{userConfig.cmake}, a configuration file where paths to HiPerLiFE, Trilinos, and other dependencies are specified. The typical installation procedure is as follows: Install HiPerLiFE following the instructions at \url{https://hiperlife.gitlab.io/hiperlife/Installation/index.html}. For our simulations, the libraries were installed on a workstation running Ubuntu 22.04.

Place the two files inside the top-level project folder {wrinkling\_epithelial\_continuum}. From the terminal, make the script executable and launch the build configuration:

chmod 777 cmake.project.ubuntu.20.04.sh
./cmake.project.ubuntu.20.04.sh


Enter the build directory and compile:

cd build
make -j4 install


This process generates the executable binary at:

/home/ubuntu/wrinkling_epithelial_continuum/source_compiled/bin/hlinflate_hold_deflate


{Running Simulations}:

To execute a simulation, we create a dedicated folder {run_simulation} that contains the mesh files and the configuration files. Meshes are provided in VTK format and typically correspond to epithelial footprints with prescribed geometries. For example, a mesh named \texttt{aspect\_26} corresponds to a circular domain with aspect ratio 26. Boundary conditions are imposed by carefully fixing the exterior nodes of the domain using elastic springs, which mimic physical anchoring constraints.

The simulation parameters, including references to the mesh files, are specified in a configuration file named {config.cfg}. This file allows the user to set model parameters, time-stepping controls, solver tolerances, and material constants.

A simulation is launched in parallel using MPI as:

mpirun -n 4 /home/ubuntu/wrinkling_epithelial_continuum/source_compiled/bin/hlinflate_hold_deflate config.cfg

Here, the option \texttt{-n 4} specifies the number of processors. This value can be adjusted according to the available computational resources and the problem size.

{Postprocessing and Visualization}:

During execution, the code generates output files in VTK format (e.g., \texttt{sol_dis*.vtk}), which store the displacement and other field variables at different time steps. These files can be directly visualized and post-processed using ParaView. To streamline postprocessing, we provide a ParaView state file (\texttt{paraview\_state.pvsm}), which loads the output files, applies predefined visualization settings, and enables rapid analysis of simulation results, including deformation fields, wrinkling patterns, and stress distributions.

{Code Availiability}:
The code is openly available at \url{https://github.com/pradeep927/mechanics_of_epithelial_domes}.
