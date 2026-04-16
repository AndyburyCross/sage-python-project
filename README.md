# Sage + Scipy Python Project (No Sage Install Required)
This project demonstrates how to run a Python project using:

- SageMath
- NumPy
- SciPy
- Matplotlib

without requiring users to install SageMath, Python, Conda or pip.

Everything runs inside Docker

## Requirements

- Docker Desktop
https://www.docker.com/products/docker-desktop

## Run the project

Open a terminal (eg Powershell) and run:

git clone https://github.com/AndyBuryCross/sage-python-project.git
cd sage-python-project
docker build -t sage-project
docker run -v {$PWD}:/app sage-project

The script will:
- Perform a symbolic computation using Sage
- Solve a numerical ODE using SciPy
- Save a plot as output.png in the project folder