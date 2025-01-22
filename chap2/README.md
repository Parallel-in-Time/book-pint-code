# Chapter 2 : Multiple Shooting Type Methods

## Section 2.1 : Idea of Nievergelt in 1964

| 📜 Description | 🧮 Matlab | 🐍 Python |
| :--- | :--- | :--- |
| implementation of the Nievergelt method | [Nievergelt.m](./sec2.1/Nievergelt.m) | [Nievergelt.py](./sec2.1/Nievergelt.py) |
| generic implementation of Forward Euler | [ForwardEuler.m](./sec2.1/ForwardEuler.m), [UForwardEuler.m](./sec2.1/UForwardEuler.m), [exampleODELinear.m](sec2.1/exampleODELinear.m) | [ForwardEuler.py](./sec2.1/ForwardEuler.py), [exampleODELinear.py](sec2.1/exampleODELinear.py) |
| Utility function to find the index of the closest interval | [FindClosestInterval.m](./sec2.1/FindClosestInterval.m) | _in `Nievergelt.py`_ |
| Example code for Nievergelt method with the linear problem | [exampleNievergelt.m](./sec2.1/exampleNievergelt.m) | [exampleNievergelt.py](./sec2.1/exampleNievergelt.py) |
| Script to generate the figure for the linear problem       | [exampleNievergeltLinearFig.m](./sec2.1/exampleNievergeltLinearFig.m) | [exampleNievergeltLinearFig.py](./sec2.1/exampleNievergeltLinearFig.py) |
| Script to generate the figure for the non-linear prolem    | [exampleNievergeltNonLinearFig.m](./sec2.1/exampleNievergeltNonLinearFig.m) | [exampleNievergeltNonLinearFig.py](./sec2.1/exampleNievergeltNonLinearFig.py) |

## Section 2.2 : Multiple Shooting Methods in Time

| 📜 Description | 🧮 Matlab | 🐍 Python |
| :--- | :--- | :--- |
| Multiple Shooting method | [MultipleShooting.m](./sec2.2/MultipleShooting.m), [CForwardEuler.m](./sec2.2/CForwardEuler.m) | [MultipleShooting.py](./sec2.2/MultipleShooting.py) |
| example for Multiple Shooting on the Lorenz system | [exampleMultipleShootingLorenz.m](./sec2.2/exampleMultipleShootingLorenz.m) | [exampleMultipleShootingLorenz.py](./sec2.2/exampleMultipleShootingLorenz.py) |
| scripts to generate the figures | [exampleMultipleShootingLorenzFig.m](./sec2.2/exampleMultipleShootingLorenzFig.m), [exampleMultipleShootingLorenzErrorFig.m](./sec2.2/exampleMultipleShootingLorenzErrorFig.m) | [figures.py](./sec2.2/figures.py) |

## Section 2.3 : The Parareal Algorithm

| 📜 Description | 🧮 Matlab | 🐍 Python |
| :--- | :--- | :--- |
| implementation of Parareal | [Parareal.m](./sec2.3/Parareal.m) | [Parareal.py](./sec2.3/Parareal.py) |
| Parareal on the Lorenz equations | [examplePararealLorenz.m](./sec2.3/examplePararealLorenz.m) | [examplePararealLorenz.py](./sec2.3/examplePararealLorenz.py) |
| error for Parareal on Lorenz equations | [examplePararealLorenzFig.m](./sec2.3/examplePararealLorenzFig.m) | [examplePararealLorenzFig.py](./sec2.3/examplePararealLorenzFig.py) |
| Backward Euler solver for the Dahlquist equation | [DahlquistBE.m](./sec2.3/DahlquistBE.m), [SDahlquistBE.m](./sec2.3/SDahlquistBE.m) | [DahlquistBE.py](./sec2.3/DahlquistBE.py) |
| Parareal on the Dahlquist equation | [examplePararealDahlquist.m](./sec2.3/examplePararealDahlquist.m) | [examplePararealDahlquist.py](./sec2.3/examplePararealDahlquist.py) |
| Backward Euler solver for the Heat equation | [HeatEquationBE.m](./sec2.3/HeatEquationBE.m), [SHeatEquationBE.m](./sec2.3/SHeatEquationBE.m) | [HeatEquationBE.py](./sec2.3/HeatEquationBE.py) |
| Parareal on the Heat equation | [examplePararealHeat.m](./sec2.3/examplePararealHeat.m) | [examplePararealHeat.py](./sec2.3/examplePararealHeat.py) |
| Backward Euler solver for the Transport equation | [TransportBE.m](./sec2.3/TransportBE.m), [STransportBE.m](./sec2.3/STransportBE.m) | [TransportBE.py](./sec2.3/TransportBE.py) |
| Parareal on the Transport equation | [examplePararealTransport.m](./sec2.3/examplePararealTransport.m) | [examplePararealTransport.py](./sec2.3/examplePararealTransport.py) |

## [Section 2.4 : Problems](./sec2.4/README.md)

[📖 Table of Content](../README.md)