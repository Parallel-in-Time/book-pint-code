# Chapter 5 : Direct Space Time Parallel Solvers

## Section 5.1 : Parallel Predictor Corrector Methods

| 📜 Description | 🧮 Matlab | 🐍 Python |
| :--- | :--- | :--- |
| Predictor-Corrector scheme from Mikanker and Liniger, sequential | [MirakerLinigerS.m](./sec5.1/MirankerLinigerS.m) | [MirakerLinigerS.py](./sec5.1/MirankerLinigerS.py) |
| Predictor-Corrector scheme from Mikanker and Liniger, parallel | [MirakerLinigerP.m](./sec5.1/MirankerLinigerP.m) | [MirakerLinigerP.py](./sec5.1/MirankerLinigerP.py) |
| Example of use on the Lorenz system | [exampleMirankerLiniger.m](./sec5.1/exampleMirankerLiniger.m) | [exampleMirankerLiniger.py](./sec5.1/exampleMirankerLiniger.py) |

## Section 5.3 : Time Parallel Cyclic Reduction

| 📜 Description | 🧮 Matlab | 🐍 Python |
| :--- | :--- | :--- |
| Cyclic Reduction for a bidiagonal system | [CyclicReduction.m](./sec5.3/CyclicReduction.m) | [CyclicReduction.py](./sec5.3/CyclicReduction.py) |
| Cyclic Reduction for a bidiagonal system | [CyclicBackSubstitution.m](./sec5.3/CyclicBackSubstitution.m) | [CyclicBackSubstitution.py](./sec5.3/CyclicBackSubstitution.py) |
| example of use on Dahlquist equation | [exampleCyclicReductionDahlquist.m](./sec5.3/exampleCyclicReductionDahlquist.m) | [exampleCyclicReductionDahlquist.m](./sec5.3/exampleCyclicReductionDahlquist.py) |

> 🔔 The code for `CyclicReduction` given in the book is incomplete, the version above is the correct one. 

## Section 5.4 : Time Parallel Methods Based on Laplace Transform

| 📜 Description | 🧮 Matlab | 🐍 Python |
| :--- | :--- | :--- |
| Time-Parallel Laplace Transform on the Heat equation | [exampleLaplaceTransform.m](sec5.4/exampleLaplaceTransform.m) | [exampleLaplaceTransform.py](sec5.4/exampleLaplaceTransform.py)

## Section 5.5 : Time Parallelization Based on Diagonalization

| 📜 Description | 🧮 Matlab | 🐍 Python |
| :--- | :--- | :--- |
| Paradiag-I with Backward Euler to solve the ODE $y'+ay(t)=f(t)$ | [ODEBEP.m](./sec5.5/ODEBEP.m) | [ODEBEP.py](./sec5.5/ODEBEP.py) |
| Paradiag-II applied on the Dahlquist problem | [exampleDirectDahlquist.m](./sec5.5/exampleDirectDahlquist.m) | [exampleDirectDahlquist.py](./sec5.5/exampleDirectDahlquist.py) |

## Section 5.6 : Time Parallelization Based on Integral Deferred Corrections

| 📜 Description | 🧮 Matlab | 🐍 Python |
| :--- | :--- | :--- |
| Lagrange Weights in IDC for $M \in \{2, 3, 4, 5, 6\}$ | [LagrangeWeights.m](./sec5.6/LagrangeWeights.m) | [LagrangeWeights.py](./sec5.6/LagrangeWeights.py) |
| Integral Deferred Corrections | [IDC.m](./sec5.6/IDC.m) | [IDC.py](./sec5.6/IDC.py) |
| example of IDC on the Dahlquist problem | [exampleIDCDahlquist.m](./sec5.6/exampleIDCDahlquist.m) | [exampleIDCDahlquist.py](./sec5.6/exampleIDCDahlquist.py) |
| Parallel Integral Deferred Corrections | [PIDC.m](./sec5.6/PIDC.m) | [PIDC.py](./sec5.6/PIDC.py) |
| Revisionist Integral Deferred Corrections | [RIDC.m](./sec5.6/RIDC.m) | [RIDC.py](./sec5.6/RIDC.py) |
| Revisionist Integral Deferred Corrections with Restarts | [RIDCRestarts.m](./sec5.6/RIDCRestarts.m) | [RIDCRestarts.py](./sec5.6/RIDCRestarts.py) |

## Section 5.7 : ParaExp

| 📜 Description | 🧮 Matlab | 🐍 Python |
| :--- | :--- | :--- |
| ParaExp used to solve the Heat equation in Matlab | [exampleParaExp.m](./sec5.7/exampleParaExp.m) | [exampleParaExp.py](./sec5.7/exampleParaExp.py) | 

[📖 Table of Content](../README.md)