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

## Section 5.4 : Time Parallel Methods Based on Laplace Transform

| 📜 Description | 🧮 Matlab | 🐍 Python |
| :--- | :--- | :--- |
| Time-Parallel Laplace Transform on the Heat equation | [exampleLaplaceTransform.m](sec5.4/exampleLaplaceTransform.m) | [exampleLaplaceTransform.py](sec5.4/exampleLaplaceTransform.py)

## Section 5.5 : Time Parallelization Based on Diagonalization

| 📜 Description | 🧮 Matlab | 🐍 Python |
| :--- | :--- | :--- |
| Paradiag-I with Backward Euler to solve the ODE $y'+ay(t)=f(t)$ | [ODEBEP.m](./sec5.5/ODEBEP.m) | [ODEBEP.py](./sec5.5/ODEBEP.py) |
| Paradiag-II applied on the Dahlquist problem | [exampleDirectDahlquist.m](./sec5.5/exampleDirectDahlquist.m) | [exampleDirectDahlquist.py](./sec5.5/exampleDirectDahlquist.py) |



[📖 Table of Content](../README.md)