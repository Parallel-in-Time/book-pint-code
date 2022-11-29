# Time Parallel Time Integration : code repository

This repository contains all the source code (Python, Matlab) associated to the book _"Time Parallel Time Integration" (Gander & Lunet 2022)_, along with corrections of the problem given at the end of each chapters.
The folder organization follows the chapter structure of the book, as follow :

- [Chapter 1 : Introduction](#chapter-1---introduction)
- [Chapter 2 : Multiple Shooting Type Methods](#chapter-2---multiple-shooting-type-methods)
- [Chapter 3 : Waveform Relaxation and Domain Decomposition](#chapter-3---waveform-relaxation-and-domain-decomposition)
- [Chapter 4 : Time Multigrid Methods](#chapter-4---time-multigrid-methods)
- [Chapter 5 : Direct Space Time Parallel Solvers](#chapter-5---direct-space-time-parallel-solvers)

---

## Chapter 1 : Introduction

- ### 1.1 : Weather Prediction as an Example

_No associated source code_

- ### 1.2 : Ordinary Differential Equations (ODEs)

_TODO_

- ### 1.3 : Partial Differential Equations (PDEs)

_TODO_

- ### 1.4 : Historical Overview

_TODO_

- ### 1.5 : Problems
    - Problem 1.1 : Study of the Lorenz equations
    - Problem 1.2 : Properties of Forward and Backward Euler
    - Problem 1.3 : Classification of PDEs
    - Problem 1.4 : Taylor expansions for Laplacian
    - Problem 1.5 : Solution of the Heat equation using Fourier's method
    - Problem 1.6 : Numerical stability condition for heat equation
    - Problem 1.7 : Implementation and testing of a heat equation solver
    - Problem 1.8 : Solution of the wave equation
    - Problem 1.9 : Implementation and testing of a wave equation solver
    - Problem 1.10 : Analysis of the Newton's method

---

## Chapter 2 : Multiple Shooting Type Methods

- ### 2.1 : Idea of Nievergelt in 1964

_TODO_

- ### 2.2 : Multiple Shooting Methods in Time

_TODO_

- ### 2.3 : The Parareal Algorithm

_TODO_

- ### 2.4 : Problems
    - Problem 2.1 : Convergence of multiple shooting for linear initial value problems
    - Problem 2.2 : Theoretical speedup of Parareal
    - Problem 2.3 : Parareal for the Lorenz equations
    - Problem 2.4 : Generating functions for Bernoulli number computations
    - Problem 2.5 : Parareal for the Dahlquist equation
    - Problem 2.6 : Convergence bound of Parareal for the heat equation
    - Problem 2.7 : Convergence bound of Parareal for the transport equation
    - Problem 2.8 : Parareal for the heat equation
    - Problem 2.9 : Parareal for the advection-diffusion-reaction equation

---

## Chapter 3 : Waveform Relaxation and Domain Decomposition

- ### 3.1 : Method of Successive Approximations

_TODO_

- ### 3.2 : Classical Waveform Relaxation

_TODO_

- ### 3.3 : Waveform Relaxation Based on Domain Decomposition

_TODO_

- ### 3.4 : Optimized Schwarz Waveform Relaxation

_TODO_

- ### 3.5 : Problems
    - Problem 3.1 : Prove the Gronwall Lemma
    - Problem 3.2 : Prove the Convolution Theorem of Laplace transforms
    - Problem 3.3 : Prove that the Laplace transform of $f'(t)$ is $s\hat{f}(s)-f(0)$.
    - Problem 3.4 : Parallel Schwarz waveform relaxation algorithm for the wave equation
    - Problem 3.5 : Min-max problem for the overlapping optimized Schwarz waveform relaxation
    - Problem 3.6 : Generic Parallel Schwarz waveform relaxation algorithm

---

## Chapter 4 : Time Multigrid Methods

- ### 4.1 : Parabolic Multigrid

_TODO_

- ### 4.2 : Time Multigrid for the Dahlquist Equation

_TODO_

- ### 4.3 : Space Time Multigrid Methods

_TODO_

- ### 4.4 : Block Iteration for iterative Parallel-in-Time methods

_TODO_

- ### 4.4 : Problems
    - Problem 4.1 : To be added

---

## Chapter 5 : Direct Space Time Parallel Solvers

- ### 5.1 : Parallel Predictor Corrector Methods

_TODO_

- ### 5.2 : Boundary Value Methods

_TODO_

- ### 5.3 : Time Parallel Time Stepping

_TODO_

- ### 5.4 : Time Parallel Cyclic Reduction

_TODO_

- ### 5.5 : Time Parallel Methods Based on Laplace Transform

_TODO_

- ### 5.6 : Time Parallelization Based on Diagonalization

_TODO_

- ### 5.7 : Revisionist Integral Deferred Correction (RIDC)

_TODO_

- ### 5.8 : ParaExp

_TODO_

- ### 5.9 : Problems
    - Problem 5.1 : To be added

---