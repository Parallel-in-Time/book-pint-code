# Time Parallel Time Integration : code repository

This repository contains all the source code (Python, Matlab) associated to the book 
_["Time Parallel Time Integration" (Gander & Lunet 2024)](https://epubs.siam.org/doi/book/10.1137/1.9781611978025)_, 
along with corrections of the problems given at the end of each chapters.
The folder organization follows the chapter structure of the book, as below.

## [Chapter 1 : Introduction](./chap1/README.md)

- 1.1 : Weather Prediction as an Example
- [1.2 : Ordinary Differential Equations (ODEs)](./chap1/README.md#section-12--ordinary-differential-equations-odes)
- [1.3 : Partial Differential Equations (PDEs)](./chap1/README.md#section-13--partial-differential-equations-pdes)
- [1.4 : Historical Overview](./chap1/README.md#section-14--historical-overview)
- [1.5 : Problems](./chap1/sec1.5/README.md)

## [Chapter 2 : Multiple Shooting Type Methods](./chap2/README.md)

- [2.1 : Idea of Nievergelt in 1964](./chap2/README.md#section-21--idea-of-nievergelt-in-1964)
- [2.2 : Multiple Shooting Methods in Time](./chap2/README.md#section-22--multiple-shooting-methods-in-time)
- [2.3 : The Parareal Algorithm](./chap2/README.md#section-23--the-parareal-algorithm)
- [2.4 : Problems](./chap2/sec2.4/README.md)

## [Chapter 3 : Waveform Relaxation and Domain Decomposition](./chap3/README.md)

- 3.1 : Method of Successive Approximations
- 3.2 : Classical Waveform Relaxation
- [3.3 : Waveform Relaxation Based on Domain Decomposition](./chap3/README.md#section-33--waveform-relaxation-based-on-domain-decomposition)
- 3.4 : Optimized Schwarz Waveform Relaxation
- 3.5 : Problems

## Chapter 4 : Time Multigrid Methods

- 4.1 : Time Parallel Time Stepping
- [4.2 : Parabolic Multigrid](./chap4/README.md#section-42--parabolic-multigrid)
- [4.3 : Time Multigrid for the Dahlquist Equation](./chap4/README.md#section-43--time-multigrid-for-the-dahlquist-equation)
- [4.4 : Space Time Multigrid Methods](./chap4/sec4.4/README.md)
- [4.5 : MultiGrid Reduction In Time](./chap4/sec4.5/README.md)
- [4.6 : Block Iteration and Generating Functions](./chap4/sec4.6/README.md)
- [4.7 : Parallel Full Approximation Scheme in Space-Time](./chap4/sec4.7/README.md)
- [4.8 : Problems](./chap4/sec4.5/README.md)

## Chapter 5 : Direct Space Time Parallel Solvers

- [5.1 : Parallel Predictor Corrector Methods](./chap5/sec5.1/README.md)
- [5.2 : Boundary Value Methods](./chap5/sec5.2/README.md)
- [5.3 : Time Parallel Cyclic Reduction](./chap5/sec5.3/README.md)
- [5.4 : Time Parallel Methods Based on Laplace Transform](./chap5/sec5.4/README.md)
- [5.5 : Time Parallelization Based on Diagonalization](./chap5/sec5.5/README.md)
- [5.6 : Time Parallelization Based on Integral Deferred Corrections](./chap5/sec5.6/README.md)
- [5.7 : ParaExp](./chap5/sec5.7/README.md)
- [5.8 : Problems](./chap5/sec5.8/README.md)

## Cite this book

```bibtex
@book{gander2024time,
	address = {Philadelphia, PA},
	author = {Gander, Martin J. and Lunet, Thibaut},
	doi = {10.1137/1.9781611978025},
	edition = {},
	eprint = {https://epubs.siam.org/doi/pdf/10.1137/1.9781611978025},
	publisher = {Society for Industrial and Applied Mathematics},
	title = {Time Parallel Time Integration},
	url = {https://epubs.siam.org/doi/abs/10.1137/1.9781611978025},
	year = {2024}
}
```