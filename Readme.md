# Git Repository for Advanced Macroeconomics II

#### **Instructor:** Sergio Ocampo

#### **Institution:** Western University

#### **Contact:** socampod@uwo.ca

This repository contains the material for the advanced macroeconomics course at Western University.<br/>
This is a second year PhD course focused on computational tools for quantitative Macroeconomics.<br/>
The course syllabus is located in the main folder.<br/>
The lecture handouts are located in the "Slides" folder.<br/>
Accompanying code is located in the "Julia_Code" folder and its subfolders. All code is in Julia and runs in version 1.5.1.<br/>
Problem sets and other assignments are located in the "Problem_Sets" folder.<br/>

## Lecture List

1. Course overview and best practices [Slides](https://github.com/ocamp020/PhD_Macro_Course_Western/blob/master/Slides/PhD_Macro_Comp_1_Handout.pdf)//[Problem Set](https://github.com/ocamp020/PhD_Macro_Course_Western/blob/master/Problem_Sets/Problem_Set_1.pdf)
  - Problem set on dynamic programming with the neoclassical growth model
2. Value function iteration [Slides](https://github.com/ocamp020/PhD_Macro_Course_Western/blob/master/Slides/PhD_Macro_Comp_2_Handout.pdf)//[Problem Set](https://github.com/ocamp020/PhD_Macro_Course_Western/blob/master/Problem_Sets/Problem_Set_2.pdf)
  - Grid search algorithm
  - Howard's policy iteration
  - MacQueen-Porteus bounds
  - Problem set on coding VFI for the neoclassical growth model
3. Interpolation [Slides](https://github.com/ocamp020/PhD_Macro_Course_Western/blob/master/Slides/PhD_Macro_Comp_3_Handout.pdf)//[Problem Set](https://github.com/ocamp020/PhD_Macro_Course_Western/blob/master/Problem_Sets/Problem_Set_3.pdf)
  - Polynomial interpolation
  - Spline interpolation (linear, cubic and shape preserving splines)
  - Curved grids
  - Problem set with applications of interpolation methods
4. Optimization [Slides](https://github.com/ocamp020/PhD_Macro_Course_Western/blob/master/Slides/PhD_Macro_Comp_4_Handout.pdf)//[Problem Set](https://github.com/ocamp020/PhD_Macro_Course_Western/blob/master/Problem_Sets/Problem_Set_4.pdf)
  - Local optimizers
  - Global optimizers
  - Sobol numbers
  - Root finding
  - VFI with continuous choice variables and continuous states
5. Integration [Slides](https://github.com/ocamp020/PhD_Macro_Course_Western/blob/master/Slides/PhD_Macro_Comp_5_Handout.pdf)//[Problem Set](https://github.com/ocamp020/PhD_Macro_Course_Western/blob/master/Problem_Sets/Problem_Set_5.pdf)
  - Monte Carlo Integration (just the concept)
  - Gaussian quadrature methods
  - Discretizing stochastic processes
    * Tauchen (1986)
    * Tauchen & Hussey (1991)
    * Rouwenhorst (1995)
    * Gaussian mixtures (the basics)
    * Extensions (Galindev & Lkhagvasuren, 2010; Fella, Gallipoli & Pan, 2019)
6. The Endogenous Grid Method (EGM) [Slides](https://github.com/ocamp020/PhD_Macro_Course_Western/blob/master/Slides/PhD_Macro_Comp_6_Handout.pdf)//[Problem Set](https://github.com/ocamp020/PhD_Macro_Course_Western/blob/master/Problem_Sets/Problem_Set_6.pdf)
  - Basic EGM: Carroll (2006)
  - EGM with labor supply: Barillas & Fernandez-Villaverde (2007)
  - EGM for non-smooth, non-concave problems: Fella (2014)
  - Envelope condition method (ECM): Maliar & Maliar (2013)
7. Recursive Competitive Equilibria (RCE) [Slides](https://github.com/ocamp020/PhD_Macro_Course_Western/blob/master/Slides/PhD_Macro_Comp_7_Handout.pdf)//[Problem Set](https://github.com/ocamp020/PhD_Macro_Course_Western/blob/master/Problem_Sets/Problem_Set_7.pdf)
  - Decentralized solution to NGM
  - Definition of RCE
  - Applications to economies with distortions (taxes/wedges)
  - Overview of Chari, Kehoe & McGrattan (2007)
  - Overview of Arellano (2008)
8. Heterogeneous Agent Models (RCE) [Slides](https://github.com/ocamp020/PhD_Macro_Course_Western/blob/master/Slides/PhD_Macro_Comp_8_Handout.pdf)//[Problem Set](https://github.com/ocamp020/PhD_Macro_Course_Western/blob/master/Problem_Sets/Problem_Set_8.pdf)
  - Hugget (1993)
  - Aiyagari (1994)
  - OLG - Aiyagari
  - Histogram method of Young (2010)
  - Fast Kernel method of Tan (2020)
9. Computational methods for continuous time models
  - Heterogeneous Agent Continuous Time model
  - Finite Difference Method
    - We will follow Ben Moll's lecture notes [Slides 1](https://benjaminmoll.com/wp-content/uploads/2019/07/Lecture1_Rochester.pdf)//[Slides 2](https://benjaminmoll.com/wp-content/uploads/2019/07/Lecture2_Rochester.pdf)//[Paper](https://benjaminmoll.com/wp-content/uploads/2019/07/HACT.pdf)//[Other resources](https://benjaminmoll.com/lectures/)
  - Markov Chain Approximation Method
    - We will follow Eslami & Phelan (2020) [Slides](https://github.com/ocamp020/PhD_Macro_Course_Western/blob/master/Slides/EslamiPhelan_Slides.pdf)//[Paper](https://github.com/ocamp020/PhD_Macro_Course_Western/blob/master/Slides/EslamiPhelan.pdf)
9. Search models I [Slides](https://github.com/ocamp020/PhD_Macro_Course_Western/blob/master/Slides/PhD_Macro_Comp_9_Handout.pdf)//[Problem Set](https://github.com/ocamp020/PhD_Macro_Course_Western/blob/master/Problem_Sets/Problem_Set_9.pdf)
  - Basic McCall (1970) model
  - Extensions to job separations, on the job search, human capital, savings
