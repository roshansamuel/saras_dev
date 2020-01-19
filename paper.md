---
title: 'Saras: A General-Purpose PDE Solver C++ for Fluid Dynamics'

tags:
  - C++
  - PDE
  - turbulence
  - fluid dynamics

authors:
  - name: Mahendra K. Verma
    orcid: 0000-0002-3380-4561
    affiliation: 2
  - name: Roshan Samuel
    orcid: 0000-0002-1280-9881
    affiliation: 1
  - name: Shashwat Bhattacharya
    orcid: 0000-0001-7462-7680
    affiliation: 1
  - name: Ali Asad
    orcid: 0000-0001-9704-6686
    affiliation: 2
  - name: Soumyadeep Chatterjee
    orcid: 0000-0001-7957-1727
    affiliation: 2

affiliations:
 - name: Department of Mechanical Engineering, Indian Institute of Technology - Kanpur
   index: 1
 - name: Department of Physics, Indian Institute of Technology - Kanpur
   index: 2

date: 15 January 2020

bibliography: paper.bib

---

# Summary

The laws that govern natural systems can often be modelled mathematically using
partial differential equations (PDEs).
More often than not, the resultant PDEs are not amenable to analytical tools for
solving them, and numerical techniques remain the only recourse in such cases.
As a result, efficiently solving PDEs numerically is a primary step in understanding
the physics of various systems.
``Saras`` is a general-purpose PDE solver written in objective C++.
In ``Saras``, the underlying mathematical constructs used to define a PDE, like
vector and scalar fields, are defined as classes.
Moreover, vector calculus operations associated with such fields, like gradient and
divergence, are defined as functions of the classes.

This design makes the code intuitive, allowing users to quickly cast PDEs into
readable code.
The initial conditions, boundary conditions, and source/forcing terms, which appear
commonly in many PDEs, are also defined as derived classes of base ``initial``,
``boundary`` and ``force`` classes.
These classes are written to be readily extensible, so that users can add custom
initial conditions, source terms, and so on.

``Saras`` also includes solvers for hydrodynamic flow, namely the incompressible
Navier-Stokes initial value problem (IVP), as well as for scalar convection,
like Rayleigh Benard Convection.
Presently, we use semi-implicit Crank-Nicholson `[@Crank:1947]` for time-advancing
the IVP.
The solver uses Marker and Cell method (MAC) `[@Harlow:PF1965]` for discretizing
the velocity and pressure field.

# Mathematics

The Navier-Stokes equations, which govern the dynamics of fluid flow, can be written as
$$ \frac{\partial \mathbf{u}}{\partial t} + \mathbf{u}\cdot\nabla\mathbf{u} = -\nabla p + \mathbf{f} + \nu\nabla^2\mathbf{u}, $$
where $\mathbf{u}$ is the velocity field, $p$ is the pressure field, $\mathbf{f}$ is
the forcing term, and $\nu$ is the kinematic viscosity of the fluid.
The fluid is assumed incompressible. Hence $\nabla\cdot\mathbf{u} = 0$, and density is
assumed constant and equal to unity.

If the velocity and pressure field at time $t = n$ are denoted as $\mathbf{u}^n$ and $p^n$
respectively, then the corresponding fields at the next time-step, $t = n+1$, namely
$\mathbf{u}^{n+1}$ and $p^{n+1}$, can be calculated as described next.
Initially, an intermediate velocity field is calculated using the known values,
$\mathbf{u}^n$ and $p^n$, as
$$\mathbf{u}^* = \mathbf{u}_{n} + \Delta t\left[\nu\nabla^2 \left( \frac{\mathbf{u}_n + \mathbf{u}^*}{2}\right) - \mathbf{u}_n.\nabla\mathbf{u}_n - \nabla p_n\right]. $$
The forcing term has been neglected for simplicity.
Note that the diffusion term (also called the viscous term) has been split into two, with half
the contribution from the velocity field at $t = n$, and the other half attributed to the
guessed velocity field, $\mathbf{u}^*$.
This results in the following implicit equation,
$$\mathbf{u}^* - \Delta t\left[\frac{\nu\nabla^2\mathbf{u}^*}{2} \right ] = \mathbf{u}_{n} + \Delta t\left[\frac{\nu\nabla^2\mathbf{u}_n}{2} - \mathbf{u}_n.\nabla\mathbf{u}_n - \nabla p_n\right]. $$
The above equation has to be solved iteratively, and this is achieved through
OpenMP parallelized Jacobi iterations.
The intermediate velocity field, $\mathbf{u}^*$, will not satisfy the continuity equation,
and requires to be corrected appropriately.
This correction is obtained from the pressure correction term, which is in turn computed from
the pressure Poisson equation,
$$\nabla^2 p^* = \frac{\nabla.\mathbf{u}^*}{\Delta t}. $$
``Saras`` uses a Geometric Multigrid library to solve the above equation [@Wesseling:MG2004, @Briggs:MG2000].
Presently the library offers the Full Multigrid (FMG) V-Cycle to solve the Poisson equation.
Other methods like F-Cycle and W-Cycle are planned updates to the library in future.



# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this: ![Example figure.](figure.png)

# Acknowledgements

We acknowledge contributions from Gaurav Gautham, Saurav Bhattacharjee, and Rishabh Sahu,
and support from Prof Fahad Anwer during the genesis of this project.

---

# References

Example paper.bib file:

@Article{Crank:1947,
    author={Crank, J. and Nicolson, P.},
    title={A practical method for numerical evaluation of solutions of partial differential equations of the heat-conduction type},
    journal={Proc. Camb. Phil. Soc.},
    year={1947},
    month={Jan},
    volume={43},
    number={1},
    pages={50--67},
    doi={10.1017/S0305004100023197},
}

@article{Harlow:PF1965,
    author = {Harlow,Francis H. and Welch,J. Eddie},
    title = {Numerical Calculation of Time‐Dependent Viscous Incompressible Flow of Fluid with Free Surface},
    journal = {Phys. Fluids},
    volume = {8},
    number = {12},
    pages = {2182-2189},
    year = {1965}
}

@book{Wesseling:MG2004,
    title={An Introduction to Multigrid Methods},
    author={Wesseling, P.},
    isbn={9781930217089},
    lccn={91024430},
    series={An Introduction to Multigrid Methods},
    year={2004},
    publisher={R.T. Edwards}
}

@book{Briggs:MG2000,
    title={A Multigrid Tutorial: Second Edition},
    author={Briggs, W.L. and Henson, V.E. and McCormick, S.F.},
    isbn={9780898714623},
    lccn={87062333},
    series={Other Titles in Applied Mathematics},
    year={2000},
    publisher={Society for Industrial and Applied Mathematics}
}
