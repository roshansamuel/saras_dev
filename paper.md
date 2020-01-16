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

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

The laws that govern natural systems are often mathematically modelled using a
set of partial differential equations (PDEs).
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
The initial conditions, boundary conditions, and source/forcing terms which appear
commonly in many PDEs, are also defined as derived classes of base ``initial``,
``boundary`` and ``force`` classes.
These classes are written to be readily extensible, so that users can add custom
initial conditions, source terms, and so on.

``Gala`` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for ``Gala`` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. ``Gala`` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the ``Astropy`` package [@astropy] (``astropy.units`` and
``astropy.coordinates``).

``Gala`` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in ``Gala`` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Mathematics

The Navier-Stokes equations, which govern the dynamics of fluid flow, can be written
as
$$ \frac{\partial \mathbf{u}}{\partial t} + \mathbf{u}\cdot\nabla\mathbf{u} = -\nabla p + \nu\nabla^2\mathbf{u} $$
where $\mathbf{u}$ is the velocity field, $p$ is the pressure field, and $\nu$ is
the kinematic viscosity of the fluid.
The fluid is assumed incompressible. Hence $\nabla\cdot\mathbf{u} = 0$, and density is
assumed constant and equal to unity.

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$


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

# References

Example paper.bib file:

@article{Pearson:2017,
  	Adsnote = {Provided by the SAO/NASA Astrophysics Data System},
  	Adsurl = {http://adsabs.harvard.edu/abs/2017arXiv170304627P},
  	Archiveprefix = {arXiv},
  	Author = {{Pearson}, S. and {Price-Whelan}, A.~M. and {Johnston}, K.~V.},
  	Eprint = {1703.04627},
  	Journal = {ArXiv e-prints},
  	Keywords = {Astrophysics - Astrophysics of Galaxies},
  	Month = mar,
  	Title = {{Gaps in Globular Cluster Streams: Pal 5 and the Galactic Bar}},
  	Year = 2017
}

@book{Binney:2008,
  	Adsnote = {Provided by the SAO/NASA Astrophysics Data System},
  	Adsurl = {http://adsabs.harvard.edu/abs/2008gady.book.....B},
  	Author = {{Binney}, J. and {Tremaine}, S.},
  	Booktitle = {Galactic Dynamics: Second Edition, by James Binney and Scott Tremaine.~ISBN 978-0-691-13026-2 (HB).~Published by Princeton University Press, Princeton, NJ USA, 2008.},
  	Publisher = {Princeton University Press},
  	Title = {{Galactic Dynamics: Second Edition}},
  	Year = 2008
}

@article{gaia,
    author = {{Gaia Collaboration}},
    title = "{The Gaia mission}",
    journal = {\aap},
    archivePrefix = "arXiv",
    eprint = {1609.04153},
    primaryClass = "astro-ph.IM",
    keywords = {space vehicles: instruments, Galaxy: structure, astrometry, parallaxes, proper motions, telescopes},
    year = 2016,
    month = nov,
    volume = 595,
    doi = {10.1051/0004-6361/201629272},
    adsurl = {http://adsabs.harvard.edu/abs/2016A%26A...595A...1G},
}

@article{astropy,
    author = {{Astropy Collaboration}},
    title = "{Astropy: A community Python package for astronomy}",
    journal = {\aap},
    archivePrefix = "arXiv",
    eprint = {1307.6212},
    primaryClass = "astro-ph.IM",
    keywords = {methods: data analysis, methods: miscellaneous, virtual observatory tools},
    year = 2013,
    month = oct,
    volume = 558,
    doi = {10.1051/0004-6361/201322068},
    adsurl = {http://adsabs.harvard.edu/abs/2013A%26A...558A..33A}
}

