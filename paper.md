---
title: 'LEIA: Lagrangian / Eulerian Interface Approximation methods'
tags:
  - Level-Set
  - Finite Volume Method
  - Unstructured mesh
  - Computational Fluid Dynamics
  - Two-phase flow
authors:
  - name: Tomislav Maric
    orcid: 0000-0001-8970-1185
    equal-contrib: true
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Julian Reitzel
    orcid: 0000-0002-3787-0283
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 1
affiliations:
 - name: Mathematical Modeling and Analysis, Technical University Darmstadt, Germany
   index: 1
date: 13 August 2017
bibliography: paper.bib

---

# Summary

In 1979 the LS method was developed by Dervieux and Thomasset and was first used
in the CFD context for the advection of two-phase flow. But beyond that, the LS method
is a general concept for describing and advecting surfaces. With several contributions
by Osher and Sethian, the LS method has been popularised and applied in the fields
of image processing, computer vision, computer graphics, computational geometry and
optimization.
Here, the Level Set method is implemented in the OpenFOAM CFD framework for the advection 
of interfaces and the simulation of two-phase flows.


[comment]: # ![Caption for example figure.](level-set-schematic.pdf){ width=50% }

\begin{figure}
    \centering
    \def\svgwidth{0.6\textwidth}
    \input{level-set-schematic.pdf_tex}\notag\\
    \caption{Adapted from Nicoguaro (2018). Level-set method. \href{https://commons.wikimedia.org/wiki/File:Level_set_method.png}{Wikimedia}. Licensed under CC-BY-4.0.}
    \label{fig:level-set}
\end{figure}

# Statement of need

`LEIA` is an OpenFOAM library that implements the level-set method. 
OpenFOAM is a highly efficient C++ framework for CFD, capable of advanced CFD features such as AMR, 
unstructured meshes and MPI parallelism to run on High Performance Computing (HPC) clusters.
The method is not included in the current version of OpenFOAM.
To the authors' knowledge, there is no published implementation for the OpenFOAM framework.
OpenFOAM has the advantage of handling all meshes topologically unstructured.
This raises the importance of this level set implementation, as most existing implementations are based on structured meshes.   

`LEIA` was designed to be used by researchers. It has already been
used in a number of scientific publications [@fricke_locally_2022].

Example scientific publications [@Pearson:2017]

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"





# Acknowledgements



# References
