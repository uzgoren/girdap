---
title: "Anisotropic h-refinement and r-refinement"
permalink: /front_adapt.html
tags: [overview]
toc: false
summary: > It’s possible for me to make a bad movie out of a good script, but I can’t make a good movie from a bad script. > George Clooney
---

George Clooney's quote for moview can be altered for numerical simulations:

> It is possible to produce inaccurate results out of a good grid, but it is almost impossible to produce accurate results out of a bad grid. > ...




Managing grid quality and size in grid-based numerical simulations is essential to obtain accurate solutions in a reasonable computation time. According to Fidkowski and Darmofal [4], the present
challenge in mesh-based simulation tools remains to be the need of robust and efficient tools for
generating and maintaining quality grids around complex geometries for transient multi-scale problems.
This is possible through adaptive solution techniques, which include (i) h-adaptation or local mesh
refinement through cell splitting/merging, (ii) r-adaptation or relocating cell vertices, and (iii) p-
adaptation or varying order of numerical approximation. All adaptation methods rely on the choice of
error indicators. Adjoint methods have provided satisfactory results for steady flow problems while
there are limited studies on its use for transient flows. Other examples include residual histograms, flow
geometry and predictable flow features. Active research areas include (i) robust anisotropic grid
adaptation techniques around complex geometries, (ii) identifying error bounds, and (iii) error
estimation for unsteady problems.



