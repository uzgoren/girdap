---
title: "Anisotropic h-refinement and r-refinement"
permalink: /front_adapt.html
tags: [overview]
toc: false
summary: Better accuracy can be achieved by using quality grids, which can be produced by local h-refinement and r-refinement using the right error estimate
---


> It’s possible for me to make a bad movie out of a good script, but I can’t make a good movie from a bad script. *George Clooney*

George Clooney's quote for moview can be altered for numerical simulations:

> It is possible to produce *inaccurate results* out of a *good grid*, but it is almost impossible to produce *accurate results* out of a *bad grid*. > ...


One of the biggest challenge in mesh-based simulation tools remains to be the need of robust and efficient tools for generating and maintaining quality grids around complex geometries for transient multi-scale problems. Adaptive solution techniques allow us to overcome these challenges. These adaptive strategies include

1. **h-adaptation**: local mesh refinement through cell splitting/merging
2. **r-adaptation**: relocating cell vertices





