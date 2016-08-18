---
title: Grid - Block1
permalink: /doc_obj_grid_block1.html
keywords: girdap, grid, 1D
tags: [grid, objects]
summary: 1D grid blocks
sidebar: mydoc_sidebar
---

## Object: `Block1`

`Block1` inherits all methods of `Grid` while it mainly focuses on generating/manipulating 1D grids formed by `Line` cells.

## Constructor

`Block1()` : generates Grid of Block1 type

`Block1({ coord1 }, { coord2 }, N)` generates a line between _coord1_ and _coord2_ which is made up of N cells; all at a level marked as 0. This means that the grid can not be further coarsened through h-refinement.

`Block1( geo, { {param0}, {param1}, ... }, nx )` generates a predefined geometry which is charaterized by parameters listed next; and using h-refinement approximately _nx_ cells are formed.

_geo_ can have two options:

* "poly": _param\#_ are simply coordinates of connected lines. This can be used to define any polygon with arbitrary number of sides when the first and last parameters point to the same coordinate. If they are different, then the polygon is left unclosed.
* "arc": requires _param0_ to be center; first number in _param1_ to be radius; first number in _param2_ is start angle in degrees and second number in _param2_ is end angle in degrees. If _param2_ is omitted, a full circle will be formed. Degrees can be placed in increasing or decreasing form, to yield counter-clockwise or clockwise trace; which alters the normal direction to point outside or inside, respectively. 

## Methods

`merge( Block1 )` merges one Block1 with another. 