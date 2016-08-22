---
title: Grid - Introduction
permalink: /doc_obj_grid.html
keywords: girdap, grid, 1D
tags: [grid, objects]
summary: Introducing surface and volume grids
sidebar: mydoc_sidebar
---

## Object: `Grid`

- `Grid` object mainly consists of two arrays of `shared_ptrs`. One is for `Vertex` for position vectors and the other for `Cell` for connectivity information. Grid does not impose a particular cell structure; it can consist of different cells. 

## Constructor

- `grid()` : forms blank grid. It can be used for variable decleration. The constructor however sets default minimum and maximum levels for h-refinement; and cfl number.

- `grid( { coord } )` : ( coord ->  x, y, z ) forms a grid with a single vertex.

- `grid( { {coord0}, {coord1}, ...} )` : forms a grid with many vertices. Connectivity information is not yet defined.

- `grid( { {coord0}, {coord1}, ...}, { {cell0}, {cell2}, ... } } )` : forms a grid with many vertices and cells. _cell0_, _cell1_ are integers pointing out the order of vertex in the first list, starting from 0.

## Methods


