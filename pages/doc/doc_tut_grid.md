---
title: Tutorial - Grids
permalink: doc_tut_grid.html
sidebar: mydoc_sidebar
tags: [tutorials, grid]
keywords: tutorials, grid generation, basics
last_updated: August 10, 2016
folder: doc
toc: false
---

## Your first grid
This example creates a grid handle and manually adds vertices and cells. h-refinement can be applied through adapt flag on the first cell. The output is written in VTK format to be visualized in additional software, i.e. Paraview.

{% highlight c++ linenos %}
Grid* grid = new Grid();
grid->addVertex({ {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0} }); 

grid->addCell( {0, 1, 2, 3} ) ; 

for (auto i =0; i<3; ++i) {
   grid->listCell[0]->adapt = {1, 1}; 
   grid->adapt(); 
   grid->writeVTK("myFirstGrid_"); 
}
  
delete(grid); 
{% endhighlight %}

Above code produces the following result:
  
{% include image.html file="tut/tut-grid01.png" alt="My first grid" caption="My first grid" %}




