---
title: Field Variables
permalink: doc_tut_var.html
sidebar: mydoc_sidebar
tags: [tutorials, var, objects]
keywords: tutorials, grid, variables, field variables
last_updated: August 10, 2016
folder: doc
toc: false
---

Variables are 1D arrays; which maintain an order same as the cells/vertices/faces indices in a grid. In addition to values, variables also contain boundary condition information. 

This example creates a 5x5 block using quads in a unit 2D domain. A new variable is defined and named as `f`; where boundary conditions on `east` and `north` are defined using Neumann and Dirichlet conditions. For `south` and `west`, default boundary condition, i.e. zero gradient, is used.

Typical loop over celllist is given in lines xx-xx, where `c->id` gives the array's index.

Grid adaptation based on gradient of the function is used to refine base grid. 

{% highlight c++ linenos %}
 double pi = 4*atan(1.0); 

 Block2* volgrid = new Block2({0,0,0}, {1,1,0}, 5, 5); 
  
 // add a new variable
 volgrid->addVar("f"); 
 auto f = volgrid->getVar("f"); // variable handle
  
 f->setBC("east", "grad", 0);   // This is the default
 f->setBC("north", "val", 1);   //   

 for (auto i=0; i < 4; ++i) { 
   for (auto c : volgrid->listCell) {
     auto x = c->getCoord(); // cell-centers
     f->set(c->id, sin(3*pi*x[0])*cos(2*pi*x[1]));
   }
   volgrid->solBasedAdapt2(volgrid->getError(f)); 
   volgrid->adapt(); 
   volgrid->writeVTK("field_"); 
 }
  
 delete(volgrid); 

{% endhighlight %}

Above code produces the following result:
  
{% include image.html file="tut/tut-var01.gif" caption="Field variable assigned by function and refined based on error. Boundary conditions are reflected to the solution field." %}



