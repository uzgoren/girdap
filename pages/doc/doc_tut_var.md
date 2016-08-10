---
title: Tutorial - Field Variables
layout: page
---

{% capture md_content %}

Variables are 1D arrays; which maintain an order same as the cells/vertices/faces indices in a grid. In addition to values, variables also contain boundary condition information. 

## Assign values defined by a function
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
  
![Result](res01.gif)

## Communication between grids
  This example generates a volume and a surface grid. Those two are connected through each other by `updateOtherVertex(..)` method. The velocity assigned to volume grid is passed to the surface grid nodes, which are formed after h-refinement. Also an indicator function is generated using the location of the surface.

{% highlight c++ linenos %}  
 Block2* volgrid = new Block2({0,0,0}, {1,1,0}, 50, 50); 

   // Velocity field
   auto uv = volgrid->getVar("u"); auto vv = volgrid->getVar("v"); 
   uv->set(1.0); // set velocity
   vv->set(-0.5); // set velocity
   // New variable at cell center
   volgrid->addVar("f"); auto f = volgrid->getVar("f"); 

   Grid* surf = new Grid(); 

   surf->addVertex( { {0.55,0.32}, {0.58,0.5}, {0.45,0.68}, {0.42,0.46} } ); 
   surf->addCell( { {0,1}, {1,2}, {2,3}, {3,0} } ); 
   // Refine cell; 
   for (auto i=0; i<4; ++i) {
     for (auto c: surf->listCell) if (c->vol().abs() > 0.02) c->adapt[0] = 1;
     surf->adapt(); 
   }
   volgrid->updateOtherVertex(surf);
   // mark location of this surface
   volgrid->indicator(surf, f);

   // Assign velocity variables to surface at vertex  
   surf->addVec("u",1);

   // Get velocity on the surface
   auto us = surf->getVar("u"); auto vs = surf->getVar("v");   
   volgrid->passVar(surf, uv, us); 
   volgrid->passVar(surf, vv, vs);   

   volgrid->writeVTK("vol"); 
   surf->writeVTK("surf"); 

   delete(volgrid); 
   delete(surf); 
{% endhighlight %}
![Result](res02.png)

{% endcapture %}
{{ md_content | markdownify }} 


