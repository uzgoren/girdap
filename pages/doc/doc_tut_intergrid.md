---
title: Communication between grids
permalink: doc_tut_intergrid.html
sidebar: mydoc_sidebar
tags: [tutorials, var, objects]
keywords: tutorials, grid, communication, field variables
last_updated: August 10, 2016
folder: doc
toc: false
---

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

{% include image.html file="tut/tut-var02.png" caption="Velocity vectors are transferred from volume grid to surface grid; and indicator function is created purely using surface grid locations." %}



