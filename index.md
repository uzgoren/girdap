---
title: girdap, Self managed grids for solving PDEs
permalink: index.html
keywords: about
tags: [overview]
sidebar: home_sidebar
summary: girdap is a c++ based object oriented library for multiphysics simulations on self-managed grids 
toc: false
---

<div class="panel-group" id="accordion">

   <!-- Flexible --> 
   <div id="flexible" class="collapse">
       <div class="panel-body">
<h2> Highly Customizable </h2>
          <blockquote>Use your own numerical algorithm or modify one available according to your needs</blockquote>
          <p>Many open source and commercial software allows an interface that allow the user to change parameters of a numerical algorithm. Even changing boundary conditions of a partial differential equation dynamically is almost impossible. `girdap` is a library equipped with objects and their methods to allow such a flexibility and allows the user to develop their own algorithms through not so complicated i.e. conditional statements and loops.</p>
          <p>At this early stage of the development, girdap will remain as a library while it is envisioned to be equipped with a REPL (read-eval-print loop) type command line interpreter, similar to Matlab, for increasing productivity. </p>
      </div>
  </div>

  <!-- accurate --> 
    <div id="adapt" class="panel-collapse collapse">
      <div class="panel-body">
<h2> Anisotropic grid adaptation </h2>

<blockquote> Better accuracy can be achieved by using quality grids, which can be produced by local h-refinement and r-refinement using the right error estimate </blockquote>

<p>Following is George Clooney's quote adapted for numerical simulations: </p>

<pre> It is possible to produce *inaccurate results* out of a *good grid*, but it is almost impossible to produce *accurate results* out of a *bad grid*. </pre>

<p>One of the biggest challenge in mesh-based simulation tools remains to be the need of robust and efficient tools for generating and maintaining quality grids around complex geometries for transient multi-scale problems. Adaptive solution techniques allow us to overcome these challenges. These adaptive strategies include </p>
<ul>
<li> **h-adaptation**: local mesh refinement through cell splitting/merging </li>
<li> **r-adaptation**: relocating cell vertices </li>
</ul>


</div>
</div>



   <!-- ooo --> 
   <div id="ooo" class="collapse">
       <div class="panel-body">
<h2> Object oriented & c++11 standards</h2>
          <blockquote>Use your own numerical algorithm or modify one available according to your needs</blockquote>
          <p>Many open source and commercial software allows an interface that allow the user to change parameters of a numerical algorithm. Even changing boundary conditions of a partial differential equation dynamically is almost impossible. `girdap` is a library equipped with objects and their methods to allow such a flexibility and allows the user to develop their own algorithms through not so complicated i.e. conditional statements and loops.</p>
          <p>At this early stage of the development, girdap will remain as a library while it is envisioned to be equipped with a REPL (read-eval-print loop) type command line interpreter, similar to Matlab, for increasing productivity. </p>
      </div>
  </div>


<figure align="center" style="1px solid #ddd">
<img class="docimage" width="833" height="576" src="{{site.baseurl}}/images/highlight.png" alt="" usemap="#Map">
<map name="Map">
    <area title="Highly Customizable" href="#flexible" class="accordion-toggle" data-toggle="collapse" data-parent="#accordion" shape="rect" coords="455,57,778,190">
    <area nametitle="Anisotropic grid refinement" href="#adapt" class="accordion-toggle" data-toggle="collapse" data-parent="#accordion" shape="rect" coords="530,234,833,353">
    <area title="Easy manage - object oriented" href="#ooo" class="noCrossRef accordion-toggle" data-toggle="collapse" data-parent="#accordion" shape="rect" coords="460,385,790,520">
    <area title="girdap" title="girdap" href="index.html" shape="rect" coords="0,0,200,180">
</map>
</figure>

</div>
<script src="{{site.baseurl}}/js/jquery.rwdImageMaps.min.js"></script>
<script>
var acc; 
$(document).ready(function(e) {
	$('img[usemap]').rwdImageMaps();
	acc = $( "area" ); 

	for (i = 0; i < acc.length; i++) {	
            $( acc[i] ).click(function() {
	      for (j = 0; j < acc.length; j++) {
	      	  if ($( $(acc[j]).attr('href') )[0].classList) {
		  $( $(acc[j]).attr('href') )[0].classList.remove("in");
		  }
	      }
            }); 
        }
});

</script>


