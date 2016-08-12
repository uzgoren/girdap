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
<div id="flexible" class="panel-collapse collapse noCrossRef">
<div class="panel-body"

> Use your own numerical algorithm or modify one available according to your needs

Many open source and commercial software allows an interface that allow the user to change parameters of a numerical algorithm. Even changing boundary conditions of a partial differential equation dynamically is almost impossible. `girdap` is a library equipped with objects and their methods to allow such a flexibility and allows the user to develop their own algorithms through not so complicated i.e. conditional statements and loops.

At this early stage of the development, girdap will remain as a library while it is envisioned to be equipped with a REPL (read-eval-print loop) type command line interpreter, similar to Matlab, for increasing productivity. 


</div>
</div>


<div>

<figure align="center" style="1px solid #ddd">
<img class="docimage" width="833" height="576" src="{{site.baseurl}}/images/highlight.png" alt="" usemap="#Map" />
<map name="Map">
    <area title="Highly Customizable" href="#flexible" class="noCrossRef accordion-toggle" data-toggle="collapse" data-parent="#accordion" shape="rect" coords="455,57,778,190" />
    <area title="Anisotropic grid refinement" href="front_adapt.html" shape="rect" coords="530,234,833,353" />
    <area title="Easy manage - object oriented" href="front_ooo.html" shape="rect" coords="460,385,790,520" />
    <area title="girdap" title="girdap" href="index.html" shape="rect" coords="0,0,200,180" />
</map>
</figure>

<script src="{{site.baseurl}}/js/jquery.rwdImageMaps.min.js"></script>
<script>
$(document).ready(function(e) {
	$('img[usemap]').rwdImageMaps();
});
</script>


