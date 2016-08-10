---
title: Frequently asked questions
permalink: doc_faq.html
sidebar: mydoc_sidebar
tags: [overview]
keywords: frequently asked questions, FAQ, question and answer
last_updated: August 10, 2016
toc: false
folder: doc

dfaq:
- q: What is girdap?
  a: _girdap_ provides building blocks for numerical simulations of complex transport equations. Tools are built around a grid management platform connected to finite volume method based differential operators. girdap’s primary goals are (1) _flexibility_ and (2) _accuracy_. <p>Flexible tools allow **researchers** to develop new numerical algorithms and **educators** to teach students existing algorithms. Tools provided aim to shift programming focus more on physics and numerical method to avoid time consuming programming details.<p> The accuracy is handled by automated grid refinement and coarsening based on the solution field. <p>_girdap_ does not target audiences who would like to get immediate results for an engineering project. It involves a numerical algorithm development phase. Those who would like to skip such a development can refer to a commercial CFD/multiphysics software. 
- q: Isn’t this done before? 
  a: Similar projects do exist. [OpenFOAM](http://www.openfoam.org/) is one example; and another one is [FEniCS](http://fenicsproject.org). Please check them out. They are really good projects. girdap just offers another flavor.
- q: What does girdap mean?
  a: _girdap_ is not an acronym. It means whirlpool in Turkish. 
- q: What operators are included?
  a: Time derivative, divergence, gradient, Laplacian, and source terms can be defined for a differential equation. 
- q: How can I install/use girdap?
  a: You can retrieve the source code (written in c++ following c++11 standards) at github. You can see examples of main_*.cpp files under src directory so that you can compile and use it in the way you like. Use the script provided in the main directory to compile. Note that our focus is mainly on the development of girdap’s skeleton rather than its transportability so you may experience problems while compiling your first application. This is of course going to change with its first release. 
---

<p>If you want to use an FAQ format, use the syntax shown on the faq.html page. Rather than including code samples here (which are bulky with a lot of nested <code>div</code> tags), just look at the source in the mydoc_faq.html theme file.</p>

<div class="panel-group" id="accordion">
{% for q in page.dfaq %}
                    <div class="panel panel-default">
                        <div class="panel-heading">
                            <h4 class="panel-title">
                                <a class="noCrossRef accordion-toggle" data-toggle="collapse" data-parent="#accordion" href="#collapse{{ forloop.index }}">
				{{ q.q }} 
				</a>
                            </h4>
                        </div>
                        <div id="collapse{{ forloop.index }}" class="panel-collapse collapse noCrossRef">
                            <div class="panel-body">
                            {{ q.a }} 
                            </div>
                        </div>
                    </div>
                    <!-- /.panel -->
{% endfor %}
</div>
<!-- /.panel-group -->

{% include links.html %}
