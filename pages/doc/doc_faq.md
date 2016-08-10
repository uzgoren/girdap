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
  a: <p><i>girdap</i> provides building blocks for numerical simulations of complex transport equations. Tools are built around a grid management platform connected to finite volume method based differential operators. girdap’s primary goals are (1) <b>flexibility</b> and (2) <b>accuracy</b>. </p><p>Flexible tools allow <b>researchers</b> to develop new numerical algorithms and <b>educators</b> to teach students existing algorithms. Tools provided aim to shift programming focus more on physics and numerical method to avoid time consuming programming details.</p><p> The accuracy is handled by automated grid refinement and coarsening based on the solution field. </p><p><i>girdap</i> does not target audiences who would like to get immediate results for an engineering project. It involves a numerical algorithm development phase. Those who would like to skip such a development can refer to a commercial CFD/multiphysics software. </p>
- q: Isn’t this done before? 
  a: <p>Similar projects do exist. <a href='http://www.openfoam.org/'>OpenFOAM</a> is one example; and another one is <a href='http://fenicsproject.org'>FEniCS</a>. Please check them out. They are really good projects. girdap just offers another flavor.</p>
- q: What does girdap mean?
  a: <p><i>girdap</i> is not an acronym. It means whirlpool in Turkish. </p>
- q: What operators are included?
  a: <p>Time derivative, divergence, gradient, Laplacian, and source terms can be defined for a differential equation. </p>
- q: How can I install/use girdap?
  a: <p>You can retrieve the source code (written in c++ following c++11 standards) at github. You can see examples of main_*.cpp files under src directory so that you can compile and use it in the way you like. Use the script provided in the main directory to compile. Note that our focus is mainly on the development of girdap’s skeleton rather than its transportability so you may experience problems while compiling your first application. This is of course going to change with its first release. </p>
---

<p>Check if your question(s) are listed below. Please send me an email or comment on Facebook page about your questions. </p>

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
