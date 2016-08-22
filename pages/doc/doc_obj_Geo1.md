---
title: Geometry - Geo1
permalink: /doc_obj_geo1.html
keywords: geometry, 1D, objects
tags: [grid, objects]
summary: 1D basic geometric definitions
sidebar: mydoc_sidebar
---

## Object\: Geo1

`Geo1` is a special class equipped with a parametric function to describe simple geometries.

The function accepts parametric variable, `s` (a `double`) to yield a vector in space (a `Vec3`). The start and end values of `s` completes geometric definition.

- `s0` \: Lower bound of the parametric variable (0 by default)
- `s1` \: Upper bound of the parametric variable (1 by default)
- `f`  \: function in terms of lambda (See c++11 tutorials like [Dr. Dobb's](http://www.drdobbs.com/cpp/lambdas-in-c11/240168241))

Proper way to define lambda is as follows:

{% highlight c++ linenos %}
function<Vec3 (double)> f = [=] (double s)->Vec3 {return Vec3(x(s), y(s), z(s));}
{% endhighlight %}

## Constructors

`Geo1()` : Creates a geometric definition using default values; (a line between $$f(0)$$ and $$f(1)$$);

`Geo1(f)` : Creates a geometric definition using a new function defined by the default limits

For example, following defines a heart shape;

{% highlight c++ linenos %}
Geo1 heart( [](double s)->Vec3 {
     		  	   auto t = 4*atan(1.0)/180*s;
			   return Vec3( (1-sin(t))*cos(t), (1-sin(t))*sin(t));
			   } );
{% endhighlight %}

`Geo1(s0, s1, f)` : Creates a geometric definition using a new function and new limits:

{% highlight c++ linenos %}
Geo1 heart(0, 8*atan(1.0), [=](double s)->Vec3 {
     	    	   	return Vec3( (1-sin(t))*cos(t), (1-sin(t))*sin(t));
		    } );
{% endhighlight %}

`heart.f(1.0)` returns a vector with components, $$x=1.2$$ and $$y=3.1$$. 

## Subclasses

Following subclasses consists of the predefined functions (and limits) for simple geometries;

- `Geo1Line(x0, x1)` defines a line between vectors (`Vec3`) $$x_0$$ and $$x_1$$
- `Geo1Circle(x0, r)` defines a full circle (in CCW direction) with a center at $$x_0$$ and a radius of $$r$$
- `Geo1Circle(x0, r, d0, d1)` defines an arc from $$\theta = d_0$$ to $$\theta = d_1$$ (in degrees) with a center at $$x_0$$ and a radius of $$r$$
- `Geo1Sine(x0, x1, a, f)` defines a sinusodial function between x0 and x1, with an amplitude of a and a frequency of f;








