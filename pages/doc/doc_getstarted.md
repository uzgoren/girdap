---
title: Getting started
keywords: about
tags: [overview]
sidebar: mydoc_sidebar
permalink: doc_getstarted.html
summary: girdap is a c++ based object oriented library for multiphysics simulations on self-managed grids.
series: "ACME series"
weight: 1.0
---

## Prepare your system

In order to compile and run girdap, it is recommended to prepare your system first and install necessary compilers and build tools for Windows, Ubuntu and Mac OS. Details can be found [here](doc_gettools.html). 

## Download or clone _girdap_

First download or clone _girdap_ from the [Github repo](https://github.com/uzgoren/girdap).

{% include note.html content="It is recommended to download `develop` branch rather than the `master` branch as many described functionalities are available in `develop` branch at this point." %}

## Build _girdap_

{% include warning.html content="Requires `cmake`" %}

Extract the package into a directory. The girdap's base folder is named as `girdap` as default. It comes with the following directory structure;
<pre>
- girdap
  + src      // *.cpp files are here; also basic main.cpp is placed here; 
  + include  // *.hpp files are here; 
  + bin      // executables are placed here after build
  + library  // a library w/o main.cpp is placed here after build
  + example  // various examples of main.cpp can be found here; 
</pre>

`example` direction contain `main.cpp` files which can utilize _girdap_'s functionality. Some of the tutorials will be placed inside this directory. Develop your own or modify one of them as your `girdap` code and place it in the root of src directory and follow the procedure below to build your code:

{% highlight bash linenos %}
cd dir_of_your_choice
tar -xzvf girdap.tar.gz
cd girdap
# Make changes to the main_xxx.cpp
# modify CMakeLists.txt to point it to main_xxx.cpp
cmake .
make
{% endhighlight %}

Now, you can run your code with the executable `girdap` which is placed in the `girdap/bin` directory.

{% highlight bash linenos %}
bin/girdap
{% endhighlight %}

Examples/tutorials can be built with the make command:

{% highlight bash linenos %}
make div          # div is the name of the example 
{% endhighlight %}

On a successful build, the executable will be placed in `girdap/bin` directory with the example's name. For the example above, executable file will be named as `div`.  

{% include links.html %}
