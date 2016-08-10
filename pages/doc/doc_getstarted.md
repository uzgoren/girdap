---
title: Getting started
keywords: about
tags: [overview]
sidebar: mydoc_sidebar
permalink: doc_getstarted.html
summary: girdap is a c++ based object oriented library for multiphysics simulations on self-managed grids. 
---

## Download or clone _girdap_

First download or clone _girdap_ from the [Github repo](https://github.com/uzgoren/girdap).

{% include note.html content="It is recommended to download `develop` branch rather than the `master` branch as many described functionalities are available in `develop` branch at this point." %}

## Build this theme

{% include warning.html content="Requires `cmake`" %}

Extract the package into a directory. The girdap's base folder is named as `girdap` as default. It does include two subdirectories; `src` and `include`. `src` contains different versions of c++'s `main()` function; each with a different purpose. These are files those utilize _girdap_'s functionality; which are available in the `include` directory. Develop your own or modify one of the `main_xxx.cpp` files as the driver file and make sure that `CMakeLists.txt` file to make sure that line at the end that starts with `add_executable` points to your driver file (`main_xxx.cpp`) in `girdap/src` folder. Now, you can go ahead with `cmake` and `make` commands as usual.

{% highlight bash linenos %}
cd dir_of_your_choice
tar -xzvf girdap.tar.gz
cd girdap
# Make changes to the main_xxx.cpp
# modify CMakeLists.txt to point it to main_xxx.cpp
cmake .
make
{% endhighlight %}

Now, you can run your code with the executable named as `girdap`:

{% highlight bash linenos %}
./girdap
{% endhighlight %}

{% include links.html %}
