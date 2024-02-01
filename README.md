<h2> RELOCDD-PY (v1.0) </h2>

<h4> <u> Versioning </u> </h4>

<p> v1.0 - Released Feb. 2024 </p>

<h4> <u> Software Description </u> </h4>

<p> <strong> RELOCDD-PY </strong> is a Python program for relocating earthquakes that employs various geometries of the double-difference (DD) method.  These include: </p>
<ol>
	<li> Event-Pair DD (Waldhauser and Ellsworth, BSSA, 2000) </li>
	<li> Station-Pair DD (Zhang et. al., GRL, 2010) </li>
	<li> Double-Pair DD (Guo and Zhang, GJI, 2016) </li>
</ol>

<h4> <u> Where to Download Code </u> </h4>

<p> This code is available to download on github: https://github.com/katie-biegel/relocDD-py </p>

<p> To download code you can also run: </p>

	cd ~/PATH/TO/PYTHON/PACKAGES
	git clone [GITURL]

<h4> <u> Dependencies and Installation </u> </h4>

<p> <strong> RELOCDD-PY </strong> can be installed on any machine that supports <em> python >= <strong> 3.5 </strong></em>. </p>

<p> <strong> RELOCDD-PY </strong> contains several self-contained software functions that may require further software installation.  Depending on your usage of the program, you may not need to install all packages. </p>

<p> To install by hand you need the following packages:

<ul>
	<li> <strong> python </strong> </li>
	<li> os </li>
	<li> numpy </li>
	<li> datetime </li>
	<li> sys </li>
	<li> subprocess </li>
	<li> matplotlib </li>
	<li> scipy </li>
	<li> glob </li>
	<li> obspy </li>
	<li> utm </li>
</ul>

<h5> <u> IF USING TO BUILD CROSS CORRELATION DT FILES: </u> </h5>

<p> For problems with large number of data, it is often best to run multi-process cross-correlations. While this software includes scripts to build cross-correlation differential time files for both single-process and multi-process workflows, the multi-process workflow does require usage of the <strong> MPI4PY </strong> package.  Installation instructions for that package can be found here: https://mpi4py.readthedocs.io/en/stable/install.html </p>

<h5> <u> IF USING FOR SYNTHETIC TESTING: </u> </h5>

<p> If you plan to use <strong> RELOCDD-PY </strong> for synthetic tests (including building the intial event catalogs), you will need to have installed the fortran program <strong> HYPOINVERSE2000 </strong> by Fred Klein. https://www.usgs.gov/software/hypoinverse-earthquake-location </p>

<p> In addition, you will have to update the <strong> utility/synthetics/hypoinv.py </strong> script <em> (line 15) </em> to include the bin path. </p>

<h4> <u> References </u> </h4>

<h5> <u> For this Software in Particular: </u> </h5>

> K.M. Biegel, & Dettmer, J. (2024). relocDD-py (v1.0). Zenodo. https://doi.org/10.5281/zenodo.10607406.

>![image](https://github.com/katie-biegel/relocDD-py/assets/32553479/ebf91cfe-c818-4e93-a6a9-aea2596a4706)

<h5> <u> For the Event-Pair DD Method: </u> </h5>

> Waldhauser, Felix, 2001, hypoDD -- A program to compute double-difference hypocenter locations, U.S. Geological Survey Open-File Report 01-113.

> Waldhauser, Felix, and William L. Ellsworth, 2000, A double-difference earthquake location algorithm: Method and application to the northern Hayward fault, California, BSSA, 90, 1353-1368.

<h5> <u> For the Station-Pair DD Method: </u> </h5>

> Zhang, Haijiang, Naudeau, Robert M., and M. Nafi Toksoz, 2010, Locating nonvolcanic tremors beneath the San Andreas fault using a station-pair double-difference location method, GRL, 37(13).

<h5> <u> For the Double-Pair DD Method: </u> </h5>

> Guo, Hao and Haijiang, Zhang, 2016, Development of double-pair double difference earthquake location algorithm for improving earthquake locations, GJI, 208(1), 333-348.

<h4> <u> Warranty and Licensing </u> </h4>

<p> This software is released opensource under the <strong> GNU General Public License </strong>.  A detailed license can be seen in <strong> LICENSE.txt </strong>. </p>

<p> Copyright (c) 2024 Katherine Biegel </p>

<p> <strong> RELOCDD-PY </strong> is free software: you may redistribute or modify it under the terms of the GNU General Public License. <strong> RELOCDD-PY </strong> is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. </p>

<h4> <u> Contributions </u> </h4>

<p> This project is an open-source project and will continue to be updated over time.  We welcome all contributions to this project! </p>

<p> <strong> <em> Last Update: Feb. 2024 </em> </strong> </p>

<h4> <u> Code Support and Contact </u> </h4>

<p> If you have any issues with installation or usage of <strong> RELOCDD-PY </strong>, please reach out! </p>

<p> Katherine Biegel (katherine.biegel@ucalgary.ca) </p>
