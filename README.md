<h2> RELOCDD-PY (v1.0) --- by Katherine Biegel (katherine.biegel@ucalgary.ca) </h2>

<h4> Versioning </h4>

<p> v1.0 - Released May 2022 </p>

<h4> Software Description </h4>

<p> <strong> RELOCDD-PY </strong> is a Python program for relocating earthquakes that employs various geometries of the double-difference (DD) method.  These include: </p>
<ol>
	<li> Event-Pair DD (Waldhauser and Ellsworth, 2000) </li>
	<li> Station-Pair DD (Zhang et. al., 2010) </li>
	<li> Double-Pair DD (Guo and Zhang, 2016) </li>
	<li> Combined-Pair DD (Biegel and Dettmer, 2022) </li>
</ol>

<h4> Where to Download Code </h4>

<p> This code is available to download on github: https://github.com/ucalgary-seismology (UPDATE LINK LATER). </p>

<p> To download code you can also run: </p>

	cd ~/PATH/TO/PYTHON/PACKAGES
	git clone [GITURL]

<h4> Dependencies and Installation </h4>

<p> <strong> RELOCDD-PY </strong> can be installed on any machine that supports <em> python >= <strong> 3.5 </strong></em>. </p>

<p> <strong> RELOCDD-PY </strong> contains several self-contained software functions that may require further software installation.  Depending on your usage of the program, you may not need to install all packages. </p>

<p> For simple installation we have included a setup file.  To install run: </p>

	cd ~/PATH/TO/RELOCDD/REPOSITORY
	sudo python3 setup.py install (--user)

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

<h5> IF USING TO BUILD CROSS CORRELATION DT FILES: </h5>

<p> For problems with large number of data, it is often best to run multi-process cross-correlations. While this software includes scripts to build cross-correlation differential time files for both single-process and multi-process workflows, the multi-process workflow does require usage of the <strong> MPI4PY </strong> package.  Installation instructions for that package can be found here: https://mpi4py.readthedocs.io/en/stable/install.html </p>

<h5> IF USING FOR SYNTHETIC TESTING: </h5>

<p> If you plan to use <strong> RELOCDD-PY </strong> for synthetic tests (including building the intial event catalogs), you will need to have installed the fortran program <strong> HYPOINVERSE2000 </strong> by Fred Klein. https://www.usgs.gov/software/hypoinverse-earthquake-location </p>

<p> In addition, you will have to update the <strong> utility/synthetics/hypoinv.py </strong> script <em> (line 15) </em> to include the bin path. </p>

<h4> Documentation and Tutorials </h4>

<p> Updated online documentation for this code can be found here: [URL] </p>

<p> A PDF version of documentation for <strong> RELOCDD-PY v1.0 </strong> can also be found
here: [URL] </p>

<h4> References </h4>

<h5> For this Software in Particular: </h5>

> (In Prep:) Biegel, K.M., and Jan Dettmer, 2022, Combined-type double-difference relocation, Journal, 
Volume, Pages.

> SOFTWARE DOI

<h5> For the Event-Pair DD Method: </h5>

> Waldhauser, Felix, 2001, hypoDD -- A program to compute double-difference hypocenter locations, U.S. Geological Survey Open-File Report 01-113.

> Waldhauser, Felix, and William L. Ellsworth, 2000, A double-difference earthquake location algorithm: Method and application to the northern Hayward fault, California, BSSA, 90, 1353-1368.

<h5> For the Station-Pair DD Method: </h5>

> Guo, Hao and Haijiang, Zhang, 2016, Development of double-pair double difference earthquake location algorithm for improving earthquake locations, GJI, 208(1), 333-348.

<h5> For the Double-Pair DD Method: </h5>

> Zhang, Haijiang, Naudeau, Robert M., and M. Nafi Toksoz, 2010, Locating nonvolcanic tremors beneath the San Andreas fault using a station-pair double-difference location method, GRL, 37(13).

<h4> Warranty and Licensing </h4>

<p> This software is released opensource under the <strong> GNU General Public License </strong>.  A detailed license can be seen in <strong> LICENSE.txt </strong>. </p>

<p> Copyright (c) 2022 Katherine Biegel </p>

<p> <strong> RELOCDD-PY </strong> is free software: you may redistribute or modify it under the terms of the GNU General Public License. <strong> RELOCDD-PY </strong> is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. </p>

<h4> Contributions </h4>

<p> This project is an open-source project and will continue to be updated over time.  We welcome all contributions to this project! </p>

<p> <strong> <em> Last Update: May 2022 </em> </strong> </p>

<h4> Code Support and Contact </h4>

<p> If you have any issues with installation or usage of <strong> RELOCDD-PY </strong>, please first check the <em> FAQ </em> section of the <strong> RELOCDD-PY </strong> documentation: [INSERT URL] </p>

<p> If issues persist, please check the "Issues" tab in this repository and you may submit an issue request as well. </p>

<p> If in doubt, please reach out! </p>

<p> Katherine Biegel (katherine.biegel@ucalgary.ca) </p>
