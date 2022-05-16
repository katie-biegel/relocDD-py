import numpy as np
from utility.universal.misc import delaz
#from memory_profiler import profile

"""
Function definitions for rayTrace subroutines taken from HypoDD v1.3
---
Includes the functions:
	direct: 	Computes upward direct ray
	vmodel: 	Extract information from vel. model inputs
	tiddid: 	Traveltime intercept and critical dist. for ray
	refract: 	Computes the refracted ray
	ttime:		Computes fast traveltime from direct and refract
	partials:	Returns traveltimes, partial derivatives for all rays
---
"""
def direct(nl, v, vsq, thk, jl, tkj, delta, depth):
	"""
	This function computes the traveltime for the direct
	(upward-departing) ray and replace the direct1 
	subroutine from HypoDD v1.3
	:::
    direct predicts the travel time, the sine of the 
    takeoff angle, and the horizontal distance of travel in
    the event layer.  
    The receiver must be located at the top of layer 1 and 
    the event must be located below layer 1.  Low velocity
    layers are permitted.
    :::
    The basic scheme adopted here is the method of false
    position.  (see acton, 1970, 'numerical methods that work,' for
    example.)  First, the characteristics of the fastest layer
    between the event and the surface are determined.
	:::
	PARAMETERS:
	nl (int) ---- Number of layers in v model
	v (float array) ---- Layer wave speeds with length nl
	vsq (float array) ---- Squares of wave speed with length nl
	thk (float array) ---- Layer thicknesses with length nl
	jl (int) ---- Index of event layer
	tkj (float) ---- Depth of event within layer jl
	delta (float) ---- Epicentral distance
	real (depth) ---- Event depth
	:::
	RETURNS:
	tdir (float) ---- Direct-ray travel time
	u (float) ---- Sine of the take-off angle
	x (float) ---- Horizontal travel distance in event layer
	:::
	"""
	delt = 0
	# Check if the event lies in the surface layers
	if jl==0:
		# Surface Layer event
		r = np.sqrt(depth*depth + delta*delta)
		tdir = r/v[0]
		u = delta/r
		x = delta
		return tdir,u,x
	# Otherwise find the fastest layer, lmax
	lmax = jl
	tklmax = tkj
	vlmax = v[jl]
	for l in range(0,jl):
		if v[l]>vlmax:
			lmax = l
			tklmax = thk[l]
			vlmax = v[l]
	# Min layer thickness
	if tklmax <= 0.05:
		tklmax = 0.05
	# Find initial bounds on sine of takeoff angle
	ua = (v[jl]/vlmax)*delta/np.sqrt(delta*delta + depth*depth)
	ub = (v[jl]/vlmax)*delta/np.sqrt(delta*delta + tklmax*tklmax)
	# Calculate Horizontal Travel Distances
	uasq = ua*ua
	ubsq = ub*ub
	# Max sine takeoff angle
	if ubsq >= 1:
		ubsq = 0.9999
	if uasq >= 1:
		uasq = 0.9999
	xa = tkj*ua/np.sqrt(1.0-uasq)
	if lmax == jl:
		xb = delta
	else:
		xb = tkj*ub/np.sqrt(1.0-ubsq)
	dela = xa
	delb = xb
	for l in range(0,jl):
		dela = dela + thk[l]*ua/np.sqrt(vsq[jl]/vsq[l]-uasq)
		ubdiv = np.sqrt(vsq[jl]/vsq[l]-ubsq)
		if ubdiv <= 1.e-20:
			ubdiv = 1.e-20
		delb = delb+thk[l]*ub/np.sqrt(vsq[jl]/vsq[l]-ubsq)
	# Loop to find the zero of delt-delta by the method of false position
	for count in range(0,25):
		if (delb-dela) < 0.02:
			x = 0.5*(xa+xb)
			u = x/np.sqrt(x*x + tkj*tkj)
			usq = u*u
			break
		else:
			x = xa+(delta-dela)*(xb-xa)/(delb-dela)
			u = x/np.sqrt(x*x + tkj*tkj)
			usq = u*u
			delt = x
			for l in range(0,jl):
				delt = delt+thk[l]*u/np.sqrt(vsq[jl]/vsq[l]-usq)
			xtest = delt-delta
			if abs(xtest) < 0.02:
				break
			if xtest < 0.0:
				xa = x
				dela = delt
			else:
				xb = x
				delb = delt
	# Calculate direct-ray travel time
	tdir = np.sqrt(x*x + tkj*tkj)/v[jl]
	for l in range(0,jl):
		tdir = tdir+thk[l]*v[jl]/(vsq[l]*np.sqrt(vsq[jl]/vsq[l]-usq))
	tdir = tdir - (u/v[jl])*(delt-delta)

	return tdir, u, x

def vmodel(nl, v, top, depth):
	"""
	Extract needed information from the velocity model
	:::
	PARAMETERS:
	nl (int) ---- Number of layers in velocity model
	v[nl] (float array) ---- Velocity in each layer
	top[nl] (float array) ---- Depth to top of layer
	depth (float) ---- Focal depth of source in km
	:::
	RETURNS:
	vsq[nl] (float array) ---- Squared velocities
	thk[nl] (float array) ---- Thicknesses of layers
	jl (int) ---- Event layer
	tkj (float) ---- Depth of event in event layer
	:::
	"""
	jl = nl-1
	for i in range(nl-1,-1,-1):
		if depth <= top[i]:
			jl = i-1
	# Depth from top of layer to source
	tkj = depth-top[jl]

	return jl,tkj

def tiddid(jl, nl, v, vsq, thk):
	"""
	This function determines the travel time intercept and 
	critical distance for a seismic ray in a layered earth model
	originating at the top of layer jl, refracting in layer m,
	and terminating at the top of layer 1.
	:::
	PARAMETERS:
	jl (int) ---- Event layer
	nl (int) ---- Number of layers
	v[nl] (float array) ---- Velocity in each layer
	vsq[nl] (float array) ---- Squared velocities
	thk[nl] (float array) ---- Thicknesses of each layer
	:::
	RETURNS:
	tid[20] (float array) ---- Travel time intercept for refraction in each layer
	did[20] (float array) ---- Critical distance each layer
	:::
	"""
	tid = np.zeros(20)
	did = np.zeros(20)	
	j1=jl+1
	for m in range(j1,nl):
		tid[m] = 0.
		did[m] = 0.
		tid1 = 0.
		tid2 = 0.
		did1 = 0.
		did2 = 0.
		for l in range(0,m):
			if vsq[m] <= vsq[l]:
				# m is a low velocity layer
				# set tid and did to large values
				tid[m] = 100000.
				did[m] = 100000.
			else:
				sqt = np.sqrt(vsq[m]-vsq[l])
				tim = thk[l]*sqt/(v[l]*v[m])
				dimm = thk[l]*v[l]/sqt
				if l < jl:
					# sum for layers above event layers
					tid1 = tid1 + tim
					did1 = did1 + dimm
				else:
					# sum for layers below and including event layer
					tid2 = tid2 + tim
					did2 = did2 + dimm
		if tid[m] != 100000.:
			tid[m] = tid1 + 2.*tid2
			did[m] = did1 + 2.*did2
	return tid,did

def refract(nl, v, vsq, thk, jl, tkj, delta):
	"""
	Find "refracted" ray with smallest travel time
	:::	
	For refracted layers in a layered earth model, refract 
	determines the fastest travel time, tref, the layer 
	in which the fastest ray is refracted, kk, the 
	critical distance for refraction in that layer, 
	didjkk, and an upper bound on delta for which a 
	direct ray can be a first arrival, xovmax. Refract 
	allows for the possibility of low velocity layers.
	:::
	PARAMETERS:
	nl (int) ---- Number of layers
	v[nl] (float array) ---- Velocity of layers
	vsq[nl] (float array) ---- Squared velocities
	thk[nl] (float array) ---- Thickness of layers
	jl (int) ---- Event layer
	tkj (float) ---- Depth of event in event layer
	delta (float) ---- Horizontal distance between event and receiver
	:::
	RETURNS:
	kk (int) ---- Refracting layer for fasted refracted ray
	tref (float) ---- Travel time of fasted refracted layer
	xovmax (float) ---- An upper bound on delta for which the 
					direct ray can be the first arrival
	:::
	"""	
	lx = int(0)
	tr = np.zeros(nl)
	tinj = np.zeros(nl)
	didj = np.zeros(nl)
	# Determine tref, kk, didjkk
	tid,did = tiddid(jl,nl,v,vsq,thk)
	tref = 100000.
	j1 = jl+1
	for m in range(j1,nl):
		if tid[m] == 100000.:
			tr[m] = 100000.
		else:
			sqt = np.sqrt(vsq[m] - vsq[jl])
			tinj[m] = tid[m] - tkj*sqt/(v[m]*v[jl])
			didj[m] = did[m] - tkj*v[jl]/sqt
			tr[m] = tinj[m] + delta/v[m]
			if didj[m] > delta:
				tr[m] = 100000.
		if tr[m] < tref:
			tref = tr[m]
			kk = m
	# If there's no refracted ray
	if tref == 100000.:
		xovmax = 100000.
		kk = 0
		return kk, tref, xovmax
	# If threre's a refracted ray, determine xovmx:
	# Find lx (the 1st layer below the event layer which is 
	# not a low velocity layer)
	m = jl+1
	if tid[m] == 100000.:
		m = m+1
	else:
		lx = m
	# Check if the event is in the first layer
	if jl == 0:
		xovmax = tinj[lx]*v[lx]*v[0]/(v[lx] - v[0])
		return kk, tref, xovmax
	m = jl

	tid[m] = 0.
	for l in range(0,m):
		if vsq[m] <= vsq[l]:
			tid[m] = 100000.
		else:
			sqt = np.sqrt(vsq[m]-vsq[l])
			tim = thk[l]*sqt/(v[l]*v[m])
			tid[m] = tid[m]+tim
	m = m-1

	if tid[m+1] < 100000.:
		jx = m+1
		xovmax = (tinj[lx]-tid[jx])*v[lx]*v[jx]/(v[lx]-v[jx])
	else:
		xovmax = tinj[lx]*v[lx]*v[0]/(v[lx]-v[0])

	return kk, tref, xovmax

def ttime(delta, depth, nl, v, top, vsq, thk):
	"""
	This function determines the fastest traveltime between
	a source at depth=depth and a receiver at distance=delta(km)
	:::
	PARAMETERS:
	delta (float) ---- Epicentral distance in km
	depth (float) ---- Focal depth of source in km
	nl (int) ---- Number of layers in velocity model
	v[nl] (float array) ---- Velocity in each layer
	top[nl] (float array) ---- Depth to top of layer
	:::
	RETURNS:
	t (float) ---- Minimum traveltime
	ain (float) ---- Angle of emergences at the source
	:::
	"""
	jl,tkj = vmodel(nl,v,top,depth)
	kk,tref,xovmax = refract(nl,v,vsq,thk,jl,tkj,delta)
	# if delta <= xovmax, call direct to find the direct
	# ray traveltime otherwise tref is the minimum traveltime
	t = tref
	if kk > 0:
		u=v[jl]/v[kk]
		ain = np.arcsin(u)*57.2958

	if delta <= xovmax:
		tdir,u,x = direct(nl,v,vsq,thk,jl,tkj,delta,depth)
		# compare traveltimes
		if tref > tdir:
			t = tdir
			ain = 180. - np.arcsin(u)*57.2958

	return t,ain

#@profile(stream=open('mem_logs/partials.mem','w+'))
def partials(nsrc, src_cusp, src_lat, src_lon, src_dep, 
			 nsta, sta_lab, sta_lat, sta_lon, 
			 mod_nl, mod_ratio, mod_v, mod_top, fn_srcpar='rayTrace.src', return_all=False):
	"""
	This function returns traveltimes (P&S) and partial time derivatives for sources
	and receiver pairs for a give homogeneous layered vmodel.
	:::
	PARAMETERS:
	nsrc (int) ---- Number of sources
	src_cusp[nsrc] (int array) ---- Integer event IDS for sources
	src_lat[nsrc] (float array) ---- Source latitudes
	src_lon[nsrc] (float array) ---- Source longitudes
	src_dep[nsrc] (float array) ---- Source depths in km
	nsta (int) ---- Number of stations
	sta_lab[nsta] (str array) ---- Station names
	sta_lat[nsta] (float array) ---- Station latitudes
	sta_lon[nsta] (float array) ---- Station longitudes
	mod_nl (int) ---- Number of layers in velocity model
	mod_ratio (float) ---- VP/VS ratio
	mod_v[mod_nl] (float array) ---- Layer P velocities (km/s)
	mod_top[mod_nl] (float array) ---- Depth to top of layer (km)
	fn_srcpar (str) ---- Source parameter file locations defaults to 'rayTrace.src'
	:::
	RETURNS:
	tmp_ttp[nsta,nsrc] (float array) ---- P traveltime for all station-event combos
	tmp_tts[nsta,nsrc] (float array) ---- S traveltime for all station-event combos
	tmp_xp[nsta,nsrc] (float array) ---- X partial derivative
	tmp_yp[nsta,nsrc] (float array) ---- Y partial derivative
	tmp_zp[nsta,nsrc] (float array) ---- Z partial derivative
	or also
	dist_all[nsta,nsrc] (float array) ---- Event-station distances
	az_all[nsta,nsrc] (float array) ---- Event-station azimuth
	ang_all[nsta,nsrc] (float array) ---- Take-off angle of ray
	:::
	"""
	# Initialise output arrays
	tmp_ttp = np.zeros((nsta,nsrc),dtype='float64')
	tmp_tts = np.zeros((nsta,nsrc),dtype='float64')
	tmp_xp = np.zeros((nsta,nsrc),dtype='float64')
	tmp_yp = np.zeros((nsta,nsrc),dtype='float64')
	tmp_zp = np.zeros((nsta,nsrc),dtype='float64')

	# Open src file
	srcpar = open(fn_srcpar,'w')

	# Make sure hypocenters don't fall on boundaries
	for i in range(0,nsrc):
		for j in range(0,mod_nl):
			if abs(src_dep[i] - mod_top[j]) < 0.0001:
				src_dep[i] = src_dep[i]-0.001 # Move by 1cm

	# Get S velocity model
	vs = np.zeros(mod_nl,dtype='float')
	for i in range(0,mod_nl):
		vs[i] = mod_v[i]/mod_ratio

	# Calculate layer thicknesses and vsq
	vsqp = np.zeros(mod_nl,dtype='float')
	vsqs = np.zeros(mod_nl,dtype='float')
	for i in range(0,mod_nl):
		vsqp[i] = mod_v[i]*mod_v[i]
		vsqs[i] = vs[i]*vs[i]

	# Layer thicknesses
	thk = np.zeros(mod_nl,dtype='float')
	for i in range(0,mod_nl-1):
		thk[i] = mod_top[i+1] - mod_top[i]

	# Compute epicentral distances, azimuths, angles of incidence,
	# and P/S trave time from sources to stations
	pi = float(3.141593) # Define to single precsision for continuity sake

	#############################
	#############################
	# THIS LOOP PUSH TO GPU
	# For fastest computation
	#############################
	#############################

	# Added for hypoinverse synthetic model generation
	if return_all:
		dist_all = np.zeros((nsta,nsrc),dtype='float')
		ang_all = np.zeros((nsta,nsrc),dtype='float')
		az_all = np.zeros((nsta,nsrc),dtype='float')

	for i in range(0,nsta):
		for j in range(0,nsrc):
			delt, dist, az = delaz(src_lat[j],src_lon[j],sta_lat[i],sta_lon[i])
			if return_all:
				dist_all[i,j] = dist
				az_all[i,j] = az
			# 1D ray tracing
			tmp_ttp[i,j],ain = ttime(dist,src_dep[j],mod_nl,mod_v,mod_top,vsqp,thk)
			tmp_tts[i,j],ain = ttime(dist,src_dep[j],mod_nl,vs,mod_top,vsqs,thk)
			if return_all:
				ang_all[i,j] = ain
			# Determine wave speed at the hypocenter
			for k in range(0,mod_nl):
				if src_dep[j] <= mod_top[k]:
					break
			# Depth Derivatives
			tmp_zp[i,j] = np.cos((ain*pi)/180.)/mod_v[k-1]
			# Epicentral Derivatives
			tmp_xp[i,j] = (np.sin((ain*pi)/180.)*np.cos(((az-90.)*pi)/180.))/mod_v[k-1]
			tmp_yp[i,j] = (np.sin((ain*pi)/180.)*np.cos((az*pi)/180.))/mod_v[k-1]

			# Write to source parameter file
			srcpar.write('%13g %16.8f %16.8f %13s %16.8f %16.8f %16.8f %16.8f\n' % 
				(src_cusp[j],tmp_ttp[i,j],tmp_tts[i,j],sta_lab[i],dist,az,ain,src_dep[j]))
	
	############################
	############################
	srcpar.close()

	if return_all: # Return all is dist/angles needed
		return tmp_ttp,tmp_tts,tmp_xp,tmp_yp,tmp_zp,dist_all,az_all,ang_all
	else:		   # Else return traveltimes and partials only
		return tmp_ttp,tmp_tts,tmp_xp,tmp_yp,tmp_zp
