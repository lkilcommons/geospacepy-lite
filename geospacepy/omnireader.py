import sys, os, copy, textwrap, datetime, subprocess, ftplib
from spacepy import pycdf
from geospacepy import special_datetime
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pp
import scipy.interpolate as interpolate
#Variables in 1 Hour CDFS
# ABS_B: CDF_REAL4 [4344]
# AE: CDF_INT4 [4344]
# AL_INDEX: CDF_INT4 [4344]
# AP_INDEX: CDF_INT4 [4344]
# AU_INDEX: CDF_INT4 [4344]
# BX_GSE: CDF_REAL4 [4344]
# BY_GSE: CDF_REAL4 [4344]
# BY_GSM: CDF_REAL4 [4344]
# BZ_GSE: CDF_REAL4 [4344]
# BZ_GSM: CDF_REAL4 [4344]
# Beta: CDF_REAL4 [4344]
# DST: CDF_INT4 [4344]
# Day: CDF_INT4 [4344]
# E: CDF_REAL4 [4344]
# Epoch: CDF_EPOCH [4344]
# F: CDF_REAL4 [4344]
# F10_INDEX: CDF_REAL4 [4344]
# HR: CDF_INT4 [4344]
# IMF: CDF_INT4 [4344]
# IMF_PTS: CDF_INT4 [4344]
# KP: CDF_INT4 [4344]
# MFLX: CDF_INT4 [4344]
# Mach_num: CDF_REAL4 [4344]
# Mgs_mach_num: CDF_REAL4 [4344]
# N: CDF_REAL4 [4344]
# PC_N_INDEX: CDF_REAL4 [4344]
# PHI-V: CDF_REAL4 [4344]
# PHI_AV: CDF_REAL4 [4344]
# PLS: CDF_INT4 [4344]
# PLS_PTS: CDF_INT4 [4344]
# PR-FLX_1: CDF_REAL8 [4344]
# PR-FLX_10: CDF_REAL4 [4344]
# PR-FLX_2: CDF_REAL4 [4344]
# PR-FLX_30: CDF_REAL4 [4344]
# PR-FLX_4: CDF_REAL4 [4344]
# PR-FLX_60: CDF_REAL4 [4344]
# Pressure: CDF_REAL4 [4344]
# R: CDF_INT4 [4344]
# Ratio: CDF_REAL4 [4344]
# Rot#: CDF_INT4 [4344]
# SIGMA-ABS_B: CDF_REAL4 [4344]
# SIGMA-B: CDF_REAL4 [4344]
# SIGMA-Bx: CDF_REAL4 [4344]
# SIGMA-By: CDF_REAL4 [4344]
# SIGMA-Bz: CDF_REAL4 [4344]
# SIGMA-N: CDF_REAL4 [4344]
# SIGMA-PHI-V: CDF_REAL4 [4344]
# SIGMA-T: CDF_REAL4 [4344]
# SIGMA-THETA-V: CDF_REAL4 [4344]
# SIGMA-V: CDF_REAL4 [4344]
# SIGMA-ratio: CDF_REAL4 [4344]
# T: CDF_REAL4 [4344]
# THETA-V: CDF_REAL4 [4344]
# THETA_AV: CDF_REAL4 [4344]
# V: CDF_REAL4 [4344]
# YR: CDF_INT4 [4344]

# Variables in 5 Minute CDFs
# <CDF:
# AE_INDEX: CDF_INT4 [8928]
# AL_INDEX: CDF_INT4 [8928]
# ASY_D: CDF_INT4 [8928]
# ASY_H: CDF_INT4 [8928]
# AU_INDEX: CDF_INT4 [8928]
# BSN_x: CDF_REAL4 [8928]
# BSN_y: CDF_REAL4 [8928]
# BSN_z: CDF_REAL4 [8928]
# BX_GSE: CDF_REAL4 [8928]
# BY_GSE: CDF_REAL4 [8928]
# BY_GSM: CDF_REAL4 [8928]
# BZ_GSE: CDF_REAL4 [8928]
# BZ_GSM: CDF_REAL4 [8928]
# Beta: CDF_REAL4 [8928]
# Day: CDF_INT4 [8928]
# E: CDF_REAL4 [8928]
# Epoch: CDF_EPOCH [8928]
# F: CDF_REAL4 [8928]
# HR: CDF_INT4 [8928]
# IMF: CDF_INT4 [8928]
# IMF_PTS: CDF_INT4 [8928]
# Mach_num: CDF_REAL4 [8928]
# Mgs_mach_num: CDF_REAL4 [8928]
# Minute: CDF_INT4 [8928]
# PC_N_INDEX: CDF_REAL4 [8928]
# PLS: CDF_INT4 [8928]
# PLS_PTS: CDF_INT4 [8928]
# PR-FLX_10: CDF_REAL4 [8928]
# PR-FLX_30: CDF_REAL4 [8928]
# PR-FLX_60: CDF_REAL4 [8928]
# Pressure: CDF_REAL4 [8928]
# RMS_SD_B: CDF_REAL4 [8928]
# RMS_SD_fld_vec: CDF_REAL4 [8928]
# RMS_Timeshift: CDF_INT4 [8928]
# SYM_D: CDF_INT4 [8928]
# SYM_H: CDF_INT4 [8928]
# T: CDF_REAL4 [8928]
# Time_btwn_obs: CDF_INT4 [8928]
# Timeshift: CDF_INT4 [8928]
# Vx: CDF_REAL4 [8928]
# Vy: CDF_REAL4 [8928]
# Vz: CDF_REAL4 [8928]
# YR: CDF_INT4 [8928]
# flow_speed: CDF_REAL4 [8928]
# percent_interp: CDF_INT4 [8928]
# proton_density: CDF_REAL4 [8928]
# x: CDF_REAL4 [8928]
# y: CDF_REAL4 [8928]
# z: CDF_REAL4 [8928]
# >
import geospacepy
localdir = geospacepy.config['omnireader']['local_cdf_dir']

class omni_downloader(object):
	def __init__(self):
		self.localdir = localdir
		self.ftpserv = 'spdf.gsfc.nasa.gov'
		self.ftpdir = '/pub/data/omni/omni_cdaweb/'	
		self.cadence_subdir = {'hourly':'hourly','5min':'hro_5min','1min':'hro_1min'}
		#Hourly CDF are every six months, 5 minute are every month as are 1 min
		self.filename_gen = {'hourly':lambda dt: 'omni2_h0_mrg1hr_%d%.2d01_v01.cdf' % (dt.year,1 if dt.month < 7 else 7),
							 '5min':lambda dt: 'omni_hro_5min_%d%.2d01_v01.cdf' % (dt.year,dt.month),
							 '1min':lambda dt: 'omni_hro_1min_%d%.2d01_v01.cdf' % (dt.year,dt.month) }
	
	def get_cdf(self,dt,cadence):
		fn = self.filename_gen[cadence](dt)
		localfn = os.path.join(self.localdir,fn)
		if not os.path.exists(localfn):
			ftp = ftplib.FTP(self.ftpserv)
			print 'Connecting to OMNIWeb FTP server %s' % (self.ftpserv)
			ftp.connect()
			ftp.login()
			#Change directory
			ftp.cwd(self.ftpdir+'/'+self.cadence_subdir[cadence]+'/%d' % (dt.year))
			print 'Downloading file %s' % (self.ftpdir+'/'+self.cadence_subdir[cadence]+'/'+'%d' % (dt.year)+'/'+fn)
			with open(localfn,'wb') as f:
				ftp.retrbinary('RETR ' + fn, f.write)
			print "Saved as %s" % (localfn)
			ftp.quit()
		return pycdf.CDF(localfn)	

class omni_derived_var(object):
	"""
	A variable that can be defined as a function
	of CDF variables, that we would like to be accessible via the same __getitem__
	interface as the variables contained in the CDF
	"""
	def __init__(self,oi):
		self.oi = oi # associated omni_interval
		self.varvals = None # To hold computed so we don't recompute without needing to
		self.attrs = dict()

class borovsky(omni_derived_var):
	"""Borovsky solar wind coupling function"""
	def __call__(self):
		#Short circut if already computed
		if self.varvals is not None:
			return self.varvals

		oi = self.oi
		#Deal with names that differ between cadences
		(densvar,vswvar) = ('N','V') if oi.cadence == 'hourly' else ('proton_density','flow_speed')
		bx,by,bz = oi['BX_GSE'],oi['BY_GSM'],oi['BZ_GSM']
		n,vsw,mach = oi[densvar],oi[vswvar],oi['Mach_num']
		bt = np.sqrt(by**2+bz**2)
		#Compute IMF clock angle
		ca = np.arccos(bz/bt)
		borovsky = 3.29e-2*(np.sin(ca/2)**2)*np.sqrt(n)*vsw**2*mach**(-.18)*np.exp(np.sqrt((mach/3.42)))
		self.attrs['CATDESC'] = 'Borovsky Solar Wind Coupling Function'
		self.attrs['UNITS'] = 'nT km/s'
		self.varvals=borovsky
		return borovsky

class newell(omni_derived_var):
	"""Newell emperical solar wind coupling function"""
	def __call__(self):
		#Short circut if already computed
		if self.varvals is not None:
			return self.varvals

		oi = self.oi
		#Deal with names that differ between cadences
		vswvar = 'V' if oi.cadence == 'hourly' else 'flow_speed'
		bx,by,bz = oi['BX_GSE'],oi['BY_GSM'],oi['BZ_GSM']
		vsw,mach = oi[vswvar],oi['Mach_num']
		bt = np.sqrt(by**2+bz**2)
		#Compute IMF clock angle
		ca = np.arccos(bz/bt)

		newell = (vsw*1000.)**(4./3)*(bt*1.0e-9)**(2./3)*(np.sin(ca/2))**(8./3);
		self.attrs['CATDESC'] = 'Newell Solar Wind Coupling Function'
		self.attrs['UNITS'] = 'm/s^(4/3) T^(2/3)'
		self.varvals=newell
		return newell

class knippjh(omni_derived_var):
	"""Knipp Joule Heating Index (Old Version)"""
	def __call__(self):
		#Short circut if already computed
		if self.varvals is not None:
			return self.varvals

		oi = self.oi
		#Computes the Joule heating index from Knipp, Tobiska, Emery,
		#Direct and Indirect Thermospheric Heating Sources For Solar Cycles 21-23,
		#Solar Physics
		dt = oi['Epoch'][:].flatten()
		doy = special_datetime.datetimearr2doy(dt).flatten()

		#Leap year modifier (adds one if this is a leap year)
		#---NOTE: this implementation does not an edge cases:
		#	1. An omni interval runs between two years one of which is a leap
		# The workaround (computing a lymod value for each dt) would degrade performance
		# and the effect is small
		# so I chose not to address it. -LMK
		if np.mod(dt[0].year,4) == 0:
			lymod = 1.
		else:
			lymod = 0.
		
		#Paper uses absolute value of pcn and dst
		PC = np.abs(self.oi['PC_N_INDEX']).flatten()
		Dst = np.abs(self.oi['DST' if oi.cadence == 'hourly' else 'SYM_H']).flatten()

		jhindex = np.zeros_like(PC)
		jhindex.fill(np.nan)

		# 4 Seasons
		annual = np.logical_or(doy > 335.+lymod, doy < 31.)
		jhindex[annual] = 24.89*PC[annual] + 3.41*PC[annual]**2 + .41*Dst[annual] + .0015*Dst[annual]**2
		
		#Winter is 21 October (294 for non-leap, 295 for leap) to 20
		#February (DOY 51 for both non-leap and leap) 
		winter = np.logical_or(doy > 294.+lymod, doy < 51.)
		jhindex[winter] = 13.36*PC[winter] + 5.08*PC[winter]**2 + .47*Dst[winter] + .0011*Dst[winter]**2
	   
		#Summer is 21 April (DOY 111, 112 leap) - 20 August (DOY 232, 233 leap)
		summer = np.logical_and(doy > 111.+lymod , doy < 232.+lymod)
		jhindex[summer] = 29.27*PC[summer] + 8.18*PC[summer]**2 - .04*Dst[summer] + .0126*Dst[summer]**2
		
		#Equinox is 21 Feb (DOY 51) - 20 Apr (DOY 110, 111 leap)
		# and 21 Aug (233, 234 leap) - 20 Oct (293 non-leap, 294 leap)
		equinox = np.logical_or(
			np.logical_and(doy > 51.+lymod , doy < 110.+lymod),
			np.logical_and(doy > 233.+lymod , doy < 293.+lymod)
			)

		jhindex[equinox] = 29.14*PC[equinox] + 2.54*PC[equinox]**2 + .21*Dst[equinox] + .0023*Dst[equinox]**2
		self.attrs['UNITS']='GW'
		self.attrs['CATDESC'] = 'Knipp Joule Heating Index'
		self.varvals = jhindex
		return jhindex

class omni_interval(object):
	def __init__(self,startdt,enddt,cadence,silent=False):
		#Just handles the possiblilty of having a read running between two CDFs 
		self.dwnldr = omni_downloader()
		self.silent = silent #No messages
		self.cadence = cadence
		self.startdt = startdt
		self.enddt = enddt
		self.cdfs = [self.dwnldr.get_cdf(startdt,cadence)]
		self.attrs = self.cdfs[-1].attrs #Mirror the global attributes for convenience
		self.transforms = dict() #Functions which transform data automatically on __getitem__
		#Find the index corresponding to the first value larger than startdt
		self.si = np.searchsorted(self.cdfs[0]['Epoch'][:],startdt)
		while self.cdfs[-1]['Epoch'][-1] < enddt:
			#Keep adding CDFs until we span the entire range
			self.cdfs.append(self.dwnldr.get_cdf(self.cdfs[-1]['Epoch'][-1]+datetime.timedelta(days=1),cadence))
		#Find the first index larger than the enddt in the last CDF
		self.ei = np.searchsorted(self.cdfs[-1]['Epoch'][:],enddt)
		if not self.silent:
			print "Created interval between %s and %s, cadence %s, start index %d, end index %d" % (self.startdt.strftime('%Y-%m-%d'),
				self.enddt.strftime('%Y-%m-%d'),self.cadence,self.si,self.ei)
		self.add_transform('KP',['hourly'],lambda x: x/10.,'Hourly Kp*10 -> Kp')
		#Implement computed variables
		self.computed = dict()
		self.computed['borovsky']=borovsky(self)
		self.computed['newell']=newell(self)
		self.computed['knippjh']=knippjh(self)
		
	def get_var_attr(self,var,att):
		"""Get a variable attribute"""
		if att in self.cdfs[-1][var].attrs:
			return self.cdfs[-1][var].attrs[att]
		else:
			return None

	def __getitem__(self,cdfvar):
		#If it's a derived variable go get it
		#with it's own __call__ method
		if cdfvar in self.computed:
			return self.computed[cdfvar]()

		#Attempt the getitem on all the cdfs in order
		if len(self.cdfs) > 1:
			ret = []
			for icdf,cdf in enumerate(self.cdfs):
				if icdf==0:
					ret.append(cdf[cdfvar][self.si:])
				elif icdf==len(self.cdfs)-1:
					ret.append(cdf[cdfvar][:self.ei])
				else:
					ret.append(cdf[cdfvar][:])
			data = np.concatenate(ret)
		else:
			data = self.cdfs[-1][cdfvar][self.si:self.ei]
		#Fix the fill values
		if np.isfinite(self.cdfs[-1][cdfvar].attrs['FILLVAL']):
			filled = data==self.cdfs[-1][cdfvar].attrs['FILLVAL']
			if np.count_nonzero(filled) > 0:
				data[filled]=np.nan
		#Check for transforms which need to be performed
		if cdfvar in self.transforms:
			transform = self.transforms[cdfvar]
			if self.cadence in transform['cadences']:
				if not self.silent:
					print "Applying transform %s to omni %s variable %s" % (transform['desc'],self.cadence,
						cdfvar)
				#print "Data before", data
				data = transform['fcn'](data)
				#print "Data after", data
		return data

	def add_transform(self,cdfvar,cadences,fcn,desc):
		"""
			Call some function to manipulate the returned data 
			whenever a particular CDF variable is '__getitem__'d
			Added to fix the obnoxious hourly 'KP' variable being Kp*10
			source of confusion.

			Arguments:
				cdfvar - str
					name of cdf variable
				cadences - list
					cadences which have this variable
				fcn - function
					function to call to manipulate data
				desc - description of manipulation performed

		"""
		self.transforms[cdfvar] = {'cadences':cadences,'fcn':fcn,'desc':desc}

	def __str__(self):
		return str(self.cdfs[0])

class omni_event(object):
	def __init__(self,startdt,enddt,label=None,cadence='5min'):
		self.interval = omni_interval(startdt,enddt,cadence)
		datetime2doy = lambda dt: dt.timetuple().tm_yday + dt.hour/24. + dt.minute/24./60. + dt.second/86400. + dt.microsecond/86400./1e6
		self.doy = special_datetime.datetimearr2doy(self.interval['Epoch'])
		self.jd = special_datetime.datetimearr2jd(self.interval['Epoch'])
		self.label = '%s-%s' % (startdt.strftime('%m-%d-%Y'),enddt.strftime('%m-%d-%Y')) if label is None else label
		self.interpolants = dict()
		self.attrs = self.interval.attrs	

	def __getitem__(self,*args):
		return self.interval.__getitem__(*args)

	def get_var_attr(self,var,att):
		"""Get a variable attribute from the last CDF in the interval"""
		return self.interval.get_var_attr(var,att)

	def interpolate(self,var,jd,**kwargs):
		"""Interpolate a variable to the julian dates in jd"""
		#Create an interpolant if we don't have one yet
		if var not in self.interpolants:
			print "No interpolant for variable %s, creating %d point interpolant" % (var,len(self.jd.flatten()))
			t,y = self.jd.flatten(),self.interval[var].flatten()
			g = np.isfinite(y)
			self.interpolants[var]=interpolate.PchipInterpolator(t[g],y[g])
		#Return the interpolated result
		return self.interpolants[var].__call__(jd,**kwargs)

	def close(self):
		"""Close the CDFs"""
		for cdf in self.interval.cdfs:
			cdf.close()

class omni_sea(object):
	def __init__(self,center_ymdhm_list,name=None,ndays=3,cadence='5min'):
		"""
		A class for superposed epoch analysis
		"""
		self.ndays = ndays
		self.name = name # Name of the list
		self.nevents = len(center_ymdhm_list)
		self.center_dts = [datetime.datetime(y,mo,d,h)+datetime.timedelta(minutes=m) for [y,mo,d,h,m] in center_ymdhm_list]
		self.center_jds = special_datetime.datetimearr2jd(self.center_dts).flatten()
		#Create an omnidata interval for each event
		self.events = [omni_event(center_dt-datetime.timedelta(days=ndays),center_dt+datetime.timedelta(days=ndays),cadence=cadence) for center_dt in self.center_dts ]
		#mirror the attributes of the first event's last CDF
		self.attrs = self.events[0].attrs

	def plot_individual(self,ax,var,show=False,cmap=None,text_kwargs=None,plot_kwargs=None):
		""" Plot all individual intervals labeled at their maxima"""
		if cmap == None:
			norm = mpl.colors.Normalize(vmin=0,vmax=len(self.center_jds))
			cmap = mpl.cm.get_cmap('hot')
		for i,(event,center_jd) in enumerate(zip(self.events,self.center_jds)):
			t,y = event.jd-center_jd, event[var]
			maxind = np.nanargmax(np.abs(y))
			#Allow for optional arguments to text and plot
			plot_kwargs = dict() if plot_kwargs is None else plot_kwargs
			ax.plot(t,y,'-',color=cmap(norm(i)),**plot_kwargs)
			text_kwargs = {'backgroundcolor':'grey','alpha':.7,'ha':'center'} if text_kwargs is None else text_kwargs
			ax.text(t[maxind],y[maxind],event.label,**text_kwargs)
		if show:
			pp.show()
			pp.pause(10)

	def get_var_attr(self,var,att):
		"""Get a variable attribute from the last CDF in the interval"""
		return self.events[0].get_var_attr(var,att)

	def plot_stats(self,ax,var,x=None,xstep=0.042,show=False,plot_events=False,**kwargs):
		"""Default to having 1 hour steps for interpolation"""
		if x is None:
			x = np.arange(-1*self.ndays,self.ndays+xstep,xstep)

		iy = np.zeros((len(self.center_jds),len(x)))
		for i,(event,center_jd) in enumerate(zip(self.events,self.center_jds)):
			t,y = event.jd-center_jd, event[var]
			iy[i,:] = event.interpolate(var,x+center_jd)

			#Plot points
			#Get len(x) random numbers from -.5 - .5
			jitter = np.random.rand(*x.flatten().shape)-.5
			#Scale
			jitter = jitter*xstep/3
			#Plot jittered points
			if plot_events:
				ax.plot(x.flatten()+jitter.flatten(),iy[i,:].flatten(),'.',
					color='b' if 'color' not in kwargs else kwargs['color'], 
					zorder=5.,alpha=.1)
		y_med,y_lq,y_uq = np.nanmedian(iy,axis=0),np.nanpercentile(iy,25,axis=0),np.nanpercentile(iy,75,axis=0)
		lab = '' if self.name is None else '%s: ' % (self.name)
		lab += 'Median %s Response' % (var)
		ax.plot(x,y_med,label=lab, linestyle='-',zorder=10,**kwargs)
		ax.plot(x,y_lq,linestyle=':',zorder=10,**kwargs)
		ax.plot(x,y_uq,linestyle=':',zorder=10,**kwargs)
		#Put units on the y axis
		un = self.get_var_attr(var,'UNITS')
		un = '' if un is None else '[%s]' % (un)
		ax.set_ylabel(var+un)
		ax.legend()
		ax.fill_between(x,y_lq,y_uq,color = 'b' if 'color' not in kwargs else kwargs['color'],alpha=.1,zorder=7)
		if show:
			pp.show()
			pp.pause(10)


if __name__ == '__main__':
	import seaborn as sns
	#Test on a few Storm Sudden Commencement Events
	storm_sudden_commencement = [
	[2006,7,27,13,53.4],
	[2006,8,7,0,35],
	[2006,9,4,0,19],
	[2007,5,7,8,25.2],
	[2007,9,27,11,50.6],
	[2007,10,25,11,35],
	[2007,11,19,18,10.3],
	[2010,4,5,8,26],
	[2010,5,28,2,57.2],
	[2011,4,6,9,33],
	[2011,6,4,20,44],
	[2011,9,9,12,42],
	[2011,9,17,3,43],
	[2011,9,26,12,34.6],
	[2011,10,24,18,31],
	[2011,11,1,9,7.3],
	[2012,2,26,21,39.2],
	[2012,4,23,3,20.2],
	[2012,7,14,18,9],
	[2012,9,3,12,13],
	[2012,10,31,15,38.2],
	[2012,11,12,23,11.4],
	[2013,5,31,16,17.6],
	[2013,10,2,1,54.6],
	[2013,10,8,20,20.4],
	[2014,2,15,13,16.4],
	[2014,6,7,16,52]
	]

	sudden_impulse = [
		[2006,4,28,1,16.8],
		[2006,7,9,21,35.8],
		[2007,7,20,6,17],
		#[2007,12,17,2,53],
		[2008,5,28,2,24.4],
		[2008,11,24,23,51.1],
		[2009,1,30,21,51.7],
		[2009,3,3,6,2.4],
		[2009,4,24,0,53.2],
		[2009,5,28,5,19],
		[2009,9,3,15,51.9],
		#[2011,3,29,16,2],
		[2011,6,17,2,39],
		[2011,7,11,8,50.4],
		#[2011,8,4,21,53.2],
		[2011,10,5,7,36],
		[2011,11,28,21,50.4],
		[2012,5,21,19,36.8],
		[2012,11,23,21,51.8],
		#[2013,1,19,17,32.4],
		[2013,2,16,12,9],
		#[2013,4,13,22,54],
		[2013,6,23,4,25.6],
		[2013,7,9,20,52],
		#[2013,8,20,22,27.8],
		[2013,12,13,13,22.8],
		[2014,1,7,15,12.2],
		[2014,1,9,20,8.2],
		#[2014,2,7,17,5],
		[2014,3,25,20,3],
		[2014,4,19,18,36.4],
		[2014,5,3,17,47],
		[2014,6,23,8,7],
		#[2014,7,14,14,31.2],
		[2014,8,19,6,57.2]
		#[2014,9,6,15,23],
		#[2014,11,1,7,5]
		]
	

	#Set defaults 
	font = {'family' : 'normal',
		'weight' : 'bold',
		'size'   : 14}

	mpl.rc('font',**font)
	mpl.rc('xtick',labelsize=12)
	mpl.rc('ytick',labelsize=12)
	mpl.rc('axes',labelsize=12)
	mpl.rc('axes',labelweight='bold')

	ndays = 5
	sea_ssc = omni_sea(storm_sudden_commencement,cadence='hourly',name='Sudden Commencement',ndays=ndays)
	sea_si = omni_sea(sudden_impulse,cadence='hourly',name='Sudden Impulse',ndays=ndays)

	f = pp.figure(figsize=(11,12))
	#a1 = f.add_subplot(2,1,1)
	a1 = f.add_subplot(6,1,1)
	a2 = f.add_subplot(6,1,2)
	a3 = f.add_subplot(6,1,3)
	a4 = f.add_subplot(6,1,4)
	a5 = f.add_subplot(6,1,5)
	a6 = f.add_subplot(6,1,6)

	pp.ion()
	#sea_ssc.plot_individual(a1,'BZ_GSM')
	sea_si.plot_stats(a1,'ABS_B',color='blue')
	sea_ssc.plot_stats(a1,'ABS_B',color='red')
	a1.xaxis.set_visible(False)

	sea_si.plot_stats(a2,'BZ_GSM',color='blue')
	sea_ssc.plot_stats(a2,'BZ_GSM',color='red')
	a2.xaxis.set_visible(False)

	sea_si.plot_stats(a3,'DST',color='blue')
	sea_ssc.plot_stats(a3,'DST',color='red')
	a3.xaxis.set_visible(False)

	sea_si.plot_stats(a4,'F10_INDEX',color='blue')
	sea_ssc.plot_stats(a4,'F10_INDEX',color='red')
	a4.set_xlabel('Days Since Zero Epoch')

	sea_si.plot_stats(a5,'AE',color='blue')
	sea_ssc.plot_stats(a5,'AE',color='red')
	a5.set_xlabel('Days Since Zero Epoch')	

	sea_si.plot_stats(a6,'PC_N_INDEX',color='blue')
	sea_ssc.plot_stats(a6,'PC_N_INDEX',color='red')
	a6.set_xlabel('Days Since Zero Epoch')	
	
	f.suptitle('Superposed Epoch Analysis of \nSudden Impulse (N=%d) and Sudden Commencement (N=%d)' % (sea_si.nevents,sea_ssc.nevents))

	print str(sea_si.events[0].interval)
	pp.show()
	pp.pause(30)
	f.savefig('/home/liamk/Desktop/omni_sea_si_ssc.png')

