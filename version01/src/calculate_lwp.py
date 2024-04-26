#This program calculates liquid water path over the domain
#from WRF output.
#Written by Chris Phillips
#
#Requirements
#python 3+
#NetCDF4
#numpy

#Importing modules
import netCDF4 as nc
import numpy

#Defining the main function
#Inputs:
#filein, string, the input file path
#Outputs, float array, numpy array of liquid water path. Dimensions are nt,nx,ny
#	were nx,ny are the WRF grid dimensions and nt is the number of time steps
def calculate_lwp(filein):
	
	#Defining constants
	Rd = 287.05 #Dry air gas constant [J/kg/K]
	Rv = 461.5 #Water vapor  gas constant [J/kg/K]
	cp = 1005 #Specicif heat of air at constant pressure [J/kg/K]
	cv = 718 #Specific heat of air at constant volume [J/kg/K]
	pres0 = 100000 #Reference pressure [Pa]

	#Reading the input files
	fn = nc.Dataset(filein, "r")
	qc = numpy.array(fn.variables["QCLOUD"][:,:,:,:])    #Cloud water mixing ratio [kg/kg]
	qv = numpy.array(fn.variables["QVAPOR"][:,:,:,:])    #Water vapor mixing ratio [kg/kg]
	pres = numpy.array(fn.variables["PB"][:,:,:,:]) #Pressure [Pa]
	theta = numpy.array(fn.variables["T"][:,:,:,:])+numpy.array(fn.variables["T00"][0]) #Potential Temperature [K]
	zh = numpy.array(fn.variables["PHB"][:,:,:,:])   #Height on model levels [m]
	fn.close()

	#Determing domain dimensions
	nt = len(qc[:,0,0,0])
	nz = len(qc[0,:,0,0])
	ny = len(qc[0,0,:,0])
	nx = len(qc[0,0,0,:])

	#Calculating air temperature at all grid points
	temp = theta*(pres/pres0)**(Rd/cp)

	#Calculating density of dry air at all grid points
	rho = pres/Rd/temp

	#Calculating liquid water path
	lwp = numpy.zeros((nt,ny,nx))
	for k in range(nz-1):
		lwp += (rho[:,k,:,:]+rho[:,k+1,:,:])/2*(qc[:,k,:,:]+qc[:,k+1,:,:])/2*(zh[:,k+1,:,:]-zh[:,k,:,:])
	
	return numpy.squeeze(lwp) #Removing extra dimensions
	
