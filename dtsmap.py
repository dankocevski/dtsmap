#!/usr/bin/env python

import os
import subprocess
import sys
import glob
import numpy
import fileinput
import time
import random
from GtApp import GtApp
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pylab as plot
from math import asin, cos, radians, sin, sqrt
import random
import getpass
import traceback
from UnbinnedAnalysis import *


##########################################################################################

def SetPfilesDirectory(pfile_dir):
    """each thread/job which uses FTOOLS must have its own
    PFILES dir"""
    
    try:
        hea_pdir = os.getenv("HEADAS")+"/syspfiles/"
        gt_pdir  = os.getenv("INST_DIR")+"/syspfiles/"
    except:
		print '\n*** Science tools has not be inialized! ***\n'
		print 'Exiting!'
		sys.exit()
	
    
    if(os.path.isdir(pfile_dir)==False):
        print "\nCreating pfile directory:\n%s" % pfile_dir        
        cmd = "mkdir -p "+pfile_dir
        os.system(cmd)

    # --- now copy all the .par files
    cmd = "cp %s/*.par %s"%(hea_pdir,pfile_dir)
    os.system(cmd)
    cmd = "cp %s/*.par %s"%(gt_pdir,pfile_dir)
    os.system(cmd)
    
    # ---Set the new environmental PFILES    
    os.putenv("PFILES",pfile_dir)

    # --- test
 #   print "Testing: temporary PFILES location "
 #   print os.listdir(pfile_dir)
       
    return pfile_dir


##########################################################################################

def AddCandidateSource(ra, dec, xmlModel):

	CandidateSource = """
<source name="CandidateSource" type="PointSource">
    <spectrum type="PowerLaw2">
      <parameter free="1" max=".10000E+07" min=".10000E-09" name="Integral" scale="1" value="1e-4"/>
      <parameter free="1" max="-1.0" min="-5.0" name="Index" scale="1.0" value="-2.1"/>
      <parameter free="0" max="200000.0" min="20.0" name="LowerLimit" scale="1.0" value="100.0"/>
      <parameter free="0" max="300000.0" min="20.0" name="UpperLimit" scale="1.0" value="1e5"/>
    </spectrum>
    <spatialModel type="SkyDirFunction">
      <parameter free="0" max="360." min="-360." name="RA" scale="1.0" value="%.3f"/>
      <parameter free="0" max="90." min="-90." name="DEC" scale="1.0" value="%.3f"/>
    </spatialModel>
  </source>
</source_library>
""" % (float(ra), float(dec))

	print 'Adding a candidate source to the xml model...'
	xmlModelContents = open(xmlModel).read()
	xmlModelContents = xmlModelContents.replace('</source_library>',CandidateSource)
	xmlModelContentsAppended = open(xmlModel,'w')
	xmlModelContentsAppended.write(xmlModelContents)
	xmlModelContentsAppended.close()
	print 'Done.'


##########################################################################################

def ModifySourceModel(xmlModel, newRa, newDec):

	xmlModelModified = xmlModel.replace('.xml','_Modified.xml')

	# Open the files
	infile = open(xmlModel,'r')
	outfile = open(xmlModelModified,'w')

	# Read in the lines from the original flare file
	lines = infile.readlines()
	
	# Loop through each line 
	doFix = False
	
	for line in lines:
		
		newline = line
	
		# Only modify the CandidateSource
		if 'name="CandidateSource"' in newline:
			print 'Modifying Component:'
			print newline
			doFix = True
			
		if 'name="RA"' in newline and doFix == True:
			print 'Changing Ra to %s' % newRa
			newline = '      <parameter free="0" max="360" min="-360" name="RA" scale="1" value="%s" />\n' % newRa

		if 'name="DEC"' in newline and doFix == True:
			print 'Changing Dec to %s' % newDec
			newline = '      <parameter free="0" max="90" min="-90" name="DEC" scale="1" value="%s" />\n' % newDec

		if '</source>' in newline:
			doFix = False

		# Write the new line to disk
		outfile.write(newline)

	print "Writing modified xml file to: %s" % xmlModelModified
	infile.close()
	outfile.close()

	# Replace the original file
	cmd = "cp %s %s; rm %s" % (xmlModelModified, xmlModel, xmlModelModified)
	print cmd
	os.system(cmd)


##########################################################################################	

def great_circle_distance(pnt1, pnt2, radius):
			""" Returns distance on sphere between points given as (latitude, longitude) in degrees. """
			lat1 = radians(pnt1[0])
			lat2 = radians(pnt2[0])
			dLat = lat2 - lat1
			dLon = radians(pnt2[1]) - radians(pnt1[1])
			a = sin(dLat / 2.0) ** 2 + cos(lat1) * cos(lat2) * sin(dLon / 2.0) ** 2
			return 2 * asin(min(1, sqrt(a))) * radius * 57.2957795


##########################################################################################	

def customRALabel(deg):
	return deg

def customDecLabel(deg):
	return deg

##########################################################################################	

def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

 ##########################################################################################	

class dtsmap(object):
	"""A distributed TS map generator.  This tool allows users to manually produce tsmaps by 
	farming each gtlike command to the batch farm at SLAC"""

	def __init__(self):

		self.scfile = None
		self.evfile = None
		self.expmap = None
		self.expcube = None
		self.cmap = None
		self.srcmdl = None
		self.irfs = None
		self.source = None
		self.optimizer = 'NewMinuit'
		self.ftol = 0.1
		self.nxdeg = 5
		self.nydeg = 5
		self.binsz = 0.2
		self.coordsys = 'CEL'
		self.xref = 0
		self.yref = 0
		self.proj = 'AIT'
		self.statistic='UNBINNED'
		self.binNumber = None
		self.outfile = "dtsmap.fits"
		self.outdir = "."
		self.maxJobs = None
		self.batch = True


	def run(self,Test=True, Verbose=False, UseGtlike=False):

		# Convert everything to a float
		ra = float(self.xref)
		dec = float(self.yref)
		nxdeg = float(self.nxdeg)
		nydeg = float(self.nydeg)
		binsize = float(self.binsz)

		# Get the current working directory
		WorkingDirectory = os.getcwd()

		# Define the output, log, and job directories
		if self.outdir == None:
			OutputDirectory = WorkingDirectory
			LogDirectory = "%s/dtsmap/" % WorkingDirectory
			JobsDirectory = "%s/dtsmap/" % WorkingDirectory
		else:
			OutputDirectory = self.outdir
			LogDirectory = self.outdir
			JobsDirectory = self.outdir

		# Define the directory where this script is located
		try:
			ScriptDirectory = os.path.dirname(os.path.realpath(__file__))
		except:
			ScriptDirectory = None

		# Creating the necessary directories
		if(os.path.isdir(OutputDirectory)==False):
			print "\nCreating Directory: " + OutputDirectory
			cmd = "mkdir " + OutputDirectory
			os.system(cmd)
			
		if(os.path.isdir(LogDirectory)==False):
			print "\nCreating Directory: " + LogDirectory
			cmd = "mkdir " + LogDirectory
			os.system(cmd)	

		if(os.path.isdir(JobsDirectory)==False):
			print "\nCreating Directory: " + JobsDirectory
			cmd = "mkdir " + JobsDirectory
			os.system(cmd)

		# Calculate tsmap grid 
		ra_min = ra - self.nxdeg/2.0
		ra_max = ra + self.nxdeg/2.0
		dec_min = dec - self.nydeg/2.0
		dec_max = dec + self.nydeg/2.0

		# Calcualate the range of ra and dec values
		ra_range = numpy.arange(ra_min,ra_max+binsize,binsize)
		dec_range = numpy.arange(dec_min,dec_max+binsize,binsize)
		xsize = len(ra_range)
		ysize = len(dec_range)

		# Make sure that we don't have any ra values below zero or greater than 360, they should wrap ra instead.
		for i in range(len(ra_range)):
			if ra_range[i] < 0:
				ra_range[i] = ra_range[i] + 360.0
			if ra_range[i] > 360:
				ra_range[i] = ra_range[i] - 360.0

		# Make sure that we don't have any dec values below or above +/- 90, they should instead wrap in both ra and dec.
		for i in range(len(dec_range)):
			if dec_range[i] < -90:
				dec_range[i] = ((dec_range[i] + 90) + 90)*-1
				ra_range[i] = ra_range[i] + 180.0
			if dec_range[i] > 90:
				dec_range[i] = 90 - (dec_range[i] - 90)
				ra_range[i] = ra_range[i] + 180


				
		print "\nTS map grid size: %s x %s" % (len(ra_range), len(dec_range))
		print "RA: %s to %s, Dec: %s to %s" % (ra_range[-1], ra_range[0], dec_range[0], dec_range[-1])
		print "Bin size: %s deg" % binsize
		print " "
		print "(%.2f, %.2f)   (%.2f, %.2f)" % (ra_range[-1], dec_range[-1], ra_range[0],  dec_range[-1])
		print " -----------------------------"
		print " |                           |"
		print " |                           |"
		print " |                           |"
		print " |                           |" 
		print " |       (%.2f, %.2f)      |" % (float(ra), float(dec))
		print " |             x             |"
		print " |                           |"
		print " |                           |"
		print " |                           |"
		print " |                           |"					
		print " -----------------------------"
		print "(%.2f, %.2f)    (%.2f, %.2f)" % (ra_range[-1], dec_range[0], ra_range[0],  dec_range[-1])


		print "\nSubmitting a total of %s Jobs" % (xsize * ysize)

		# Mark the start time
		startTime = time.time()
		
		# Loop through all combinations of Ra and Dec and submit a job for each one
		RaDecPairs = {}
		binNumbers = []
		binNumber = 0

		# Keep track of the jobs submitted. This allows the user to specify how many jobs should be running at any given time.
		JobsInQueue = 0

		for raStep in ra_range:
			for decStep in dec_range:

				# Check to see if the number of jobs in the queue is greater than the max.  If so, wait for the jobs to leave the queue before submitting more.
				while JobsInQueue >= int(self.maxJobs):

					print "\nMaximum number of submitted jobs (%s) reached.  Waiting..." % self.maxJobs
					print "Total number of remaining jobs to submit: %s" % remainingJobs

					# Wait 30 seconds before polling the job statuses
					time.sleep(60)

					# Get the number of jobs actively in the queue
					command = "bjobs -g %s | wc" % JobsDirectory	
					process = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
					lines = process.stdout.readlines()
					
					# Extract the number of jobs that are running
					JobsInQueue = int(lines[0].split()[0])
						
				# Setup the job
				logfile = LogDirectory + "/dtsmap_bin%s.log" % binNumber

				# Make the job name
				try:
					sourceNumber = int(self.outdir.split('/')[-1].replace('dtsmap_Source',''))
					jobName = "dts_S%sb%s" % (sourceNumber, binNumber)
				except:
					jobName = "dtsmap_b%s" % binNumber
				
				# Construct the command		
				command = """dtsmap.py scfile=%s evfile=%s expmap=%s expcube=%s cmap=%s srcmdl=%s irfs=%s source=%s optimizer=%s ftol=%s nxdeg=%s nydeg=%s binsz=%s coordsys=%s xref=%s yref=%s proj=%s statistic=%s binNumber=%s outfile=%s outdir=%s PerformAnalysis=1""" % (self.scfile, 
				 			 	self.evfile, self.expmap, self.expcube, self.cmap, self.srcmdl, self.irfs, self.source, self.optimizer, self.ftol, self.nxdeg, 
				 			 	self.nydeg, self.binsz, self.coordsys, raStep, decStep, self.proj, self.statistic, binNumber, self.outfile, OutputDirectory)

				# Specify where to find the python script
				if ScriptDirectory != None:
					command = ScriptDirectory + "/" + command 

				# Construct the process call
				process = 'bsub -o ' + logfile + ' -J ' + jobName + ' -q xlong -R rhel60 -g ' + JobsDirectory + ' "' + command + '"'
			
				# Display the command
				if Verbose == True:
					print process

				# Start the process
				if Test == False:
					if self.batch == True:
						subprocess.call(process, shell=True)
					else:
						os.system(command)

				# Put the ra, dec, and bin info in a dictionary for easy retrieval later
				RaDecPairs[binNumber] = [raStep,decStep]

				# Increment the bin number
				binNumbers.append(binNumber)
				binNumber = binNumber + 1

				# Increment the bin number
				JobsInQueue = JobsInQueue + 1

				# Get the number of remaining jobs
				remainingJobs = ((xsize * ysize) - binNumber)


		print "\nTotal number of jobs submitted: %s" % binNumber
		print "All jobs submitted."
		
		nJobs = -1
		while nJobs != 0:		

			# Determine the number of jobs still running
			command = "bjobs -g %s | wc" % JobsDirectory	
			process = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
			lines = process.stdout.readlines()
		
			# Display the results and the elapsed time since the start of the analysis
			if 'No unfinished job found' in lines[0]:
				nJobs = 0
				print '\nAll Jobs Complete.'

				# Display the elapsed time
				elapsedTimeSeconds = time.time() - startTime
				elapsedTimeMinutes = elapsedTimeSeconds/60.0
				elapsedTimeHours = elapsedTimeMinutes/60.0

				# Modify the units depending on the amount of time that has elapsed
				if elapsedTimeMinutes > 60:				
					print "Total Elapsed Time: %.2f Hours" % elapsedTimeHours
				else:
					print "Total Elapsed Time: %.2f Minutes" % elapsedTimeMinutes

				# Make sure we actually break out of this while loop!
				break

			else:
				
				# Get the number of jobs in the queue from the stdout
				nJobs = int(lines[0].split()[0])
				print "\nJobs Remaining: %s" % nJobs
				
				# Display the elapsed time
				elapsedTimeSeconds = time.time() - startTime
				elapsedTimeMinutes = elapsedTimeSeconds/60.0
				elapsedTimeHours = elapsedTimeMinutes/60.0

				if elapsedTimeMinutes > 60:				
					print "Elapsed Time: %.2f Hours" % elapsedTimeHours
				else:
					print "Elapsed Time: %.2f Minutes" % elapsedTimeMinutes 

			# Wait 120 seconds before polling the job statuses
			time.sleep(120)

		if UseGtlike == True:
			# Gather the results (gtlike)
			command = "grep 'TS value' %s/likelihoodResults_bin*.txt" % LogDirectory	
			process = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
			lines = process.stdout.readlines()
		else:
			# Gather the results (pylikelihood)
			command = "grep 'TS = ' %s/dtsmap_bin*.log" % LogDirectory	
			process = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
			lines = process.stdout.readlines()

		# Read in the results
		binNumbersReadIn = []
		TSs = []
		for line in lines:
			if ('No match' in line) or ('grep' in line):
				print "Error: *** No available likelihood results ***"
				print "Check your science tools setup and try again..."
				print "Exiting.\n"

				Results = {'MaxTS':None, 'MaxRA':None, 'MaxDec':None, 'Error':None, 'Index':None, 'IndexError':None, 'Flux':None,'FluxError':None, 'MaxBin':None}
				return Results
			else:

				if UseGtlike == True:

					# For use with gtlike results
					binNumber = int(line.split()[0].split('/')[-1].replace('.txt:\'TS','').replace('likelihoodResults_bin','').strip())
					TS = float("%.2f" % float(line.split()[2].strip().replace("'","").replace(',','')))
					binNumbersReadIn.append(binNumber)
					TSs.append(TS)

				else:

					# For use with pylikelihood results
					binNumber = int(line.split()[0].split('/')[-1].replace('.log:TS','').replace('dtsmap_bin','').strip())
					TS = float("%.2f" % float(line.split()[2].strip().replace("'","").replace(',','')))
					binNumbersReadIn.append(binNumber)
					TSs.append(TS)


		# Put the results in a dictionary for easy retrieval later
		LikelihoodResults = {key:value for key, value in zip(binNumbersReadIn,TSs)}

		# Make an matrix to store the values.  Saving any missing values as NaNs
		TSMap = numpy.zeros(shape=(xsize, ysize))
		binNumber = 0
		xStep = 0
		yStep = 0
		for raStep in ra_range:
			yStep = 0
			for decStep in dec_range:
				if binNumber in LikelihoodResults:
					TSMap[xStep][yStep] = LikelihoodResults[binNumber]
					pass
				else:
					print "Bad bin: %s" % binNumber
					TSMap[xStep][yStep] = numpy.nan
					pass
				if TSMap[xStep][yStep] < -1:
					TSMap[xStep][yStep] = numpy.nan
					pass
				binNumber = binNumber + 1
				yStep = yStep + 1
				pass
			xStep = xStep + 1
			pass


		# # Loop through the results matrix and fill any missing values with interpolations from nearby pixels
		# binNumber = 0
		# xStep = 0
		# yStep = 0
		# for raStep in ra_range:
		# 	yStep = 0
		# 	for decStep in dec_range:
		# 		if numpy.isnan(TSMap[xStep][yStep]) == True:

		# #			try:
		# #				TSMap[xStep][yStep] = numpy.mean([TSMap[xStep-1,yStep-1],TSMap[xStep-1,yStep],TSMap[xStep-1,yStep+1],TSMap[xStep,yStep-1],TSMap[xStep,yStep+1],TSMap[xStep+1,yStep-1],TSMap[xStep+1,yStep],TSMap[xStep+1,yStep+1]])
		# 			try:					
		# 				if numpy.isnan(TSMap[xStep-1][yStep]) == False and  numpy.isnan(TSMap[xStep+1][yStep]) == False:
		# 					TSMap[xStep][yStep] = (TSMap[xStep-1][yStep] + TSMap[xStep+1][yStep]) / 2.0
		# 				elif numpy.isnan(TSMap[xStep][yStep-1]) == False and  numpy.isnan(TSMap[xStep][yStep+1]) == False:
		# 					TSMap[xStep][yStep] = (TSMap[xStep][yStep-1] + TSMap[xStep][yStep+1]) / 2.0
		# 				else:
		# 					TSMap[xStep][yStep] = numpy.nan
		# 			except:
		# 				TSMap[xStep][yStep] = numpy.nan

		# 		yStep = yStep + 1
		# 		pass
		# 	xStep = xStep + 1

		# Finding the maximum TS
		MaxTS = numpy.nanmax(TSMap)
		MaxBin = numpy.nanargmax(TSMap)

		# Check to see if a maximum couldn't be found.
		if numpy.isnan(MaxBin) == True:
			print "\nAnalysis Complete."
			print "Maximum TS: %s" % 'None'
			print "Coordinates: RA = %s, Dec = %s" % ('NA', 'NA')
			print "*** Unable to locate maximum TS ***"

			print '\nCleaning up...'
			# Cat the gtlike log files
			cmd = "cat %s/likelihoodResults_bin*.txt > dtsmap_LikelihoodResults.txt" % self.outdir
#			print cmd
			os.system(cmd)

			# Cat the log files
			cmd = "cat %s/dtsmap_bin*.log > dtsmap_LikelihoodFits.log" % self.outdir
#			print cmd
			os.system(cmd)

			# Erase the individual xml files
			cmd = "rm %s/ModelSource_bin*.xml" % self.outdir
#			print cmd
			os.system(cmd)

			print 'Done.'

			Results = {'MaxTS':None, 'MaxRA':None, 'MaxDec':None, 'Error':None, 'Index':None, 'IndexError':None, 'Flux':None,'FluxError':None, 'MaxBin':None}
			return Results

		else:

			MaxRa = RaDecPairs[MaxBin][0]
			MaxDec = RaDecPairs[MaxBin][1]

		# Define default spectral parameters
		index = 'NA'
		indexError = 'NA'
		flux = 'NA'
		fluxError = 'NA'

		if UseGtlike == True:

			# Extract the fit parameters for the bin with the maximum TS (gtlike)
			MaxBinFile = "%s/dtsmap_bin%s.log" % (LogDirectory, MaxBin)
			for line in fileinput.input([MaxBinFile]):
				if 'Index:' in line: 
					lineContents = line.split()	
					index = lineContents[1]
					indexError = lineContents[3]
				if 'Flux:' in line: 
					lineContents = line.split()	
					flux = lineContents[1]
					fluxError = lineContents[3]
					break
			fileinput.close()

		else:

			# Extract the spectral fit parameters for the bin with the maximum TS (pyLikelihood)
			MaxBinFile = "%s/dtsmap_bin%s.log" % (LogDirectory, MaxBin)
			for line in fileinput.input([MaxBinFile]):
				if 'Flux =' in line: 
					lineContents = line.split()	
					flux = float(lineContents[2])
					fluxError = float(lineContents[4])
				if 'Index =' in line: 
					lineContents = line.split()	
					index = float(lineContents[2])
					indexError = float(lineContents[4])
			fileinput.close()


		# Rotate and flip the matrix in order to it to match ra and dec ordering conventions
		TSMapRotated = numpy.rot90(TSMap)
		TSMapFlippedUp = numpy.flipud(TSMapRotated)
		TSMapFlippedLeft2Right = numpy.fliplr(TSMapFlippedUp)
		MaxRa = RaDecPairs[numpy.nanargmax(TSMap)][0]
		MaxDec = RaDecPairs[numpy.nanargmax(TSMap)][1]

		# Import the basemap module
		sys.path.append("/nfs/slac/g/ki/ki08/kocevski/LATBA/lib/python_rhel6-64/")
		from mpl_toolkits.basemap import Basemap

		# Create the figure
		fig = plot.figure()
		ax = fig.add_subplot(111)

		# Create a base map on which to plot the results
		#m = Basemap(llcrnrlon=ra_range[-1], llcrnrlat=dec_range[0], urcrnrlon=ra_range[0], urcrnrlat=dec_range[-1], projection='laea', lon_0 = ra, lat_0 = dec,resolution ='l',area_thresh=1000.)# ,celestial=True )
		m = Basemap(height=5.5e5,width=5.5e5, projection='laea', lon_0 = ra*-1, lat_0 = dec, resolution ='l', area_thresh=1000., celestial=True )

		# Set the plot limits (in map coordinates)
		xMin, yMin = m(ra_range[0], dec_range[0])
		xMax, yMax = m(ra_range[-1], dec_range[-1])
		m.lonmin = xMin
		m.lonmax = xMax
		m.latmin = yMin
		m.latmax = yMax

		# Plot the matrix as an image
		extent=[xMax, xMin, yMin, yMax]
		#m.imshow(TSMapFlippedLeft2Right, origin='lower', extent=extent)
		m.imshow(TSMapFlippedLeft2Right, origin='lower')

		# Setup the map grid
		m.drawmapboundary(fill_color='#ffffff')
		m.drawparallels(numpy.arange(181)-90,labels=[1,0,0,0], fmt=customDecLabel, color='grey', linewidth=0)
		m.drawmeridians(numpy.arange(0,360,1),labels=[0,0,0,1], fmt=customRALabel, color='grey', linewidth=0)
		m.suppress_ticks = False
		m.fix_aspect = False
	
		# Force the aspect ratio to be 1:1
		try:
			forceAspect(ax,aspect=1)
			#extent=[xMax, xMin, yMin, yMax]
			#ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2])))
		except Exception, message:
			print traceback.format_exc()

		# Add a color bar
		colorBarObject = plot.colorbar()
		colorBarObject.set_label('TS')

		# Setup the plot
		plot.xlabel('RA')
		plot.ylabel('Dec')
		plot.gca().xaxis.labelpad = 20
		plot.gca().yaxis.labelpad = 20

		# Get the ra and dec of the max TS in map coordinates
		mMaxRa, mMaxDec = m(MaxRa, MaxDec)
		mRa, mDec = m(float(ra), float(dec))

		# Define the 68%, 90%, and 95% confidence levels
		MaxTS_68CL = (MaxTS-2.30)
		MaxTS_90CL = (MaxTS-4.61)
		MaxTS_95CL = (MaxTS-5.99)

		# Add contour lines
		try:
			contourObject = plot.contour(TSMapFlippedLeft2Right, levels=[MaxTS_68CL,MaxTS_90CL,MaxTS_95CL], origin='lower', extent=[ra_range[-1], ra_range[0], dec_range[0], dec_range[-1]])

			# Extract the 68% contour interval data
			contourLine68 = numpy.array(contourObject.allsegs[0])[0]
			contourLine68_ra = contourLine68[:,0]
			contourLine68_dec = contourLine68[:,1]

			# Extract the 90% contour interval data	
			contourLine90 = numpy.array(contourObject.allsegs[1])[0]
			contourLine90_ra = contourLine90[:,0]
			contourLine90_dec = contourLine90[:,1]

			# Extract the 95% contour interval data	
			contourLine95 = numpy.array(contourObject.allsegs[2])[0]
			contourLine95_ra = contourLine95[:,0]
			contourLine95_dec = contourLine95[:,1]		

			# Convert the contour line data into map coordinates
			mContourLine68_ra, mContourLine68_dec = m(contourLine68_ra, contourLine68_dec)
			mcontourLine90_ra, mcontourLine90_dec = m(contourLine90_ra, contourLine90_dec)
			mcontourLine95_ra, mcontourLine95_dec = m(contourLine95_ra, contourLine95_dec)

			# Plot the contour lines in map coordinates
			plot.plot(mContourLine68_ra, mContourLine68_dec, linestyle='dotted', color='black')	
			plot.plot(mcontourLine90_ra, mcontourLine90_dec, linestyle='dotted', color='black')	
			plot.plot(mcontourLine95_ra, mcontourLine95_dec, linestyle='dotted', color='black')	

			# Find the spherical distance from the contour interval data to the starting ra and dec
			sphericalDistances68 = []
			for contour_ra, contour_dec in zip(contourLine68_ra,contourLine68_dec):
				sphericalDistance = great_circle_distance([MaxRa,MaxDec],[contour_ra,contour_dec],1)
				sphericalDistances68.append(sphericalDistance)
				pass

			sphericalDistances90 = []
			for contour_ra, contour_dec in zip(contourLine90_ra,contourLine90_dec):
				sphericalDistance = great_circle_distance([MaxRa,MaxDec],[contour_ra,contour_dec],1)
				sphericalDistances90.append(sphericalDistance)
				pass

			sphericalDistances95 = []
			for contour_ra, contour_dec in zip(contourLine95_ra,contourLine95_dec):
				sphericalDistance = great_circle_distance([MaxRa,MaxDec],[contour_ra,contour_dec],1)
				sphericalDistances95.append(sphericalDistance)
				pass		

			# Find the median of the spherical distances.  This will be considered as the error in the localization
			sphericalDistances68 = numpy.array(sphericalDistances68)
			errorRadius68 = numpy.median(sphericalDistances68)
			sphericalDistances90 = numpy.array(sphericalDistances90)
			errorRadius90 = numpy.median(sphericalDistances90)
			sphericalDistances95 = numpy.array(sphericalDistances95)
			errorRadius95 = numpy.median(sphericalDistances95)

			# Annotate the plot
			m.scatter(mMaxRa, mMaxDec, marker='x', s=75, facecolors='none', edgecolors='w')
			m.scatter(mRa, mDec, marker='+', s=75, facecolors='none', edgecolors='w')

			# Add a plot legend
			plot.annotate("Max TS = %s\nRA = %s, Dec = %s\nError = +/-%0.3f (95%% CL)" % (MaxTS, MaxRa, MaxDec, errorRadius95), xy=(0,0), xycoords='axes points', xytext=(10,10), textcoords='offset points', ha='left', va='bottom', bbox=dict(boxstyle='round,pad=0.2', fc='w', alpha=0.3))#, arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=1e-10',color='b'))

			# Print the results
			print "\nAnalysis Complete."
			print "Maximum TS: %s" % MaxTS
			print "Coordinates: RA = %s, Dec = %s" % (MaxRa, MaxDec)
		#	print "Error Radius (68%%): %.4f deg" % errorRadius68
			print "Error Radius (90%%): %.4f deg" % errorRadius90
			print "Error Radius (95%%): %.4f deg" % errorRadius95

			print ''
			print 'Best Fit Parameters:'
			print "Index = %.2f +/- %.2f" % (index, indexError)
			print "Flux = %.2e +/- %.2e" % (flux, fluxError)

		except:

			# Set the error radius to None
			errorRadius95 = None

			# Annotate the plot
			m.scatter(mMaxRa, mMaxDec, marker='x', s=75, facecolors='none', edgecolors='w')
			m.scatter(mRa, mDec, marker='+', s=75, facecolors='none', edgecolors='w')

			# Add a plot legend
			plot.annotate("Max TS = %s\nRA=%s, Dec=%s\nError = NA" % (MaxTS, MaxRa, MaxDec), xy=(0,0), xycoords='axes points', xytext=(10,10), textcoords='offset points', ha='left', va='bottom', bbox=dict(boxstyle='round,pad=0.2', fc='w', alpha=0.3))#, arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=1e-10',color='b'))

			# Print the results
			print "\nAnalysis Complete."
			print "Maximum TS: %s" % MaxTS
			print "Coordinates: RA = %s, Dec = %s" % (MaxRa, MaxDec)
		#	print "Error Radius (68%%): NA"
			print "Error Radius (90%%): NA"
			print "Error Radius (95%%): NA"

			print ''
			print 'Best Fit Parameters:'
			print "Index = %s +/- %s" % (index, indexError)
			print "Flux = %s +/- %s" % (flux, fluxError)			


	#	plot.show()
		print "\nSaving ts map: %s" % self.outfile
		plot.savefig(self.outfile, bbox_inches='tight', dpi=100)
		plot.close()

		# Saving the max bin source model
		print "Saving optimized xml model: %s/ModelSource_MaxTS.xml" % (self.outdir)
		cmd = "cp %s/ModelSource_bin%s.xml %s/ModelSource_MaxTS.xml" % (self.outdir, MaxBin, self.outdir)
		os.system(cmd)

		# Saving the max bin log file
		print "Saving optimized dtsmap log: %s/dtsmap_MaxTS.log" % (self.outdir)
		cmd = "cp %s/dtsmap_bin%s.log %s/dtsmap_MaxTS.log" % (self.outdir, MaxBin, self.outdir)
		os.system(cmd)

		if UseGtlike == True:

			# Saving the max bin results file (for use with gtlike)
			print "Saving optimized likelihood results: %s/likelihoodResults_MaxTS.txt" % (self.outdir)
			cmd = "cp %s/likelihoodResults_bin%s.txt %s/likelihoodResults_MaxTS.txt" % (self.outdir, MaxBin, self.outdir)
			os.system(cmd)

			# Cat the gtlike log files (for use with gtlike)
			cmd = "cat %s/likelihoodResults_bin*.txt > dtsmap_likelihoodResults.log" % self.outdir
			print cmd
			os.system(cmd)

		else:

			# Write the max bin results to a file
			likelihoodResults = "%s/likelihoodResults_MaxTS.txt" % (self.outdir)
			output = open(likelihoodResults, 'w')
			output.write("Maximum TS Likelihood Results\n")
			output.write("MaxTS: %.2f\n" % MaxTS)
			if errorRadius95 == None:
				output.write("RaDec: %.3f %.3f +/- NA\n" % (MaxRa, MaxDec))
			else:
				output.write("RaDec: %.3f %.3f +/- %.3f\n" % (MaxRa, MaxDec, errorRadius95))
			output.write("PhotonIndex: %.2f +/- %.2f\n" % (index, indexError))
			output.write("PhotonFlux: %.2e +/- %.2e\n" % (flux,fluxError))
			output.close()


		# Cat the log files
		# cmd = "cat %s/dtsmap_bin*.log > %s/dtsmap.log" % (self.outdir, self.outdir)
		# print cmd
		# os.system(cmd)

		# Remove the individual log files
		# cmd = "rm %s/dtsmap_bin*.log" % (self.outdir)
		# print cmd
		# os.system(cmd)		

		# Erase the individual xml files
		# cmd = "rm %s/ModelSource_bin*.xml" % self.outdir
		# print cmd
		# os.system(cmd)

		print 'Done.'

		# Save some selected information in a dictionary
		Results = {'MaxTS':MaxTS, 'MaxRA':MaxRa, 'MaxDec':MaxDec, 'Error':errorRadius95, 'Index':index, 'IndexError':indexError, 'Flux':flux,'FluxError':fluxError, 'MaxBin':MaxBin}

		return Results


##########################################################################################

	def PerformLikelihoodAnalysis(self):

		print "\nPerforming likelihood analysis on position: ra=%s, dec=%s" % (self.xref, self.yref)

		# Wait a random amount of time between 1 and 5 minutes before starting in order to not crash the asf/nsf disks at SLAC
		waitTime = random.random()*300
		time.sleep(waitTime)

		# Defind the scratch directory
		JobID = os.environ.get('LSB_JOBID')
		Username = getpass.getuser()
		ScratchDirectory = "/scratch/%s/%s/" % (Username, JobID)

		# Define the pfile directory
		if JobID == 'None':
			PFILESDirectory = "%s/pfiles_%s/" % (self.outdir, self.binNumber)	
		else:
			PFILESDirectory = "%s/pfiles/" % ScratchDirectory

		# Create the output directory if it doesn't already exist
	  	if(os.path.isdir(self.outdir)==False):
	  		print "\n >> Creating Directory: " + self.outdir
	  		cmd = "mkdir " + self.outdir
	  		os.system(cmd)

		# Define where to save the results
		likelihoodResults = '%s/likelihoodResults_bin%s.txt' % (self.outdir, self.binNumber)

		# Remove any pre-existing pfiles
		if(os.path.isdir(PFILESDirectory)==True):
	  		cmd = "rm -r %s" % PFILESDirectory
	  		os.system(cmd)			

		# Set the new pfiles directory
		SetPfilesDirectory(PFILESDirectory)

		# Make a copy of the source model
		xmlModelWithPutativeSource = '%s/ModelSource_bin%s.xml' % (self.outdir, self.binNumber)
		cmd = "cp " + self.srcmdl + " " + xmlModelWithPutativeSource
		print cmd
	  	os.system(cmd)

		# Add a putative point source at the requested location
#		AddCandidateSource(self.xref, self.yref, xmlModelWithPutativeSource)
		ModifySourceModel(xmlModelWithPutativeSource, self.xref, self.yref)

		# # Import the necessary gtapps	
		# gtlike = GtApp('gtlike')

		# # Run the likelihood analysis
		# print '\nPerforming the likelihood fit:'		
		# gtlike.run(statistic=self.statistic,
		# 				scfile=self.scfile,
		# 				evfile=self.evfile,
		# 				expmap=self.expmap,
		# 				expcube=self.expcube,
		# 				srcmdl=xmlModelWithPutativeSource,
		# 				irfs=self.irfs,
		# 				optimizer=self.optimizer,
		# 				results=likelihoodResults,
		# 				plot='no',
		# 				save='yes')

		# Setup the unbinned likelihood object
		print '\nPerforming the likelihood fit:'
		try:
		
			obs = UnbinnedObs(self.evfile,self.scfile,expMap=self.expmap,expCube=self.expcube,irfs=self.irfs)
		
			# Define the likelihood object
			#like = UnbinnedAnalysis(obs,xmlModelWithPutativeSource,optimizer=self.optimizer)
			like = UnbinnedAnalysis(obs,xmlModelWithPutativeSource,optimizer='MINUIT')
				
			# Setup the likelihood parameters
			Source = 'CandidateSource'
			Integral = like.par_index(Source, 'Integral')
			Index = like.par_index(Source, 'Index')
			LowerLimit = like.par_index(Source, 'LowerLimit')
			UpperLimit = like.par_index(Source, 'UpperLimit')
		
			# Setup the likelihood bounds
			like[Integral].setScale(1e-3)
			like[Index].setBounds(-5, -0.5)
			# like[LowerLimit] = emin
			# like[UpperLimit] = emax
		
			# Perform the likelihood fit
			#optObject = pyLike.NewMinuit(like.logLike)				  
			#like.fit(verbosity=0,covar=True,tol=0.02,optObject=optObject)
			like.fit(verbosity=1,covar=True,tol=1e-10,optimizer='MINUIT', optObject=None)
	
			# Extract the best fit index
			IndexValue = like[Index].value()
			IndexError = like[Index].error()

			# Extract the best fit flux
			FluxValue = like.flux(Source, emin=100, emax=3e5)
			FluxError = like.fluxError(Source, emin=100, emax=3e5)

			# Extract likelihood fit results
			print '\nLikelihood Results:'
			print like.model[Source]
			print "TS = %s" % like.Ts(Source)
			print "Flux = %s +/- %s" % (FluxValue, FluxError)
			print "Index = %s +/- %s" % (IndexValue, IndexError)

			# Save the xml file
			like.writeXml(xmlFile=xmlModelWithPutativeSource)
			
		except Exception, message:
			print traceback.format_exc()	


		# Cleaning up
		if(os.path.isdir(ScratchDirectory)==True):
			cmd = "rm -R %s" % ScratchDirectory
			print ""
			print cmd
	  		os.system(cmd)
			


##########################################################################################

if __name__ == '__main__':

	errorMessage = """
	Distributed TS Map v1.0 
	Support Contact: Daniel Kocevski (daniel.kocevski@nasa.gov)

	Usage:  dtsmap.py scfile=scfile, evfile=evfile, expmap=expmap, expcube=expcube, cmap=cmap, srcmdl=srcmdl, irfs=irfs, source=source, 
	optimizer=optimizer, ftol=ftol, nxdeg=xdeg, nydeg=nydeg, binsz=binsz, coordsys=coordsys, xref=xref, yref=yref, proj=proj, outfile=outfile, outdir=outdir

	scfile 		->	Spacecraft FT2 file
	evfile 		->	Fermi-LAT FT1 file
	expmap 		->	Exposure Map
	expcube 	->	Live time cube
	cmap 		->	Counts map
	srcmdl 		->	Likelihood xml model
	irfs 		->	instrument response function
	irfs 		->	instrument response function
	source 		->	Name of the source of interest
	optimizer	->	Likelihood fit optimizer (e.g. NewMinuit)
	ftol		->	Likelihood fit tolerance
	nxdeg		->	Number of degres the TS map will space in the x (RA) direction
	nydeg		->	Number of degres the TS map will space in the y (Dec) direction
	binsz		->	TS map resolution in degrees
	coordsys	->	Coordinate system (e.g. CEL)
	xref		->	Center RA of the TS map
	xref		->	enter Dec of the TS map
	proj		->	TS map projection (e.g. AIT)
	outfile		->	Output TS map filename
	outdir		->	Writable directory in which to store intermediate products


	"""

	if len(sys.argv) > 1:

		# Extact the keywords
		kwargs = {}
		for keyword in sys.argv:
			if '=' in keyword:
				key, value = keyword.split('=', 1)
				kwargs[key] = value

		# Perform the specified action
		if 'PerformAnalysis' not in kwargs:

			dtsmapObject = dtsmap()
			dtsmapObject.scfile = kwargs['scfile']
			dtsmapObject.evfile = kwargs['evfile']
			dtsmapObject.expmap = kwargs['expmap']
			dtsmapObject.expcube = kwargs['expcube']
			dtsmapObject.cmap = kwargs['cmap']
			dtsmapObject.srcmdl = kwargs['srcmdl']
			dtsmapObject.irfs = kwargs['irfs']
			dtsmapObject.source = kwargs['source']
			dtsmapObject.optimizer = kwargs['optimizer']
			dtsmapObject.ftol = kwargs['ftol']
			dtsmapObject.nxdeg = kwargs['nxdeg']
			dtsmapObject.nydeg = kwargs['nydeg']
			dtsmapObject.binsz = kwargs['binsz']
			dtsmapObject.coordsys = kwargs['coordsys']
			dtsmapObject.xref = kwargs['xref']
			dtsmapObject.yref = kwargs['yref']
			dtsmapObject.proj = kwargs['proj']
			dtsmapObject.outfile = kwargs['outfile']
			dtsmapObject.outdir = kwargs['outdir']
#			dtsmapObject.binNumber = kwargs['binNumber']
		
			# Submit the individual jobs
			dtsmapObject.run(Test=True)
		
		else:

			dtsmapObject = dtsmap()
			dtsmapObject.scfile = kwargs['scfile']
			dtsmapObject.evfile = kwargs['evfile']
			dtsmapObject.expmap = kwargs['expmap']
			dtsmapObject.expcube = kwargs['expcube']
			dtsmapObject.cmap = kwargs['cmap']
			dtsmapObject.srcmdl = kwargs['srcmdl']
			dtsmapObject.irfs = kwargs['irfs']
			dtsmapObject.source = kwargs['source']
			dtsmapObject.optimizer = kwargs['optimizer']
			dtsmapObject.ftol = kwargs['ftol']
			dtsmapObject.nxdeg = kwargs['nxdeg']
			dtsmapObject.nydeg = kwargs['nydeg']
			dtsmapObject.binsz = kwargs['binsz']
			dtsmapObject.coordsys = kwargs['coordsys']
			dtsmapObject.xref = kwargs['xref']
			dtsmapObject.yref = kwargs['yref']
			dtsmapObject.proj = kwargs['proj']
			dtsmapObject.outfile = kwargs['outfile']
			dtsmapObject.outdir = kwargs['outdir']
			dtsmapObject.binNumber = kwargs['binNumber']

			# Run the likelihood analysis
			dtsmapObject.PerformLikelihoodAnalysis()

	else:	

 		print errorMessage

	 	sys.exit()
