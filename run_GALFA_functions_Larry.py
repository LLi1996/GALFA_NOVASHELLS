"""
This Python module is for the extraction of spectrum, integrated intensity and all kinds of maps from a standard GALFA_HI cube.
This module will 


This module is written by Larry with input from Yong.
08/04/2015
"""

import os
import sys
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import ConnectionPatch
import astropy.units as u #currently not used

### CUBE AND FOLDER STUFF ####################################
def GALFAHI_cubeinfo(fitsfile):
    """
    this function reads and extract information form a GALFA HI data cube
    and returns the data, header, velocity ticks, DEC ticks, RA ticks and velocity channel separation 
    # the majority of this funciton was taking from Yong
    # note I made the variables global for future use in this modeule
    """
    global data
    global vlsr 
    global delv
    global raSize
    global decSize
    global minRA
    global maxRA
    global minDEC
    global maxDEC
    global delta

    data = fits.getdata(fitsfile)
    header = fits.getheader(fitsfile)
    
    # RA and DEC are in degrees!
    ra = header['CRVAL1'] + header['CDELT1'] * (np.arange(data.shape[2])+1 - header['CRPIX1'])             
    dec = header['CRVAL2'] + header['CDELT2'] * (np.arange(data.shape[1])+1 - header['CRPIX2'])  
    # vlsr and delv in km/s
    vlsr = (header['CRVAL3'] + header['CDELT3'] * (np.arange(data.shape[0])+1 - header['CRPIX3']))*(10**(-3))
    delv = header['CDELT3']*(10**(-3))           

    # this is all getting properties of the various fields
    raSize = ra.size
    decSize = dec.size
    minRA = ra[raSize-1]
    maxRA = ra[0]
    minDEC = dec[0]
    maxDEC = dec[decSize-1]
    delta = 1.0/60
    # essentially the converion between degrees and arcMin, since the pixels are all 1 arcMin
    # note 1/60 is slightly lower (by 3.3*10^-8 degrees) than the actual separation in both RA and DEC
    # ^which is fine because over the entire cube, 512 pixels, the difference is only around 1.7*10^-5 degrees
    
    return data, header, vlsr, dec, ra, delv

def getTargetCubeinfo(objectRA, objectDEC):
    """
    this function finds the (neariest) cube which the target belongs to
    and returns the RA & DEC of that cube
    # note this takes in the object's coordinates, RA & DEC, in degrees
    """
    cubeSep = 8
    cubeRAstep = 4
    cubeDECstep = 2.35
    
    # calculates the RA and DEC values of the GALFA cube
    this_cube_ra = ((objectRA-cubeRAstep) // 8) * 8 + cubeRAstep
    this_cube_dec = ((objectDEC-cubeDECstep) // 8) * 8 + cubeDECstep
    if ((objectRA-cubeRAstep) % cubeSep > 4):
        this_cube_ra += 8
    if ((objectDEC-cubeDECstep) % cubeSep > 4):
        this_cube_dec += 8
        
    return this_cube_ra, this_cube_dec

def setTargetInfo(name='', raCoord=0, decCoord=0, gamma=0, gammaE=0, c=''):
    """
    this functions adds a few global variables to this module
    """
    global objectName
    global objectRA
    global objectDEC
    global objectGamma
    global objectGammaError
    global cycle

    objectName = name
    objectRA = raCoord
    objectDEC = decCoord
    objectGamma = gamma
    objectGammaError = gammaE
    cycle = c

def checkTargetFolder(objectName, cycleNum):
    """
    this function checks to see if there is a target directory,
    if not, this function creates the target directory
    """
    dirName = os.path.dirname('Results/' + objectName + '/' + cycleNum + '/')
    if not os.path.exists(dirName):
        os.makedirs(dirName)


### PLOTTING SETTINGS ####################################
def setVelocityRange(minVelocity, maxVelocity):
    """
    this function sets the plotting velocity range
    and returns the velocity indices of the input range
    # takes in the desired minimum and maximum velocity"
    """
    minVindex = 1024 - (0 - minVelocity)//delv
    maxVindex = 1024 + (maxVelocity)//delv

    return minVindex, maxVindex

def degreeCoordinatesToCubeIndices(inRA, inDEC):
    """
    this function takes in the RA&DEC coordinates in degrees
    and returns the RA&DEC coordinates in GALFA cube indices (for that specific cube)
    """
    raIndex = (raSize-1) - (inRA - minRA)//delta
    decIndex = (inDEC - minDEC)//delta

    return raIndex, decIndex

def setPlotRange(centerRA, centerDEC, raRange, decRange):
    """
    this function sets the plotting ranges
    and returns the RA & DEC and RA & DEC indices of the input RA & DEC range
    # taks in the center and the xy (RA & DEC) range
    # note graphs below are usually centered on the object, AKA objectRA & objectDEC = centerRA & centerDEC
    # note raRange and decRange are expected in pixels (1 pixel = 1 arcMin in our case)
    # note minRA, minDEC, raSize & delta (used below) are properties of the cube
    """
    # calculates the range of RA and DEC that will be included
    minGraphRA = centerRA - (raRange / 2) * delta
    maxGraphRA = centerRA + (raRange / 2) * delta
    minGraphDEC = centerDEC - (decRange / 2) * delta
    maxGraphDEC = centerDEC + (decRange / 2) * delta
    
    # convert the above into indexes
    minRAindex, minDECindex = degreeCoordinatesToCubeIndices(minGraphRA, minGraphDEC)
    maxRAindex, maxDECindex = degreeCoordinatesToCubeIndices(maxGraphRA, maxGraphDEC)
    
    return maxGraphRA, minGraphRA, minGraphDEC, maxGraphDEC, maxRAindex, minRAindex, minDECindex, maxDECindex


### SPECTRUM EXTRACTION ####################################
def getSpectrum(objectRA, objectDEC, boxSize=0, box = 'no', maxRAindex=0, minRAindex=0, minDECindex=0, maxDECindex=0, region = 'no'):
    """
    this function extracts the full (~~~ -700 to 700 km/s) spectrum for a specific point or a region
    and returns that spectrum
    # note this by default returns the spectrum of a point
    # note boxSize is expected in pixels (and 1 pixel = 1 arcMin for us)
    # note that, if box or region == 'yes', output spectrum is normalized by the number of pixels in the region
    """
    if box == 'yes':
        plotRADECindices = setPlotRange(objectRA, objectDEC, boxSize, boxSize)
        spectrum = data[:, 
                        plotRADECindices[6]:plotRADECindices[7], 
                        plotRADECindices[4]:plotRADECindices[5]].sum(axis=2).sum(axis=1)
        spectrum /= (boxSize*boxSize)
    elif region == 'yes':
        spectrum = data[:, 
                        minDECindex:maxDECindex, 
                        maxRAindex:minRAindex].sum(axis=2).sum(axis=1)
        spectrum /= ((maxDECindex-minDECindex)*(minRAindex-maxRAindex))
    else:
        objectRAindex, objectDECindex = degreeCoordinatesToCubeIndices(objectRA, objectDEC)
        spectrum = data[:, objectDECindex, objectRAindex]

    return spectrum

def getOnOffSpectrum(objectRA, objectDEC, boxSize=4):
    """
    this function extracts the full (~~~ -700 to 700 km/s) on-off spectrum centered around the target
    and returns the on-off spectrum and all of the off spectra, in that order
    # note the off spectrum is created by averaging the spectra of the eight boxes that surrounds the center box
    # note the boxes have length = boxSize and boxSize is expected in pixels (and again for us, 1 pixel = 1 arcMin)
    # note here the default boxSize is 1 GALFA beam 
        ^which I suppose isn't really all that useful but I'm assuming there'll always be an input I guess
    """
    raDecIndices = setPlotRange(objectRA, objectDEC, boxSize, boxSize)
    onSpectrum = getSpectrum(objectRA, objectDEC, boxSize=boxSize, box = 'yes')
    offSpectrum = np.zeros(2048)
    surroundingSpectrum = np.zeros((9,2048))
    goodSpectrum = 0

    for x in range(0,9):
        surroundingSpectrum[x] = getSpectrum(objectRA, objectDEC,
                                             maxRAindex=raDecIndices[4]+(boxSize*(x%3-1)), minRAindex=raDecIndices[5]+(boxSize*(x%3-1)),
                                             minDECindex=raDecIndices[6]+(boxSize*(1-x/3)), maxDECindex=raDecIndices[7]+(boxSize*(1-x/3)),
                                             region = 'yes')
        # this prevents the inclusion of NaNs
        if (x!=4):
            if np.isnan(surroundingSpectrum[x].sum()):
                print "Warning: getOnOffSpectrum: NaN in off spectrum " + str(x)
            else:
                offSpectrum += surroundingSpectrum[x]
                goodSpectrum += 1

    print goodSpectrum
    offSpectrum /= goodSpectrum
    onOffSpectrum = onSpectrum - offSpectrum
    
    return onOffSpectrum, surroundingSpectrum


### INTEGRATED INTENSITY ####################################
def getIntegratedIntensity(minVindex, maxVindex, intAxis = 0, vbinwidth=.18, noNegValues='yes'):
    '''
    this function calculates the integrated intensity for a given range of velocities
    and returns the full (for the entire cube) integrated intensity map
    # note this takes in velocity indicies not velocity
    # note intAxis: the axis on which the vel is on, intAxis = 0 by default
    # note (minV, maxV): the integration range     
    # note this is largely taken from Yong
    # note idk what vbinwidth is
    '''    
    if maxVindex-minVindex<1:
        print 'Velocity Channel mistake!'
        return 0
    else:
        intIntensity = np.sum(data[minVindex:maxVindex, :, :], axis=intAxis)*vbinwidth

    intIntensity *= 1.823*10**18    # in cm-2, assuming HI is optically thin. Draine, ISM book 10, page 74 # T km/s

    # this sets all the non-positive values to the smallist positive value
    # ^just setting them to NaNs seem to mess with saving the figure to PDF
    if noNegValues == 'yes':
        if np.isnan(intIntensity.sum()) or np.nanmin(intIntensity) <10**18:
            intIntensity[np.where(intIntensity<10**18)] = np.nan 
            minPosInt = np.nanmin(intIntensity)
            intIntensity[np.isnan(intIntensity)] = 10**18

    return intIntensity


### PLOTTING SPECTRUM ####################################
def plotSingleSpectrum(objectName, spectrum, minVzoom=0, maxVzoom=0, setVzoomRange = 'no', spectrumType = 'default', closefig = 'yes'):
    """
    this funciton plots a figure of the target's spectrum and saves in the the target's foldder
    the figure will contain two graphs:
    TOP: same as TOP but with a different intensity range that's centered around 0
    MIDDLE: the full spectrum with the target's systemic velocity marked
    BOTTOM: the zoomed in spectrum centered on the target's systemic velocity 
    # note this takes in the full specturm
    # note the default range of the zoomed in spectrum is systemic veloicty +- 75 km/s, unless otherwise stated
    # note default prints regular spectrum, input 'on-off' to change captions
    """
    if setVzoomRange == 'no':
        minVzoom = objectGamma - 75
        maxVzoom = objectGamma + 75
    elif setVzoomRange == 'yes':
        if minVzoom>maxVzoom:
            print "Illegal veloicty indices in plotSingleSpectrum call!"
            return 0
    else:
        print "plotSingleSpectrum call setVzoomRange field takes in 'yes' or 'no'"
        return 0

    minVindex, maxVindex = setVelocityRange(minVzoom, maxVzoom)
        
    if spectrumType == 'default':
        addtionalCaption = ''
        fileName = objectName + "_spectrum.pdf"
        yZoom = .3
    elif spectrumType == 'on-off':
        addtionalCaption = ' On-off'
        fileName = objectName + "_onOff_spectrum.pdf"
        yZoom = .25
        
    fig = plt.figure(figsize=(11.6,7.25))
    gs = gridspec.GridSpec(3,4, wspace=0, hspace=.51)
    ax2 = plt.subplot(gs[0,0:4]) #^yzoom
    ax1 = plt.subplot(gs[1,0:4]) #all velcoity
    ax3 = plt.subplot(gs[2,0:4]) #select velocity around gamma
    
    # yLowerBound & yUpperBound are here adjusted to the to the max of the spectrum
    yLowerBoundFactor = -0.09
    yUpperBoundFactor = 0.18
    yLowerBound = (np.max(spectrum)-np.min(spectrum)) * yLowerBoundFactor
    yUpperBound = (np.max(spectrum)-np.min(spectrum)) * yUpperBoundFactor


    ax1.plot(vlsr, spectrum)    
    ax1.set_xlim(vlsr[0], vlsr[-1])
    ax1.set_ylim(np.min(spectrum)+yLowerBound, np.max(spectrum)+yUpperBound)
    # plots the box where we zoom in
    ax1.plot([minVzoom, minVzoom, maxVzoom, maxVzoom, minVzoom],
             [np.min(spectrum)+yLowerBound*.6, np.max(spectrum)+yUpperBound*.1, np.max(spectrum)+yUpperBound*.1, 
              np.min(spectrum)+yLowerBound*.6, np.min(spectrum)+yLowerBound*.6], 
             color = 'black', linewidth = 0.81)
    ax1.text(minVzoom, np.min(spectrum)+yLowerBound*.5 , "Zoomed in below", fontsize = 3.9)
    # plots the expected velocity as a vertical line with notation and horizontal line indicating 0
    ax1.plot([vlsr[0],vlsr[-1]], [0,0], color = 'black', linewidth = 0.21)
    ax1.plot([objectGamma,objectGamma], [np.min(spectrum)+yLowerBound,np.max(spectrum)+yUpperBound],
             color = 'black', linewidth = 0.18)
    ax1.text(objectGamma, np.max(spectrum)+yUpperBound*.3,
             "Systemic Velocity = " + str(objectGamma) + " +-" + str(objectGammaError) + " km/s", fontsize = 6, ha = 'center')
    # titles and axis labels
    ax1.set_title("Full" + addtionalCaption + " Spectrum from GALFA_HI", fontsize = 10)
    ax1.set_xlabel("Velocity (km/s)", fontsize = 9)
    ax1.set_ylabel("Tb (K)", fontsize = 9)
    

    # same as ax1 but with different y scales
    ax2.plot(vlsr, spectrum)
    ax2.set_xlim(vlsr[0], vlsr[-1])
    ax2.set_ylim(-yZoom, yZoom)
    # plots the box where we zoom in
    ax2.plot([minVzoom, minVzoom], [-yZoom, yZoom], color = 'black', linewidth = 1)
    ax2.plot([maxVzoom, maxVzoom], [-yZoom, yZoom], color = 'black', linewidth = 1)
    ax2.text(minVzoom, -yZoom*.9, "Zoomed in below", fontsize = 3.9)

    # vertical line indicaiton gamma and horizontal line indicating 0
    ax2.plot([objectGamma,objectGamma], [-yZoom, yZoom], color = 'black', linewidth = 0.18)
    ax2.plot([vlsr[0],vlsr[-1]], [0,0], color = 'black', linewidth = 0.21)
    # titles axis labels
    ax2.set_title("Full" + addtionalCaption + " Spectrum Zoomed-in Around 0", fontsize = 10)
    ax2.set_xlabel("Velocity (km/s)", fontsize = 9)
    ax2.set_ylabel("Tb (K)", fontsize = 9)

    # yLowerBound & yUpperBound are here adjusted to the to the max of the spectrum
    yLowerBound = (np.max(spectrum[minVindex:maxVindex])-np.min(spectrum[minVindex:maxVindex])) * yLowerBoundFactor
    yUpperBound = (np.max(spectrum[minVindex:maxVindex])-np.min(spectrum[minVindex:maxVindex])) * yUpperBoundFactor
    

    # zoomed in
    ax3.plot(vlsr[minVindex:maxVindex], spectrum[minVindex:maxVindex])
    ax3.set_xlim(vlsr[minVindex], right=vlsr[maxVindex])
    ax3.set_ylim(np.min(spectrum[minVindex:maxVindex])+yLowerBound, np.max(spectrum[minVindex:maxVindex])+yUpperBound)
    # plots the expected velocity as a vertical line with notation and horizontal line indicating 0
    ax3.plot([vlsr[minVindex],vlsr[maxVindex]], [0,0], color = 'black', linewidth = 0.21)
    ax3.plot([objectGamma, objectGamma], 
             [np.min(spectrum[minVindex:maxVindex])+yLowerBound, np.max(spectrum[minVindex:maxVindex])+yUpperBound],
             color = 'black', linewidth = 0.18)
    ax3.text(objectGamma, np.max(spectrum[minVindex:maxVindex])+yUpperBound*.3, 
             "Systemic Velocity = " + str(objectGamma) + " +-" + str(objectGammaError) + " km/s", fontsize= 6, ha = 'center')
    # titles and axis labels 
    ax3.set_title(addtionalCaption + " Spectrum Zoomed-in Around the Systemic Velocity", fontsize = 10)
    ax3.set_xlabel("Velocity (km/s)", fontsize = 9)
    ax3.set_ylabel("Tb (K)", fontsize = 9)


    ### adding the zoom lines
    coordsA = 'data'
    coordsB = 'data'
    xyA = (minVzoom, np.min(spectrum)+yLowerBound*.6)
    xyB = (vlsr[minVindex], np.max(spectrum[minVindex:maxVindex])+yUpperBound)
    xyC = (maxVzoom, np.min(spectrum)+yLowerBound*.6)
    xyD = (vlsr[maxVindex], np.max(spectrum[minVindex:maxVindex])+yUpperBound)

    con1 = ConnectionPatch(xyA, xyB, coordsA = coordsA, coordsB = coordsB, axesA = ax1, axesB = ax3)
    con2 = ConnectionPatch(xyC, xyD, coordsA = coordsA, coordsB = coordsB, axesA = ax1, axesB = ax3)

    ax1.add_artist(con1)
    ax1.add_artist(con2)

    plt.suptitle(objectName + addtionalCaption + " Spectrum", fontsize = 15)
    fig.savefig("Results/" + objectName + "/" + cycle + "/" + fileName)
    if closefig == 'yes':
        plt.close()

def plotNearbySpectrum(objectName, centerRA, centerDEC, boxSize, minV=0, maxV=0, setVrange = 'no', closefig = 'yes'):
    """
    this function plots the target spectrum and the spectrum of the nearby (the 8 immidiate) squares
    the center square is where the target is
    # note the default sets the plotting velocity range to be gamma +- 100 km/s, unless setVRange = 'yes'
    """    
    if setVrange == 'no':
        minVindex, maxVindex = setVelocityRange(objectGamma-100, objectGamma+100)
    elif setVrange == 'yes':
        minVindex, maxVindex = setVelocityRange(minV, maxV)
    else:
        print "plotNearbySpectrum call setVrange field takes in 'yes' or 'no'"
    
    onOffSpectrum, allSpectrum = getOnOffSpectrum(centerRA, centerDEC, boxSize)
    maxInt = 1.15*(np.max(allSpectrum)-np.min(allSpectrum))

    fig = plt.figure(figsize=(6,4))
    grid = gridspec.GridSpec(3,3,wspace=0.06,hspace=0.06)
    
    for x in range(0,9):
        ax = plt.Subplot(fig,grid[x])
        ax.plot(vlsr[minVindex:maxVindex], allSpectrum[x,minVindex:maxVindex])
        # plots the systemic velocity line in the center graph
        if (x == 4):
            ax.plot([objectGamma,objectGamma],[-5,np.max(allSpectrum[x, minVindex:maxVindex])+5],
                    color = 'black', linewidth = 0.25)
            ax.text(objectGamma, np.max(allSpectrum[x, minVindex:maxVindex])*1.05, 
                    "Systemic Velocity = " + str(objectGamma) + " +-" + str(objectGammaError) + " km/s", fontsize= 5.1, ha = 'center')
        # sets and changes all the axes, ticks and labels
        ax.set_xlim(vlsr[minVindex], vlsr[maxVindex])
        ax.set_ylim(bottom=-5, top=maxInt)
        ax.tick_params(axis='y', which='major', labelsize = 7.2)
        ax.tick_params(axis='x', which='major', labelsize = 6.3)
        if (x>=0 and x<=5):
            ax.set_xticklabels([])
        if (x!=0 and x!=3 and x!=6):
            ax.set_yticklabels([])
        if (x == 3):
            ax.set_ylabel("Tb (K)", fontsize = 10)
        if (x ==7):
            ax.set_xlabel("Velocity (km/s)", fontsize = 10)
        fig.add_subplot(ax)

    plt.suptitle(objectName + " & nearby squares", fontsize=12)

    fileName = objectName + "_nearby_spectrums.pdf"
    fig.savefig("Results/" + objectName + "/" + cycle + "/" + fileName)
    if closefig == 'yes':
        plt.close()


### PLOTTING VELOCITY CHANNEL MAP ####################################
def makeVelocityChannelMaps(objectName, centerRA, centerDEC, raPlotRange, decPlotRange, vSeparation=5, closefig = 'yes'):
    """
    this funciton plots the target's velocity channel maps 
    the figure will be centered on the target's systemic velocity, AKA the middle graph have vlsr = gamma
    # note raPlotRange and decPlotRange are expected in pixels (again, 1 pixel = 1 arcMin for us)
    # note the default channel separation is 5 km/s
    """
    plotRADECindices = setPlotRange(centerRA, centerDEC, raPlotRange, decPlotRange)
    
    fig = plt.figure(figsize=(9.8,7))
    grid = gridspec.GridSpec(5,7,wspace=0.0,hspace=0.0)

    graphV = objectGamma + vSeparation * (np.arange(35)-17)
    vIndex = 1024 + (graphV)//delv

    # this picks out the absolute max in the range of cubes, used later for scaling the integrated intensity maps
    maxIntensity = np.max(data[vIndex[0]:vIndex[-1],
                               plotRADECindices[6]:plotRADECindices[7], 
                               plotRADECindices[4]:plotRADECindices[5]])

    for x in range(0,35):
        ax = plt.Subplot(fig,grid[x])
        plotData = data[vIndex[x], 
                        plotRADECindices[6]:plotRADECindices[7], 
                        plotRADECindices[4]:plotRADECindices[5]]
        ax.scatter(objectRA, objectDEC, c='red', marker = '+')
        ax.imshow(plotData, extent=[plotRADECindices[0], plotRADECindices[1], plotRADECindices[2], plotRADECindices[3]], 
            origin='lower', cmap = 'Greys', vmax = maxIntensity)
        if (x == 17):
            ax.text(plotRADECindices[0]-delta,plotRADECindices[2]+delta, str(graphV[x]) + " km/s", fontsize = 8.1, weight ='bold', color = 'blue')
        else:
            ax.text(plotRADECindices[0]-delta,plotRADECindices[2]+delta, str(graphV[x]) + " km/s", fontsize = 8.1, color = 'blue')
        ax.tick_params(axis = 'x', which = 'major', labelsize = 5.1)
        ax.tick_params(axis = 'y', which = 'major', labelsize = 7.2)
        if (x>=0 and x<=27):
            ax.set_xticklabels([])
        if (x!=0 and x!=7 and x!=14 and x!=21 and x!=28):
            ax.set_yticklabels([])   
        if (x == 31):
            ax.set_xlabel("RA(J2000)", fontsize = 12)
        if (x == 14):
            ax.set_ylabel("DEC(J2000)", fontsize = 12)
        fig.add_subplot(ax)


    plt.suptitle(objectName + " channel map from GALFA_HI data", fontsize=15)

    fileName = objectName + "_channel_map.pdf"
    fig.savefig("Results/" + objectName + "/" + cycle + "/" + fileName)
    if closefig == 'yes':
        plt.close()


### PLOTTING INTEGRATED INTENSITY MAP ####################################
def makeIntegratedIntensityMap(objectName, centerRA, centerDEC, raPlotRange, decPlotRange, minV=0, maxV=0, setVrange = 'no', closefig = 'yes'):
    """
    this function plots the 2D image (integrated intensity) of a GALFA cube
    # note raPlotRange and decPlotRange are expected in pixels (again, 1 pixel = 1 arcMin for us)
    # note minV and maxV are for setting the velocity range over which integrated intensity is calculated
    # note default ^range would be gamma +- vRange km/s
    """
    if setVrange == 'no':
        vRange = 15
        minV = objectGamma - vRange
        maxV = objectGamma + vRange
    elif setVrange == 'yes':
        if minV>maxV:
            print "invalid velocity range"
    else:
        print "makeIntegratedIntensityMap call setVrange field takes in 'yes' or 'no'"
        return 0
        
    minVindex, maxVindex = setVelocityRange(minV, maxV)
    plotRADECindices = setPlotRange(centerRA, centerDEC, raPlotRange, decPlotRange)
    intensityMap = getIntegratedIntensity(minVindex, maxVindex)
    intensityMap = intensityMap[plotRADECindices[6]:plotRADECindices[7], plotRADECindices[4]:plotRADECindices[5]]
    if np.min(intensityMap) <= 0:
        indGreater = np.where(intensityMap>0)
        indZero = np.where(intensityMap==0)
        indLess = np.where(intensityMap<0)
        intensityMap[indGreater] = np.log10(intensityMap[indGreater])
        intensityMap[indLess] = np.log10(-1*intensityMap[indLess])
    else:
        intensityMap = np.log10(intensityMap)

    contourLevelNum = 7
    contourLevel = np.arange(contourLevelNum) * (np.max(intensityMap)-np.min(intensityMap))/(contourLevelNum-1)
    contourLevel += np.min(intensityMap)

    #fig = plt.figure(figsize=(7.2,4.8))
    fig = plt.figure(figsize=(7,5))

    
    plt.imshow(intensityMap, 
               extent=[plotRADECindices[0], plotRADECindices[1], plotRADECindices[2], plotRADECindices[3]],
               origin='lower', cmap = 'Greys')
    plt.colorbar(shrink = 0.8, ticks = contourLevel, format = '%.2f', label = 'log10(#/cm^2)')
    plt.suptitle(objectName + " Integrated Intensity Map" + " (" + str(minV) + " to " + str(maxV) + " km/s)", fontsize = 12.9)
    # axis and labels 
    plt.xlim(plotRADECindices[0], plotRADECindices[1])
    plt.ylim(plotRADECindices[2], plotRADECindices[3])
    plt.ylabel("DEC (J2000)")
    plt.xlabel("RA (J2000)")
    # plots star marker("+")
    plt.scatter(objectRA, objectDEC, c='red', marker = '+')
    # plots contour
    cp = plt.contour(intensityMap,
                     extent=[plotRADECindices[0], plotRADECindices[1], plotRADECindices[2], plotRADECindices[3]],
                     colors = 'black', linewidths = 0.3, levels = contourLevel)
    plt.clabel(cp,inline=1, fontsize = 3.6, fmt = '%.1f')
    # plots scale (1 GALFA beam = 4' = 4 pixels)
    plt.plot([plotRADECindices[0] - (raPlotRange/30.0)*delta, plotRADECindices[0] - (raPlotRange/30.0 + 4)*delta],
             [plotRADECindices[2] + (decPlotRange/30.0)*delta, plotRADECindices[2] + (decPlotRange/30.0)*delta],
             color = 'b', linewidth = 1.5)
    plt.text(plotRADECindices[0] - (raPlotRange/30.0)*delta, 
             plotRADECindices[2] + (decPlotRange/30.0)*1.5*delta,
             "Beam=4'", fontsize = 6.6, color = 'blue')

    fileName = objectName + "_image(" + str(raPlotRange/60.) + "*" + str(decPlotRange/60.) + ").pdf"
    fig.savefig("Results/" + objectName + "/" + cycle + "/" + fileName)

    if closefig == 'yes':
        plt.close()


### PLOTTING _ALL MAP ####################################
def makeCombinedFigure(objectName, smallBoxSize, bigBoxSize, minIntV=0, maxIntV=0, setIntensityVrange = 'no', minV = 0, maxV = 0, setVrange = 'no', closefig = 'yes'):
    """
    this functions makes a combined output figure
    includes: integrated intensity maps, target spectrum and on-off spectrum
    this plots two sets of graphs: smallbox and bigbox, with box refering to the region which we extract the spectrum
    this is done because we don't know how big the target nova shells are
    each set of graphs contains: 
    an integrated intensity map, 
    an on spectrum and the zoom-in of that on spectrum, 
    an on-off spectrum and the zoom-in of that on-off spectrum
    # note smallBoxSize and bigBoxSize are expected in pixels (again, 1 pixel = 1 arcMin for us)
    # note minIntV and maxI tV sets the velocity range over which we compute the integrated intensity 
    # ^note the default range is objectGamma +- vRange km/s
    # note minV and maxV sets the zoomed in velocity range for target's spectrum (and on-off spectrum)
    # ^note the default range is objectGamma +- 50 km/s
    """
    if setIntensityVrange == 'no':
        vRange = 15
        minIntV = objectGamma - vRange
        maxIntV = objectGamma + vRange
    elif setIntensityVrange == 'yes':
        if minIntV>maxIntV:
            print "invalid int velocity range"
    else:
        print "makeCombinedFigure call setIntensityVrange field takes in 'yes' or 'no'"
    if setVrange == 'no':
        minV = objectGamma - 50
        maxV = objectGamma + 50
    elif setVrange == 'yes':
        if minV>maxV:
            print "invalid velocity range"
    else:
        print "makeCombinedFigure call setVrange field takes in 'yes' or 'no'"
        
    minVindex, maxVindex = setVelocityRange(minV, maxV) #this is for the spectrum
    minIntVindex, maxIntVindex = setVelocityRange(minIntV, maxIntV) #this is for the integrated intensity map
    fullIntensityMap = getIntegratedIntensity(minIntVindex, maxIntVindex)
    
    fig = plt.figure(figsize=(12.6,12.6))
    gs = gridspec.GridSpec(8,8, wspace=.6, hspace=.45)
    
    ######### loops through two sets of figures #########
    for x in range(0,2):
        
        ax1 = plt.subplot(gs[0+x*4:4+x*4,0:4]) #map
        ax12 = plt.subplot(gs[0+x*4,4:8]) #spectrum
        ax13 = plt.subplot(gs[1+x*4,4:8]) #on-off spectrum
        ax121 = plt.subplot(gs[2+x*4,4:8]) #spectrum vZoom
        ax131 = plt.subplot(gs[3+x*4,4:8]) #on-off spectrum vZoom
        
        if x == 0:
            graphBoxSize = smallBoxSize*5
            boxSize = smallBoxSize
        elif x ==1:
            graphBoxSize = bigBoxSize*3
            boxSize = bigBoxSize

        ### plotting the intensity map
        plotRADECindices1 = setPlotRange(objectRA, objectDEC, graphBoxSize, graphBoxSize)
        intensityMap = fullIntensityMap[plotRADECindices1[6]:plotRADECindices1[7], plotRADECindices1[4]:plotRADECindices1[5]]
        intensityMap = np.log10(intensityMap)
        ax1.imshow(intensityMap, 
                   extent=[plotRADECindices1[0], plotRADECindices1[1], plotRADECindices1[2], plotRADECindices1[3]],
                   origin='lower', cmap = 'Greys')
        # plotting the target marker
        ax1.scatter(objectRA, objectDEC, c='red', marker = '+')
        # plotting the scale, 1 GALFA beam = 4'
        ax1.plot([plotRADECindices1[0]-(graphBoxSize/30.0)*delta, plotRADECindices1[0]-(graphBoxSize/30.0 + 4)*delta],
                 [plotRADECindices1[2]+(graphBoxSize/30.0)*delta, plotRADECindices1[2]+(graphBoxSize/30.0)*delta],
                 color = 'b', linewidth = 1.5)
        ax1.text(plotRADECindices1[0] - (graphBoxSize/30.0)*delta, 
                 plotRADECindices1[2] + (graphBoxSize/30.0)*1.5*delta,
                 "Beam=4'", fontsize = 6.6, color = 'blue')
        # axis and labels 
        ax1.set_xlim(plotRADECindices1[0], plotRADECindices1[1])
        ax1.set_ylim(plotRADECindices1[2], plotRADECindices1[3])
        ax1.set_xlabel("RA (J2000)")
        ax1.set_ylabel("DEC (J2000)")
        if x==0:
            ax1.set_title("Integrated Intensity", fontsize = 12.9)
        # resetting things for spectrum
        plotRADECindices2 = setPlotRange(objectRA, objectDEC, boxSize, boxSize)
        # box over which the target spectrum is extracted
        ax1.plot([plotRADECindices2[0],plotRADECindices2[0],plotRADECindices2[1],plotRADECindices2[1],plotRADECindices2[0]],
                 [plotRADECindices2[2],plotRADECindices2[3],plotRADECindices2[3],plotRADECindices2[2],plotRADECindices2[2]],
                 color = 'black', linewidth = 0.9)
        ax1.text(plotRADECindices2[1], plotRADECindices2[2], str(boxSize) + "' ", va='top', ha = 'right')
        ax1.text(plotRADECindices1[1], plotRADECindices1[2], str(graphBoxSize) + "'.", va='bottom', ha = 'right')


        ### plotting the full on spectrum 
        spectrum = getSpectrum(objectRA, objectDEC, boxSize, box='yes')
        # set the upper and lower bounds
        yLowerBoundFactor = -0.09
        yUpperBoundFactor = 0.18
        yLowerBound = (np.max(spectrum)-np.min(spectrum)) * yLowerBoundFactor
        yUpperBound = (np.max(spectrum)-np.min(spectrum)) * yUpperBoundFactor
        ax12.plot(vlsr, spectrum)    
        ax12.set_xlim(vlsr[0], vlsr[-1])
        ax12.set_ylim(np.min(spectrum)+yLowerBound, np.max(spectrum)+yUpperBound)
        # plots and labels the box where we zoom in
        ax12.plot([minV, minV, maxV, maxV, minV],
                  [np.min(spectrum)+yLowerBound*.6, np.max(spectrum)+yUpperBound*.1, np.max(spectrum)+yUpperBound*.1, 
                  np.min(spectrum)+yLowerBound*.6, np.min(spectrum)+yLowerBound*.6],
                  color = 'black', linewidth = 0.81)
        ax12.text(minV, np.min(spectrum)+yLowerBound*.5 , "Zoomed in below", fontsize = 3)
        # plots the expected velocity as a vertical line with notation and horizontal line indicating 0
        ax12.plot([vlsr[0],vlsr[-1]], [0,0], color = 'black', linewidth = 0.3)
        ax12.plot([objectGamma,objectGamma], [np.min(spectrum)+yLowerBound,np.max(spectrum)+yUpperBound],
                 color = 'black', linewidth = 0.3)
        # titles and axis labels
        ax12.tick_params(axis = 'both', which = 'major', labelsize = 6)
        ax12.set_ylabel("Tb (K)", fontsize = 7)
        ax12.set_title("Spectrum", fontsize = 9)
                

        ### plotting the zoomed in spectrum
        ax121.plot(vlsr[minVindex:maxVindex], spectrum[minVindex:maxVindex])
        # reset the upper and lower bounds
        yLowerBoundZoom = (np.max(spectrum[minVindex:maxVindex])-np.min(spectrum[minVindex:maxVindex])) * yLowerBoundFactor
        yUpperBoundZoom = (np.max(spectrum[minVindex:maxVindex])-np.min(spectrum[minVindex:maxVindex])) * yUpperBoundFactor
        ax121.set_xlim(vlsr[minVindex], right=vlsr[maxVindex])
        ax121.set_ylim(np.min(spectrum[minVindex:maxVindex])+yLowerBoundZoom, np.max(spectrum[minVindex:maxVindex])+yUpperBoundZoom)
        # plots the expected velocity as a vertical line with notation and horizontal line indicating 0
        ax121.plot([vlsr[minVindex],vlsr[maxVindex]], [0,0], color = 'black', linewidth = 0.3)
        ax121.plot([objectGamma, objectGamma], 
                   [np.min(spectrum[minVindex:maxVindex])+yLowerBoundZoom, np.max(spectrum[minVindex:maxVindex])+yUpperBoundZoom],
                   color = 'black', linewidth = 0.3)
        ax121.text(objectGamma, np.max(spectrum[minVindex:maxVindex])+yUpperBoundZoom*.3, 
                 "Systemic Velocity = " + str(objectGamma) + " +-" + str(objectGammaError) + " km/s", fontsize= 5.4, ha = 'center')
        # titles and axis labels 
        ax121.tick_params(axis = 'both', which = 'major', labelsize = 6)
        ax121.set_ylabel("Tb (K)", fontsize = 7)
        ax121.set_title("Zoomed-in Spectrum", fontsize = 8.1)


        ### plotting the full on-off spectrum
        onOffSpectrum, offSpectrum = getOnOffSpectrum(objectRA, objectDEC, boxSize)
        ax13.plot(vlsr, onOffSpectrum)   
        # reset the upper and lower bounds
        yLowerBoundOnOff = (np.max(onOffSpectrum)-np.min(onOffSpectrum)) * yLowerBoundFactor
        yUpperBoundOnOff = (np.max(onOffSpectrum)-np.min(onOffSpectrum)) * yUpperBoundFactor 
        ax13.set_xlim(vlsr[0], vlsr[-1])
        ax13.set_ylim(np.min(onOffSpectrum)+yLowerBoundOnOff, np.max(onOffSpectrum)+yUpperBoundOnOff)
        # plots the labels the box where we zoom in
        ax13.plot([minV, minV, maxV, maxV, minV],
                  [np.min(onOffSpectrum)+yLowerBoundOnOff*.6, np.max(onOffSpectrum)+yUpperBoundOnOff*.1, np.max(onOffSpectrum)+yUpperBoundOnOff*.1, 
                  np.min(onOffSpectrum)+yLowerBoundOnOff*.6, np.min(onOffSpectrum)+yLowerBoundOnOff*.6],
                  color = 'black', linewidth = 0.81)
        ax13.text(minV, np.min(onOffSpectrum)+yLowerBoundOnOff*.5 , "Zoomed in below", fontsize = 3)
        # plots the expected velocity as a vertical line with notation and horizontal line indicating 0
        ax13.plot([vlsr[0],vlsr[-1]], [0,0], color = 'black', linewidth = 0.3)
        ax13.plot([objectGamma,objectGamma], [np.min(onOffSpectrum)+yLowerBoundOnOff,np.max(onOffSpectrum)+yUpperBoundOnOff],
                 color = 'black', linewidth = 0.3)
        # titles and axis labels
        ax13.tick_params(axis = 'both', which = 'major', labelsize = 6)
        ax13.set_ylabel("Tb (K)", fontsize = 7)
        ax13.set_title("On-off Spectrum", fontsize = 9)


        ### plotting the zoomed in on-off spectrum
        ax131.plot(vlsr[minVindex:maxVindex], onOffSpectrum[minVindex:maxVindex])
        yLowerBoundOnOffZoom = (np.max(onOffSpectrum[minVindex:maxVindex])-np.min(onOffSpectrum[minVindex:maxVindex])) * yLowerBoundFactor
        yUpperBoundOnOffZoom = (np.max(onOffSpectrum[minVindex:maxVindex])-np.min(onOffSpectrum[minVindex:maxVindex])) * yUpperBoundFactor 
        ax131.set_xlim(vlsr[minVindex], right=vlsr[maxVindex])
        ax131.set_ylim(np.min(onOffSpectrum[minVindex:maxVindex])+yLowerBoundOnOffZoom, np.max(onOffSpectrum[minVindex:maxVindex])+yUpperBoundOnOffZoom)
        # plots the expected velocity as a vertical line with notation and horizontal line indicating 0
        ax131.plot([vlsr[minVindex],vlsr[maxVindex]], [0,0], color = 'black', linewidth = 0.3)
        ax131.plot([objectGamma, objectGamma], 
                 [np.min(onOffSpectrum[minVindex:maxVindex])+yLowerBoundOnOffZoom, np.max(onOffSpectrum[minVindex:maxVindex])+yUpperBoundOnOffZoom],
                 color = 'black', linewidth = 0.3)
        ax131.text(objectGamma, np.max(onOffSpectrum[minVindex:maxVindex])+yUpperBoundOnOffZoom*.3, 
                 "Systemic Velocity = " + str(objectGamma) + " +-" + str(objectGammaError) + " km/s", fontsize= 5.4, ha = 'center')
        # titles and axis labels 
        ax131.tick_params(axis = 'both', which = 'major', labelsize = 6)
        ax131.set_ylabel("Tb (K)", fontsize = 7)
        ax131.set_title("Zoomed-in On-off Spectrum", fontsize = 8.1)


        ### adding the zoom lines
        coordsA = 'data'
        coordsB = 'data'
        xyA = (plotRADECindices2[1], plotRADECindices2[3])
        xyB = (vlsr[0], np.max(spectrum)+yUpperBound)
        xyC = (plotRADECindices2[1], plotRADECindices2[2])
        xyD = (vlsr[0], np.min(onOffSpectrum)+yLowerBoundOnOff)

        con1 = ConnectionPatch(xyA, xyB, coordsA = coordsA, coordsB = coordsB, axesA = ax1, axesB = ax12)
        con2 = ConnectionPatch(xyC, xyD, coordsA = coordsA, coordsB = coordsB, axesA = ax1, axesB = ax13)

        ax1.add_artist(con1)
        ax1.add_artist(con2)


    plt.suptitle(objectName, fontsize=15)
    fileName = objectName + "_all_maps.pdf"
    fig.savefig("Results/" + objectName + "/" + cycle + "/" + fileName)
    if closefig == 'yes':
        plt.close()














