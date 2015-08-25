"""
This Python script will go through all the targrets in the select text file and produce ......

This script is written by Larry with input from Mary, Jeno, Josh and Yong
08/04/2015
"""
import run_GALFA_functions_Larry as pf # as in "Plotting Functions"
import numpy as np
import datetime

############################## RUNNING ##############################
# this will read in all the targets and print all maps
# coordinate file info
coordinateFileName = 'currentTargets.txt'
coordinateDir = 'coordinateFiles/'
coordinateFileDir = coordinateDir+  coordinateFileName
cycle = 'cycle1'

objectInfoList = np.genfromtxt(coordinateFileDir, names=True, dtype=None)
numTargets = objectInfoList.size

print "\n######### RUNNING run_all_CVs_GALFA_Larry #########"
print "######### " + str(numTargets) + " targets found #########"

processed = 0

for x in range(0,numTargets):
    # getting object information
    thisObjectInfo = objectInfoList[x]
    objectName = thisObjectInfo[0]
    objectRA = thisObjectInfo[1]
    objectDEC = thisObjectInfo[2]
    objectGamma = thisObjectInfo[3]
    objectGammaError = thisObjectInfo[4]
    pf.setTargetInfo(name=objectName, raCoord=objectRA, decCoord=objectDEC, gamma=objectGamma, gammaE=objectGammaError, c=cycle)

    print "### Target: " + objectName

    # checking if folder exist
    pf.checkTargetFolder(objectName, cycle)
    # getting the cube information and get the data from the cube
    cubeFileInfo = pf.getTargetCubeinfo(objectRA, objectDEC)

    # prints out info into the text file
    targetFileName = objectName + '_targetFile.txt'
    f = open('Results/' + objectName + '/' + cycle + '/' + targetFileName, 'w')
    f.write("Target name: " + objectName + "\n")
    f.write("This is " + cycle  + "\n")
    f.write("Target coordinates: \n")
    f.write("    RA: " + str(objectRA) + " deg \n")
    f.write("    DEC: " + str(objectDEC) + " deg \n")
    f.write("GALFA cube coordinates: " + str(cubeFileInfo) + "\n\n")

    # skipping cubes that I know are 'broken'
    # 292, 34.35: found when V587_Lyr crashed the script 
    if int(cubeFileInfo[0]) == 292 and int(cubeFileInfo[1]*100) == 3435:
        print "known bad cube, skipping this target... \n"
        f.write("Target lies in known bad cube, " + cycle + " skipped to the next target at " + str(datetime.datetime.today()))
        continue

    if x == 0:
        cube_name = 'GALFAdata/data/GALFA_HI_RA+DEC_%06.2f+%05.2f_W.fits' % (cubeFileInfo[0], cubeFileInfo[1])
        cubeData = pf.GALFAHI_cubeinfo(cube_name)
        print "### reading from GALFA cube %06.2f, %05.2f" % (cubeFileInfo[0], cubeFileInfo[1])
    elif x>0:
        lastCubeFileInfo = pf.getTargetCubeinfo(objectInfoList[x-1][1], objectInfoList[x-1][2])
        if cubeFileInfo != lastCubeFileInfo:
            cube_name = 'GALFAdata/data/GALFA_HI_RA+DEC_%06.2f+%05.2f_W.fits' % (cubeFileInfo[0], cubeFileInfo[1])
            cubeData = pf.GALFAHI_cubeinfo(cube_name)
            print "### reading from GALFA cube %06.2f, %05.2f" % (cubeFileInfo[0], cubeFileInfo[1])
        else:
            print "### same cube as last target"

    # takes out the different slices of cubedata, remember data:[v,dec,ra]
    data = cubeData[0]
    vlsr = cubeData[2]
    dec = cubeData[3]
    ra = cubeData[4]
    delv = cubeData[5]

    # this checks how far from the cube edge the targets are
    raInd, decInd = pf.degreeCoordinatesToCubeIndices(objectRA, objectDEC)
    up = dec.size - decInd
    down = decInd
    left = raInd
    right = ra.size - raInd
    print left, right, up, down
    f.write("Dist from edge of cube: " + str([left, right, up, down]) + "\n\n")

    if up>=45 and down>=45 and left>=45 and right>=45:
        ############################## PLOTTING CALLS ##############################
        # spectrum
        spectrum = pf.getSpectrum(objectRA, objectDEC, boxSize=16, box='yes')
        pf.plotSingleSpectrum(objectName, spectrum)
        f.write("Target spectrum plotted at " + str(datetime.datetime.today()) + "\n")

        # on-off spectrum
        onOffSpectrum = pf.getOnOffSpectrum(objectRA, objectDEC, boxSize=16)
        pf.plotSingleSpectrum(objectName, onOffSpectrum[0], spectrumType='on-off')
        f.write("On-off spectrum plotted at " + str(datetime.datetime.today()) + "\n")

        # nearby spectrum
        pf.plotNearbySpectrum(objectName, objectRA, objectDEC, boxSize=8)
        f.write("Spectra of the 9 nearby squares plotted at " + str(datetime.datetime.today()) + "\n")

        # integrated intensity map
        pf.makeIntegratedIntensityMap(objectName, objectRA, objectDEC, 90, 90)
        f.write("Integrated intensity map made at " + str(datetime.datetime.today()) + "\n")

        # channel map
        pf.makeVelocityChannelMaps(objectName, objectRA, objectDEC, 90, 90, vSeparation=5)
        f.write("Velocity channel maps made at " + str(datetime.datetime.today()) + "\n")

        # _ALL map
        pf.makeCombinedFigure(objectName, 6, 30)
        f.write("_all map made at " + str(datetime.datetime.today()) + "\n")

        # _s1 map
        pf.makeS1Map(objectName, objectRA, objectDEC, 30, 90)
        f.write("_s1 map made at " + str(datetime.datetime.today()) + "\n")

        processed += 1
        f.write("\n" + cycle + " completed for " + objectName + " at " + str(datetime.datetime.today()) + "\n")
        print "######### " + objectName + " " + cycle + " complete #########" + " (" + str(x+1) + " out of " + str(numTargets) + " targets) \n"
    else:
        print "######### " + objectName + " " + cycle + " incomplete. Too close to edge of cube"
        f.write("\n" + cycle + " incomplete for " + objectName + " at " + str(datetime.datetime.today()) + "\n")
        f.write("Target too close to the edge of the cube")

    f.close()

print "Done, " + str(processed) + " out of " + str(numTargets) + " targets fully processed. \n"

