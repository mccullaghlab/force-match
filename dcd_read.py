#!/Users/martinmccullagh/anaconda/bin/python


import scipy
import matplotlib.pyplot as plt
import numpy
import MDAnalysis

configFile = 'forcematchFiles.config'
psf = None
forceDcd = None
coordDcd = None
force = None
coord = None
debug = False

# default coordinate max/min and binsize
rMin = 3
rMax = 25
rBinSize = 0.1
nBins = int((rMax - rMin) / rBinSize)

# average force array
fMagAvg = numpy.zeros(nBins, dtype=float)
# count array
fMagCount = numpy.zeros(nBins, dtype=int)

print("T:",len(fMagAvg),nBins)

# Set debug mode from config file
def setDebug(cF):
    global debug
    txt = open(cF, 'r')
    line = txt.next()
    while line != "END CONFIG\n":
        if line == "DEBUG MODE: ON\n":
            debug = True
        line = txt.next()

# Get name of PSF file from config file
def getPsf(cF):
    global psf
    print('\n\t***DCD Analysis***')
    if debug:
        print('\t\tDebug Mode ON')
    else:
        print('\t\tDebug Mode OFF')
    txt = open(cF, 'r')
    while psf is None:
        line = txt.next()
        if line == 'PSF FILE:\n':
            psf = txt.next()[:-1]
            if debug:
                print('PSF File: %s' % psf)
        elif line == 'END CONFIG FILE\n':
            print('No PSF file found in config.')
            break

# Get name of Force DCD files from config file
def getForceDCDs(cF):
    global forceDcd
    forceDcd = []
    txt = open(cF, 'r')
    while len(forceDcd) == 0:
        line = txt.next()
        if line == 'FORCE DCD FILES:\n':
            line = txt.next()
            while line != '\n':
                if line == 'END CONFIG\n':
                    print('NO FORCE DCD FILES FOUND')
                forceDcd.append(line[:-1])
                line = txt.next()
            if debug:
                print('Force DCD files found: %s' % forceDcd)

# Get name of Coordinate DCD files from config file
def getCoordDCDs(cF):
    global coordDcd
    coordDcd = []
    txt = open(cF, 'r')
    while len(coordDcd) == 0:
        line = txt.next()
        if line == 'COORD DCD FILES:\n':
            line = txt.next()
            while line != '\n':
                if line == 'END CONFIG\n':
                    print('NO FORCE DCD FILES FOUND')
                coordDcd.append(line[:-1])
                line = txt.next()
            if debug:
                print('Coordinate DCD files found: %s\n' % coordDcd)

# Set coordinate max/min and binsize
def getCoordBounds(cF):
    txt = open(cF,'r')
    line = txt.next()
    length1 = len("COORDINATE MIN: ")
    length2 = len("COORDINATE MAX: ")
    length3 = len("BIN SIZE: ")

    global rMin
    global rMax
    global rBinSize

    # scan config file for coord and bin values
    while line != "END CONFIG\n":
        line = txt.next()
        if line[:length1] == "COORDINATE MIN: ":
            rem = -1 * (len(line) - length1)
            rMin = int(line[rem:-1])
        elif line[:length2] == "COORDINATE MAX: ":
            rem = -1 * (len(line) - length2)
            rMax = int(line[rem:-1])
        elif line[:length3] == "BIN SIZE: ":
            rem = -1 * (len(line) - length3)
            rBinSize = float(line[rem:-1])
    global nBins
    global fMagAvg
    global fMagCount
    nBins = int((rMax - rMin) / rBinSize)
    # average force array
    fMagAvg = numpy.zeros(nBins, dtype=float)
    # count array
    fMagCount = numpy.zeros(nBins, dtype=int)
    if debug:
        print "--Dimensions of System--"
        print "\tCoord Min: {}".format(rMin)
        print "\tCoord Max: {}".format(rMax)
        print "\tBin Size: {}".format(rBinSize)
        print "\tNumber of bins:", nBins

# Define subset of data without solvent
def parseWater():
    # select all atoms that are not water or hydrogen
    if debug:
        print("\nReading DCD data:\n\tParsing out WATER...\n")
    global ionsCoord
    global ionsForce
    ionsCoord = coord.selectAtoms("not name H1 and not name H2 and not name OH2")
    ionsForce = force.selectAtoms("not name H1 and not name H2 and not name OH2")

# Initialize MD Analysis
def initMDA():
    # start MDAnalysis with a universal
    # force universe
    global force
    global coord
    global debug

    force = MDAnalysis.Universe(psf, forceDcd)
    # coordinate universe
    coord = MDAnalysis.Universe(psf, coordDcd)

    # Truncate solvent out of the data
    parseWater()

    # Print log data
    printLogData(debug)

# Print log data
def printLogData(d):
    if d:
        global ionsCoord
        global ionsForce
        # print list of atoms being force matched
        print "DCD coord universe after course-graining:", len(ionsCoord), "atom(s)"
        for i in range(0, len(ionsCoord)):
            print "\t", ionsCoord[i]
        print "Force DCD universe after course-graining:", len(ionsForce), "atom(s)"
        for i in range(0, len(ionsForce)):
            print "\t", ionsForce[i]
        # some general log info

        print "\nNumber of time steps in coordinate trajectory:", len(coord.trajectory)
        print "Number of time steps in force trajectory:", len(force.trajectory)
        if len(coord.trajectory) == len(force.trajectory):
            print("Coord and Force time step counts MATCH\n")
        else:
            print("Coord and Force time step counts DO NOT MATCH\nCheck .config file for incorrect .dcd or .force.dcd file names.")

# Iterate over all time steps and perform MD calculations
def iterateAndCompute():

    # choose a time step to be referenced as an example when debugging
    exampleTs = int(len(coord.trajectory) / 2)
    # junk variable for log purposes
    j = 0

    # Debugging reference time step
    if debug:
        print("Referencing time step {} for logging purposes.\n".format(exampleTs))

    # initialize distance vector
    r = [0., 0., 0.]

    # get box dimensions (does not change between time steps)
    box = coord.trajectory.ts.dimensions[:3]



    # iterate over coordinate trajectory
    for ts in coord.trajectory:

        # iterate through all combinations of coord ion pairs
        for a in range(0, len(ionsCoord)):
            for b in range(a, len(ionsCoord)):
                # don't compute the force of an atom on itself
                if (a != b):
                    # compute 3D distance vector, r
                    for i in range(0, 3):
                        r[i] = ionsCoord[a].pos[i] = ionsCoord[b].pos[i]

                    # check which image is closest
                    if r[i] > box[i] / 2:
                        r[i] -= box[i]
                    elif r[i] < -box[i] / 2:
                        r[i] += box[i]
                    # compute the magnitude of the distance vector (aka the distance)
                    magR = numpy.linalg.norm(r)
                    # compute the normal vector in the distance direction
                    rHat = r / magR
                    # debug print step, (atom A, atom B, Distance from A->B)
                    # print("%i, %i, %f\n" % (a, b, magR))

                    # Now work in Force Universe
                    if (ts.frame == exampleTs) & debug & (j == 0):
                        print("Force computations for timestep {}:\n".format(exampleTs))
                        j = 1
                    ionAForce = ionsForce[a].pos
                    if (ts.frame == exampleTs) & debug:
                        print("Atom pair (%d, %d)" % (a, b))
                        print("\tAtom A force traj: {}".format(ionAForce))
                    ionBForce = ionsForce[b].pos
                    if (ts.frame == exampleTs) & debug:
                        print("\tAtom B force traj: {}".format(ionBForce))

                    # compute the projection of the force in the distance direction
                    AMagF = numpy.linalg.norm(numpy.dot(ionAForce, rHat))
                    BMagF = numpy.linalg.norm(numpy.dot(ionBForce, rHat))

                    if (ts.frame == exampleTs) & debug:
                        print("\t\tMagnitude of Force of atom A: {}".format(AMagF))
                        print("\t\tMagnitude of Force of atom B: {}".format(BMagF))

                    # First compute array index of magR
                    magRbin = int((magR - rMin) / rBinSize)
                    if (ts.frame == exampleTs) & debug:
                        print("\tMagnitude of R = {} -> Bin: {}".format(magR, magRbin))

                    # add to average force array
                    if 0 <= magRbin < nBins:
                        fMagAvg[magRbin] += AMagF
                        fMagAvg[magRbin] += BMagF
                        # add to counts
                        fMagCount[magRbin] += 2

    bins = numpy.arange(nBins)
    makeScatter(bins, "R", fMagCount, "F Magnitude", "F Mag vs Distance R")


# display 2D scatter plot of some data
def makeScatter(x, xLabel, y, yLabel, plotTitle):
    fig = plt.figure()
    plt.scatter(x, y, 30)
    fig.suptitle(plotTitle, fontsize=20)
    plt.xlabel(xLabel, fontsize=18)
    plt.ylabel(yLabel, fontsize=18)
    plt.show()

'''
    Left off after cleaning up the iterateAndCompute() function.
    Need to verify that the avg force array data is computed correctly, and
    attempt to visualize it.
    Then continue building functions for the remaining free-floating code.
    Also research non-bonding energy and hamiltonians.

'''

# Read config setting for debug mode
setDebug(configFile)
# Get name of PSF file from config file
getPsf(configFile)
# Get names of Force DCD files from config file
getForceDCDs(configFile)
# Get names of Coord DCD files from config file
getCoordDCDs(configFile)
# Define coordinate min/max and bin size
getCoordBounds(configFile)
# Initialize MD Analysis
initMDA()
# Iterate over time steps, and perform MD calculations
iterateAndCompute()









'''
# complete averaging and pring
for i in range(0, nBins):
    if fMagCount[i] == 0:
        print i * rBinSize, 0.0, fMagCount[i]
    else:
        print i * rBinSize, fMagAvg[i] / float(fMagCount[i]), fMagCount[i]'''
