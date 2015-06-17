__author__ = 'Greg Poisson'

import matplotlib.pyplot as plt
import numpy
import MDAnalysis


# FILE VARIABLES
configFile = 'forcematchFiles3.config'
psf = None
pdb = None
forceDcd = None
coordDcd = None
param = None
force = None
coord = None

debug = False



# DEFUALT GLOBAL VARIABLES
rMin = 0
rMax = 25
binSize = 0.1
binCount = 0
elecPermittivity = 8.854e-12      # electrical permittivity of vacuum C^2/Jm
epsilon = []        # list of epsilon values for all non-solvent atoms
lj_rMin = []          # list of rMin values for all non-solvent atoms
exTs = 0            # example time step used for log purposes

# PLOT DATA VARIABLE
plots = []



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
                print('PSF File: {}'.format(psf))
        elif line == 'END CONFIG FILE\n':
            print('No PSF file found in config.')
            break

# Get name of PDB file from config file
def getPdb(cF):
    global pdb
    txt = open(cF, 'r')
    while pdb is None:
        line = txt.next()
        if line == 'PDB FILE:\n':
            pdb = txt.next()[:-1]
            if debug:
                print("PDB File: {}".format(pdb))
            elif line == 'END CONFIG FILE\n':
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
                print('Force DCD files: {}'.format(forceDcd))

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
                    print('NO FORCE DCD FILES FOUND IN CONFIG')
                coordDcd.append(line[:-1])
                line = txt.next()
            if debug:
                print('Coordinate DCD files: {}'.format(coordDcd))

# Get name of the parameter file from config file
def getParam(cF):

    global param
    txt = open(cF, 'r')
    while param is None:
        line = txt.next()
        if line == 'PARAMETER FILE:\n':
            line = txt.next()
            while line != '\n':
                if line == 'END CONFIG\n':
                    print('NO PARAMETER FILE FOUND IN CONFIG')
                param = line[:-1]
                line = txt.next()
            if debug:
                print('Parameter file: {}\n'.format(param))

# Set coordinate max/min and binsize
def getCoordBounds(cF, pdb):
    txt = open(cF,'r')
    txt2 = open(pdb, 'r')
    line = txt.next()
    dims = txt2.next().split()
    dims = dims[1], dims[2], dims[3]
    length1 = len("MIN DISTANCE: ")
    length2 = len("MAX DISTANCE: ")
    length3 = len("BIN SIZE: ")

    global rMin, rMax, binSize, binCount

    # scan config file for coord and bin values
    while line != "END CONFIG\n":
        line = txt.next()
        if line[:length1] == "MIN DISTANCE: ":
            rem = -1 * (len(line) - length1)
            rMin = int(line[rem:-1])
        elif line[:length2] == "MAX DISTANCE: ":
            rem = -1 * (len(line) - length2)
            rMax = int(line[rem:-1])
        elif line[:length3] == "BIN SIZE: ":
            rem = -1 * (len(line) - length3)
            binSize = float(line[rem:-1])

    binCount = int((rMax - rMin)/binSize)

    if debug:
        print "--Dimensions of System--"
        print "\tTotal System Dimensions: {} A x {} A x {} A".format(dims[0], dims[1], dims[2])
        print "\tMin Interparticle Distance Considered: {} A".format(rMin)
        print "\tMax Interparticle Distance Considered: {} A".format(rMax)
        print "\tBin Size Used: {}".format(binSize)
        print "\tBin Count: {}".format(binCount)

# Define subset of data without solvent
def parseWater():
    # select all atoms that are not water or hydrogen
    if debug:
        print("\n--Reading DCD Data--\n\t- Parsing out WATER...\n")
    global ionsCoord
    global ionsForce
    ionsCoord = coord.selectAtoms("not name H1 and not name H2 and not name OH2")
    ionsForce = force.selectAtoms("not name H1 and not name H2 and not name OH2")

# Initialize MD Analysis
def initMDA():
    # start MDAnalysis with a universal
    # force universe
    global force, coord, debug

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
        global ionsCoord, ionsForce, exTs
        print "--Simulation Log Info--"
        # print list of atoms being force matched
        print "DCD coord universe after course-graining:", len(ionsCoord), "atom(s)"
        for i in range(0, len(ionsCoord)):
            print "\t", ionsCoord[i]
        print "Force DCD universe after course-graining:", len(ionsForce), "atom(s)"
        for i in range(0, len(ionsForce)):
            print "\t", ionsForce[i]

        # define example timestep
        exTs = len(coord.trajectory) / 2
        # some general log info
        print "\nNumber of time steps in coordinate trajectory:", len(coord.trajectory)
        print "Number of time steps in force trajectory:", len(force.trajectory)

        if len(coord.trajectory) == len(force.trajectory):
            print("Coord and Force time step counts MATCH\n")
        else:
            print("Coord and Force time step counts DO NOT MATCH\nCheck .config file for incorrect .dcd or .force.dcd file names.")

# Iterate through all pairs of particles in all simulations,
#    identifying each pair of particles, performing computations,
#    and storing the results in a data set
def iterate():
    global plots
    if debug:
        print "-- Iterating through all particle pairs in first time step to establish pair types"
    for ts in coord.trajectory:                 # Iterate through all time steps
        for a in ionsCoord:                     # First particle
            for b in ionsCoord:                 # Second particle
                if a.number != b.number:        # Ensure particles don't have same index
                    if ts.frame == 1:               # Determine particle pairs
                        if debug:
                            print "  Identified a {}-{} pair.".format(a.name, b.name)
                        if pairIsUnique(a, b):
                                plots.append("{} {}".format(a.name, b.name))
                                plots.append(buildBlankDataSet(a, b))
                                computeCoordData(a, b)
                                computeForceData(a, b)
                        else:
                            computeCoordData(a, b)
                            computeForceData(a, b)

                    else:
                        computeCoordData(a, b)
                        computeForceData(a, b)

        if ts.frame == exTs:                        # Print a sample of the data we've computed
            exampleTimestepDebug()

    averageAll()
    zeroToNan()

# Identifies current particle pair, returns True if pair is unique (doesn't yet exist in plots array), False otherwise
def pairIsUnique(a, b):
    pair = "{} {}".format(a.name, b.name)
    pairFlipped = "{} {}".format(b.name, a.name)
    if pair in plots:
        print "\t{}-{} data set found, submitting new measurements.".format(a.name, b.name)
        return False
    elif pairFlipped in plots:
        print "\t{}-{} data set found, submitting new measurements.".format(b.name, a.name)
        return False
    else:
        if debug:
            print "\tBuilding new data set for all {}-{} pairs. Including new measurements.".format(a.name, b.name)
        return True

# Returns the index of the data set of a given pair of atoms, assuming it exists in the plots array
def findPair(a, b):
    pair = "{} {}".format(a.name, b.name)
    pairFlipped = "{} {}".format(b.name, a.name)
    if pair in plots:
        return plots.index(pair) + 1
    elif pairFlipped in plots:
        return plots.index(pairFlipped) + 1

# Builds a blank data set in which to accumulate data
def buildBlankDataSet(a, b):
    coul = numpy.zeros(binCount)
    lJ = numpy.zeros(binCount)
    totForce = numpy.zeros(binCount)
    dataCounts = numpy.zeros(binCount)  # Keep track of the number of measurements made
                                        # in each bin, for averaging purposes

    dataSet = [coul, lJ, totForce, dataCounts]
    return dataSet

# Performs a series of coordinate-related computations on one pair of particles for one time step and returns the data
def computeCoordData(a, b):
    global plots
    magR = computeMagR(a, b)            # Determine distance between particles
    if rMin <= magR < rMax:             # If distance is within specified range
        binNo = int(magR / binSize)     #   determine the appropriate bin number
        dataSet = findPair(a, b)
        if binNo < binCount:
            plots[dataSet][0][binNo] += computeCoulombic(a, b, magR)
            plots[dataSet][1][binNo] += computeLJ(a, b, magR)
            plots[dataSet][len(plots[dataSet]) - 1][binNo] += 1     # Increase measurement count by 1

# Takes atoms in coord file and returns corresponding atoms in force file
def findAtomsForce(a, b):
    particleA = None
    particleB = None

    for i in ionsForce:
        if i.number == a.number:
            particleA = i
        elif i.number == b.number:
            particleB = i
    return [particleA, particleB]

# Performs a series of force computations on one pair of particles for one time step
def computeForceData(a, b):
    global plots
    magR = computeMagR(a, b)
    if rMin <= magR < rMax:
        binNo = int(magR / binSize)
        if binNo < binCount:
            dataSet = findPair(a, b)
            forceAtoms = findAtomsForce(a, b)
            plots[dataSet][2][binNo] += computeTotalForce(magR, a, b, forceAtoms[0], forceAtoms[1])

# Converts all running sums to averages
def averageAll():
    global plots
    for i in range(0, len(plots)):                  # Sort through all particle pairings
        if i % 2 == 1:
            for f in range(0, len(plots[i])-1):
                for c in range(0, len(plots[i][f])):   # Sort through all data sets
                    if plots[i][-1:][0][c] != 0:
                        plots[i][f][c] /= plots[i][-1:][0][c]  # Divide running sum by measurement count

# Sets all uncomputed zeros to NaN
def zeroToNan():
    global plots
    plotsLength = len(plots)
    setLength = len(plots[1])
    subsetLength = len(plots[1][0])
    for set in range(1, plotsLength):
        if set % 2 == 1:        # Skip over name elements, go into data elements
            for e in range(0, setLength - 1):       # Iterate through data sets
                for m in range(0, subsetLength - 1):       # Iterate through measurement counts
                    if plots[set][setLength - 1][m] == 0:          # If there are zero measurements taken for a bin...
                        for q in range(0, setLength - 1):
                            plots[set][q][m] = numpy.nan           # Set that bin = NaN for all subsets

# Extract LJ parameters from param file
def defineEpsilonSigma(paramFile):
    global epsilon, lj_rMin
    print "-- Obtaining epsilon/(rMin/2) values from parameter file."
    txt = open(paramFile, 'r')
    line = txt.next().split()
    while line[0] != 'END':
        for a in ionsCoord:
            if a.name == line[0]:
                if a.name not in epsilon:
                    epsilon.append(a.name)
                    epsilon.append(line[2])
                    lj_rMin.append(a.name)
                    lj_rMin.append(line[3])
                    print "{} Epsilon: {}\trMin/2: {}".format(a.name, epsilon[-1:][0], lj_rMin[-1:][0])
        line = txt.next().split()
        while (len(line) == 0):
            line = txt.next().split()
    print "\n"

# Use coord dcd file to determine magnitude of distance R between two particles
def computeMagR(a, b):
    global magR
    r = a.position - b.position
    return numpy.linalg.norm(r)

# Use coord dcd file to determine coulombic potential of an atom in eV
def computeCoulombic(a, b, magR):
    global elecPermittivity
    q1 = a.charge * 1.60217646e-19              # get particle charges (in units C)
    q2 = b.charge * 1.60217646e-19
    magR *= 1e-10       # obtain the distance we just computed (in units m)
    p = (q1 * q2)/(4 * numpy.pi * elecPermittivity * magR)
    p /= 4184          # Convert from J to kcal/mol
    p *= 6.022e23
    return p

# Use coord dcd file to determine LJ potential of an atom
def computeLJ(a, b, magR):
    epsA = 0
    epsB = 0
    rmA = 0
    rmB = 0

    # Define a single epsilon given the two particles
    for i in range(0, len(epsilon)):
        if epsilon[i] == a.name:
            epsA = float(epsilon[i + 1])
        elif epsilon[i] == b.name:
            epsB = float(epsilon[i + 1])
    eps = numpy.sqrt(epsA * epsB)

    # Define a single rMin given the two particles
    for i in range(0, len(lj_rMin)):
        if lj_rMin[i] == a.name:
            rmA = float(lj_rMin[i + 1])
        elif lj_rMin[i] == b.name:
            rmB = float(lj_rMin[i + 1])
    rm = (rmA / 2) + (rmB / 2)

    h = rm / magR
    s = numpy.power(h, 12) - (2 * numpy.power(h, 6))
    lj = eps * s
    lj /= 4184
    lj *= 6.022e23

    #print " LJ ---- {} {} {} {} {}".format(eps, rm, h, s, lj)

    return lj

# Use force dcd file to examine total force interactions on particles
def computeTotalForce(magR, a, b, fA, fB):
    forceA = fA.position
    forceB = fB.position

    r = a.position - b.position
    rHat = r / magR

    avgForce = (forceA - forceB) / 2
    MagAvgF = numpy.linalg.norm(numpy.dot(avgForce, rHat))
    return MagAvgF

# Print debug info at example timestep
def exampleTimestepDebug():
    global plots
    if debug:
        print "\n-- Sample of Data"
        print "  Time Step: {}".format(exTs)
        print "  Collecting data on {} unique paring(s) of particles:".format(len(plots)/2)
        for i in range (0, len(plots)):
            if i % 2 == 0:
                print "\t{}".format(plots[i])
        print "\n\tThis process takes some time, please wait."

# Draw subplots for each data set
def plotAllData(color='r'):
    global plots
    plt.close()

    for pair in range(0, len(plots)):
        if pair % 2 == 0:
            xAxis = numpy.arange(len(plots[pair+1][0])) * binSize

            # Plot Coulombic Potential and regression line
            f, (ax0, ax1, ax2) = plt.subplots(3, sharex=True)
            ax0.scatter(xAxis, plots[pair+1][0], c=color, s=15)
            ax0.set_title("Coulombic Potential")
            ax0.set_ylabel("Energy (kcal/mol)", fontsize=10)
            coef0 = numpy.polyfit(xAxis, plots[pair+1][0], 3)
            poly0 = numpy.poly1d(coef0)
            ax0fit = poly0(xAxis)
            ax0.plot(ax0fit)

            # Plot Lennard-Jones Potential and regression line
            ax1.scatter(xAxis, plots[pair+1][1], c=color, s=15)
            ax1.set_title("Lennard-Jones Potential")
            ax1.set_ylabel("Energy (kcal/mol)", fontsize=10)
            coef1 = numpy.polyfit(xAxis, plots[pair+1][1], 5)
            poly1 = numpy.poly1d(coef1)
            ax1fit = poly1(xAxis)
            ax1.plot(ax1fit)

            # Plot Total Force and regression line
            ax2.scatter(xAxis, plots[pair+1][2], c=color, s=15)
            ax2.set_title("Total Force")
            ax2.set_ylabel("Avg Force", fontsize=10)
            coef2 = numpy.polyfit(xAxis, plots[pair+1][2], 5)
            poly2 = numpy.poly1d(coef2)
            ax2fit = poly2(xAxis)
            ax2.plot(ax2fit)

            plt.xlabel("Intermolecular distance (Angstroms)", fontsize=10)
            plt.suptitle("{} Interactions".format(plots[pair]))

    plt.show()
    '''
    # create figure with subplots
    global numPlots
    f, (ax0, ax1, ax2) = plt.subplots(numPlots, sharex=True)
    f.suptitle("MD Analysis")

    ax0.scatter(x0, y0, c=color, s=15)
    ax0.set_title("Coulombic Potential Energy")
    coef0 = numpy.polyfit(x0, y0, fitDegree)
    poly0 = numpy.poly1d(coef0)
    ys0 = poly0(x0)
    ax0.plot(x0, ys0)

    ax1.scatter(x1, y1, c=color, s=15)
    ax1.set_title("Lennard-Jones Potential Energy")
    coef1 = numpy.polyfit(x1, y1, fitDegree)
    poly1 = numpy.poly1d(coef1)
    ys1 = poly1(x1)
    ax1.plot(x1, ys1)

    ax2.scatter(x2, y2, c=color, s=15)
    ax2.set_title("Total Force")
    coef2 = numpy.polyfit(x2, y2, fitDegree)
    poly2 = numpy.poly1d(coef2)
    ys2 = poly2(x2)
    ax2.plot(x2, ys2)

    plt.xlabel("Intermolecular distance", fontsize=10)
    plt.show()
    '''

# main program
def main():
    # access global var for config file
    global configFile, pdb

    # Read config setting for debug mode
    setDebug(configFile)

    # Get name of PSF file from config file
    getPsf(configFile)

    # Get name of PDB file from config file
    getPdb(configFile)

    # Get names of Force DCD files from config file
    getForceDCDs(configFile)

    # Get names of Coord DCD files from config file
    getCoordDCDs(configFile)

    # Get name of Parameter file from config file
    getParam(configFile)

    # Define coordinate min/max and bin size
    getCoordBounds(configFile, pdb)

    # Initialize MD Analysis
    initMDA()

    # Define epsilon and sigma values for particles of interest
    defineEpsilonSigma(param)

    # Iterate over time steps, and perform MD calculations
    iterate()

    # Generate figures and plots
    plotAllData()

# Main program code
main()
