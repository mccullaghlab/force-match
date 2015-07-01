__author__ = 'Greg Poisson'

import matplotlib.pyplot as plt
import numpy
import MDAnalysis
import time


# FILE VARIABLES
configFile = 'forcematchFiles.config'
psf = None
pdb = None
forceDcd = None
coordDcd = None
param = None
force = None
coord = None
temperature = None
dims = None

debug = False
junkCounter = 0     # counter used for debugging


# DEFUALT GLOBAL VARIABLES
rMin = 0
rMax = 25
binSize = 0.1
binCount = 0
elecPermittivity = 8.854e-12      # electrical permittivity of vacuum C^2/Jm
boltzmann = 1.9872041e-3          # boltzmann constant in kcal/(K mol)
epsilon = []        # list of epsilon values for all non-solvent atoms
lj_rMin = []          # list of rMin values for all non-solvent atoms
exTs = 0            # example time step used for log purposes

# PLOT DATA VARIABLE
plots = []
plotOrder = []

titles = [
    "Coulombic Potential",
    "Lennard-Jones Potential",
    "Total Force",
    "Integrated Force",
    "Radial Distribution Frequency",
    "Free Energy",
    "Probability Density Distribution of Distance Frequency",
    "Frequency Distribution of Particle Distance"
]

yLabels = [
    "Energy (kcal/mol)",
    "Energy (kcal/mol)",
    "Avg Force",
    "Free Energy",
    "g(r)",
    "Free Energy \n(kcal/mol)",
    "Probability",
    "Occurrences"
]


'''
    DATA ORGANIZATION:
        Plots -- List of name / dataset pairs
        Data -- Ordered list of data sets. Real-time data populates the front of the list,
                while post-process data populates the end of the list.
                Plots may be drawn in any order, but they are built and stored as follows:
        Data set order:

            (Active measurements:)
            [0] -- Coulombic Potential
            [1] -- Lennard-Jones Potential
            [2] -- Total Force

            (Post-process measurements:)
            [-5] -- Integrated Force
            [-4] -- Radial Distribution Frequency
            [-3] -- Free Energy
            [-2] -- Probability Density Distribution of Distance Frequency
            [-1] -- Frequency of Particle Distance

            EX: Plots-  0       1       2       3       4       5
                     SOD-SOD   Data  SOD-CLA   Data  CLA-CLA   Data
                             [Coul]          [Coul]          [Coul]
                             [LJ]            [LJ]            [LJ]
                             [TF]            [TF]            [TF]
                             [IF]            [IF]            [IF]
                             [RDF]           [RDF]           [RDF]
                             [FE]            [FE]            [FE]
                             [Prob Dist]     [Prob Dist]     [Prob Dist]
                             [Freq]          [Freq]          [Freq]

'''



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
    global dims
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

# Get temperature of system from config file
def getTemp(cF):
    global temperature
    txt = open(cF, 'r')
    while temperature is None:
        line = txt.next().split()
        if len(line) > 1:
            if (line[0] == "SYSTEM") & (line[1] == "TEMPERATURE:"):
                temperature = float(line[2])
                if debug:
                    print "\tSystem temperature: {} K".format(temperature)
            elif line[0] == "END":
                break

# Define which plots to draw, and in which order
def getPlotOrder(cF):
    global plotOrder
    txt = open(cF, 'r')
    line = txt.next().split()
    while len(line) == 0:
        line = txt.next().split()
    read = False
    while line[0] != "END":
        if (line[0] == "PLOTS") & (line[1] == "PLOT") & (line[2] == "ORDER"):   # Ignore config lines before plot order
            read = True
            line = txt.next().split()
        if read:                        # Begin reading plot order section of config
            count = int(line[-1:][0])   # Read user-defined number
            plotOrder.append(count)
            line = txt.next().split()
            if len(line) == 0:
                break
        else:
            line = txt.next().split()
            while len(line) < 3:
                line = txt.next().split()

            if line[0] == "END":
                break

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
    # start MDAnalysis with a universal force universe
    global force, coord, dims, rMax, binCount, debug

    m = False       # Swtiched to True if user requests an rMax value greater than the system allows
    force = MDAnalysis.Universe(psf, forceDcd)
    # coordinate universe
    coord = MDAnalysis.Universe(psf, coordDcd)
    dims = [coord.dimensions[0], coord.dimensions[1], coord.dimensions[2]]
    rMaxLimit = numpy.sqrt((dims[0]**2) + (dims[1]**2) + (dims[2]**2))
    if rMax > rMaxLimit:
        rMax = rMaxLimit
        m = True

    binCount = int((rMax - rMin)/binSize)

    if debug:
        print "--Dimensions of System--"
        print "\tTotal System Dimensions: {} A x {} A x {} A".format(dims[0], dims[1], dims[2])
        print "\tMin Interparticle Distance Considered: {} A".format(rMin)
        if m:
            print "\tMax Interparticle Distance Considered: {} A\tThe requested rMax value was bigger than the simulation box size permits, and was truncated".format(rMax)
        else:
            print "\tMax Interparticle Distance Considered: {} A".format(rMax)
        print "\tBin Size Used: {}".format(binSize)
        print "\tBin Count: {}".format(binCount)

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

    postProcess()

# Identifies current particle pair, returns True if pair is unique (doesn't yet exist in plots array), False otherwise
def pairIsUnique(a, b):
    pair = "{} {}".format(a.name, b.name)
    pairFlipped = "{} {}".format(b.name, a.name)
    if pair in plots:
        if debug:
            print "\t{}-{} data set found, submitting new measurements.".format(a.name, b.name)
        return False
    elif pairFlipped in plots:
        if debug:
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
    integF = numpy.zeros(binCount)

    freeEnergy = numpy.zeros(binCount)
    rdf = numpy.zeros(binCount)
    probDensity = numpy.zeros(binCount) # Probability density of data counts
    dataCounts = numpy.zeros(binCount)  # Keep track of the number of measurements made
                                        # in each bin, for averaging purposes

    dataSet = [coul, lJ, totForce, integF, freeEnergy, rdf, probDensity, dataCounts]
    return dataSet

# Performs a series of coordinate-related computations on one pair of particles for one time step and returns the data
def computeCoordData(a, b):
    global plots, plotOrder
    magR = computeMagR(a, b)            # Determine distance between particles
    if rMin <= magR < rMax:             # If distance is within specified range
        binNo = int(magR / binSize)     #   determine the appropriate bin number
        dataSet = findPair(a, b)
        if binNo < binCount:
            if plotOrder[0] > 0:        # If user chose to plot coulombic data, compute the data
                #print "Coulombic plot: plots[{}][0][{}]".format(dataSet, binNo)
                plots[dataSet][0][binNo] += computeCoulombic(a, b, magR)
            if plotOrder[1] > 0:        # if user chose to plot lennard-jones data, compute the data
                #print "LJ plot: plots[{}][1][{}]".format(dataSet, binNo)
                plots[dataSet][1][binNo] += computeLJ(a, b, magR)

            #print "Incrementing count: plots[{}][{}][{}]".format(dataSet, len(plots[dataSet])-1, binNo)
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

    if (plotOrder[2] > 0) | (plotOrder[3] > 0):
        magR = computeMagR(a, b)
        if rMin <= magR < rMax:
            binNo = int(magR / binSize)
            if binNo < binCount:
                dataSet = findPair(a, b)
                forceAtoms = findAtomsForce(a, b)
                plots[dataSet][2][binNo] += computeTotalForce(magR, a, b, forceAtoms[0], forceAtoms[1])

# Integrates a set of force data for a given pair of atoms
def integrateForce():
    global integF
    for a in ionsCoord:
        for b in ionsCoord:
            if a.number != b.number:
                index = findPair(a, b)
                sum = 0

                # Integrate the force data array, store in integrated force data array
                for tf in range(0, len(plots[index][2]) - 1):
                    tf = len(plots[index][2]) - 1 - tf
                    sum += plots[index][2][tf] * binSize
                    plots[index][len(plots[index])-5][tf] = sum

# Perform post-data mining calculations
def postProcess():
    averageAll()

    if plotOrder[len(plotOrder)-5] > 0:
        integrateForce()
    if (plotOrder[len(plotOrder)-4] > 0) | (plotOrder[len(plotOrder)-3] > 0):
        rdf()
    if plotOrder[len(plotOrder)-3] > 0:
        freeEnergy()
    if plotOrder[len(plotOrder)-2] > 0:
        distanceDistribution()
    # The last element in plotOrder is distance frequency, which is automatically populated

    zeroToNan()

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

# Convert distance frequency to distance probability distribution
def distanceDistribution():
    global plots
    lastSet = len(plots)                             # index of last plots[] element (a set)
    countSubsetIndex = len(plots[1])-1               # index of last subset in a set, containing
                                                     #     observation counts at each bin length
    probDistIndex = len(plots[1])-2                  # index of subset given for storing probability density data

    for set in range(0, lastSet):
        if set % 2 == 1:                             # set will reference all interaction pair datasets
            n = plots[set][countSubsetIndex]/numpy.sum(plots[set][countSubsetIndex])
            decimal = n / binSize
            plots[set][probDistIndex] = decimal

# Extract LJ parameters from param file
def defineEpsilonSigma(paramFile):
    global epsilon, lj_rMin
    if debug:
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
                    if debug:
                        print "{} Epsilon: {}\trMin/2: {}".format(a.name, epsilon[-1:][0], lj_rMin[-1:][0])
        line = txt.next().split()
        while (len(line) == 0):
            line = txt.next().split()
    if debug:
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
    global junkCounter
    epsA = None
    epsB = None
    rmA = None
    rmB = None

    # Define a single epsilon given the two particles
    for i in range(0, len(epsilon)):
        if (i % 2 == 0) & (epsA == None) & (epsilon[i] == a.name):
            epsA = float(epsilon[i + 1])
    for i in range(0, len(epsilon)):
        if (i % 2 == 0) & (epsB == None) & (epsilon[i] == b.name):
            epsB = float(epsilon[i + 1])

    if epsA == None:
        epsA = 0
    if epsB == None:
        epsB = 0
    eps = numpy.sqrt(epsA * epsB)       # Epsilon, the well depth

    # Define a single rMin given the two particles
    for i in range(0, len(lj_rMin)):
        if (rmA == None) & (lj_rMin[i] == a.name):
            rmA = float(lj_rMin[i + 1])     # rMin/2 value for particle A
    for i in range(0, len(lj_rMin)):
        if (rmB == None) & (lj_rMin[i] == b.name):
            rmB = float(lj_rMin[i + 1])     # rMin/2 value for particle B

    # Catch cases to prevent compiler complaining
    if rmA == None:
        rmA = 0
    if rmB == None:
        rmB = 0

    rm = (rmA) + (rmB)              # rMin,i,j, value of V at rm is epsilon

    lj = eps * (numpy.power((rm/magR), 12) - (2 * numpy.power((rm/magR), 6)))

    junkCounter = 1
    return lj

# Use force dcd file to examine total force interactions on particles
def computeTotalForce(magR, a, b, fA, fB):
    forceA = fA.position
    forceB = fB.position

    r = a.position - b.position
    rHat = r / magR

    avgForce = (forceA - forceB) / 2
    MagAvgF = numpy.dot(avgForce, rHat)
    return MagAvgF

'''
    RDF function should approach 1 as r >> 0.
    Currently, it appraches 0 as r >> 0.

'''
# Determine radial distribution frequency data
def rdf():
    for set in range(0, len(plots)):
        if set % 2 == 1:
            for bin in range(1, binCount):
                radius = bin*binSize
                deltaR = binSize

                volume = 4 * numpy.pi * radius**2 * deltaR

                N = float(len(plots)/2)
                V = float(dims[0]) * float(dims[1]) * float(dims[2])
                rho = float(N/V)

                Pr = plots[set][len(plots[set])-1][bin]/numpy.sum(plots[set][len(plots[set])-1])     # counts(r) / r
                #g = Pr / (volume/V)
                hr = plots[set][len(plots[set])-1][bin]
                g = hr / (volume * rho * len(coord.trajectory))
                plots[set][len(plots[set])-4][bin] = g

# Compute Free Energy data from Radial Distribution Data
def freeEnergy():
    for set in range(0, len(plots)):
        if set % 2 == 1:
            for bin in range(0, binCount):
                rdf = plots[set][len(plots[set])-4][bin]
                if rdf != 0.0:
                    logGr = numpy.log(rdf)         # Get probability and take log
                    fe = -boltzmann * temperature * logGr          # Compute free energy
                    plots[set][len(plots[set])-3][bin] = fe

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
            axes = []       # List of the requested plot order

            for ax in range(1, len(plotOrder)+1):
                # go through plotOrder and arrange these plots.
                if ax in plotOrder:
                    axes.append(plotOrder.index(ax))


            # Create plots and labels
            if len(axes) == 0:        # If zero plots are requested, end the program
                break
            elif len(axes) > 1:       # Multiple plots
                f, yAxes = plt.subplots(len(axes), sharex=True)

                # Draw all requested plots
                for i in range(0, len(axes)):
                    yAxes[i].scatter(xAxis, plots[pair+1][axes[i]], c=color, s=15)
                    yAxes[i].plot(xAxis, plots[pair+1][axes[i]])
                    yAxes[i].grid(True)
                    yAxes[i].set_title(titles[axes[i]])
                    yAxes[i].set_ylabel(yLabels[axes[i]], fontsize=10)
                    #coef0 = numpy.polyfit(xAxis, plots[pair+1][0], 3)
                    #poly0 = numpy.poly1d(coef0)
                    #ax0fit = poly0(xAxis)
                    #ax0.plot(ax0fit)

                plt.suptitle("{} Interactions".format(plots[pair]))

            elif len(axes) == 1:           # Single plot
                plt.scatter(xAxis, plots[pair+1][axes[0]], c=color, s=15)
                plt.plot(xAxis, plots[pair+1][axes[0]])
                plt.grid(True)
                plt.suptitle("{} -- {} Interactions".format(titles[axes[0]], plots[pair]))
                plt.ylabel(yLabels[axes[0]], fontsize=10)

            plt.xlabel("Intermolecular Distance (Angstroms)", fontsize=10)

    plt.show()

# main program
def main():

    # access global var for config file
    global configFile, pdb
    start = time.time()

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

    # Get temperature from config file
    getTemp(configFile)

    # Get plot order from config file
    getPlotOrder(configFile)

    # Initialize MD Analysis
    initMDA()

    # Define epsilon and sigma values for particles of interest
    defineEpsilonSigma(param)

    # Iterate over time steps, and perform MD calculations
    iterate()

    end = time.time()
    t = end - start
    print "\nTotal running time: {} sec".format(t)

    # Generate figures and plots
    plotAllData()


# Main program code
main()
