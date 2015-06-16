__author__ = 'Greg Poisson'

import matplotlib.pyplot as plt
import numpy
import MDAnalysis


# FILE VARIABLES
configFile = 'forcematchFiles2.config'
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
elecPermittivity = 8.854e-12      # electrical permittivity of vacuum C^2/Jm
epsilon = []        # list of epsilon values for all non-solvent atoms
sigma = []          # list of sigma values for all non-solvent atoms
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

    global rMin
    global rMax

    # scan config file for coord and bin values
    while line != "END CONFIG\n":
        line = txt.next()
        if line[:length1] == "MIN DISTANCE: ":
            rem = -1 * (len(line) - length1)
            rMin = int(line[rem:-1])
        elif line[:length2] == "MAX DISTANCE: ":
            rem = -1 * (len(line) - length2)
            rMax = int(line[rem:-1])

    if debug:
        print "--Dimensions of System--"
        print "\tTotal System Dimensions: {} A x {} A x {} A".format(dims[0], dims[1], dims[2])
        print "\tMin Interparticle Distance Considered: {} A".format(rMin)
        print "\tMax Interparticle Distance Considered: {} A".format(rMax)

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

######  Structure of computed data

#    Using NaCl simulation as an example:

#                           Na-Na                         Na-Cl                            Cl-Cl
#    "plots"      |-----------|---------------|-------------|---------------|----------------|---
#              SOD SOD   /  /   \   \      SOD CLA     /  /   \   \      CLA CLA        /  /   \   \
#                        |  |   |   |                  |  |   |   |                     |  |   |   |
#                        R Coul LJ  TotForce           R Coul LJ  TotForce              R Coul LJ  TotForce


# Iterate through all pairs of particles in all simulations,
#    identifying each pair of particles, performing computations,
#    and storing the results in a data set
def iterate():
    global plots
    if debug:
        print "-- Iterating through all particle pairs in first time step"
    for ts in coord.trajectory:                 # Iterate through all time steps
        for a in ionsCoord:                     # First particle
            for b in ionsCoord:                 # Second particle
                if a.number != b.number:
                    if ts.frame == 1:               # Determine particle pairs
                        if debug:
                            print "  Identified a {}-{} pair.".format(a.name, b.name)
                        if pairIsUnique(a, b):
                                plots.append("{} {}".format(a.name, b.name))
                                plots.append(assembleDataSet(a, b))
                        else:
                            addData(assembleDataSet(a, b), plots[findPair(a, b)])
                    else:
                        addData(assembleDataSet(a, b), plots[findPair(a, b)])
        if ts.frame == exTs:                        # Print a sample of the data we've computed
            exampleTimestepDebug(a, b)

# Identifies current particle pair, returns True if pair is unique (doesn't yet exist in plots array), False otherwise
def pairIsUnique(a, b):
    pair = "{} {}".format(a.name, b.name)
    pairFlipped = "{} {}".format(b.name, a.name)
    if pair in plots:
        print "\t{}-{} data set already found.".format(a.name, b.name)
        return False
    elif pairFlipped in plots:
        print "\t{}-{} data set already found.".format(a.name, b.name)
        return False
    else:
        if debug:
            print "\tBuilding new data set for all {}-{} pairs.".format(a.name, b.name)
        return True

# Returns the index of the data set of a given pair of atoms, assuming it exists in the plots array
def findPair(a, b):
    pair = "{} {}".format(a.name, b.name)
    pairFlipped = "{} {}".format(b.name, a.name)
    if pair in plots:
        return plots.index(pair) + 1
    elif pairFlipped in plots:
        return plots.index(pairFlipped) + 1

# Adds a set of data to a pre-existing data set
def addData(newData, currentData):
    for i in range (0, len(currentData)):
        currentData[i].append(newData[i][0])

# Performs a series of computations on one pair of particles for one time step and returns the data
def assembleDataSet(a, b):
    magR = []
    coul = []
    lJ = []
    totForce = []
    dataSet = [magR, coul, lJ, totForce]
    dataSet[0].append(computeMagR(a, b))              # Compute distance between particles
    dataSet[1].append(computeCoulombic(a, b, dataSet[0][-1:][0]))         # Compute Coulombic potential of both particles
    dataSet[2].append(computeLJ())                    # Compute Lennard-Jones potential of both particles
    dataSet[3].append(computeTotalForce())            # Compute total force on both particles
    return dataSet

# Extract LJ parameters from param file
def defineEpsilonSigma(paramFile):
    global epsilon, sigma
    txt = open(paramFile, 'r')
    line = txt.next().split()
    while ((len(epsilon) < len(ionsCoord)) | (len(sigma) < len(ionsCoord))) & (line[0] != 'END'):
        for i in range(len(epsilon), len(ionsCoord)):
            if ionsCoord[i].name == line[0]:
                # add epsilons in the order atoms are stored in ionsCoord
                if i == (len(epsilon)):
                    epsilon.append(float(line[2]))
                    sigma.append(float(line[3]))

        line = txt.next()
        while line == '\n':
            line = txt.next()
        line = line.split()

# Use coord dcd file to determine magnitude of distance R between two particles
def computeMagR(a, b):
    global magR
    r = a.position - b.position
    return numpy.linalg.norm(r)

# Use coord dcd file to determine coulombic potential of an atom in eV
def computeCoulombic(a, b, r):
    global magR, elecPermittivity
    q1 = a.charge * 1.60217646e-19              # get particle charges (in units C)
    q2 = b.charge * 1.60217646e-19
    r *= 1e-10       # obtain the distance we just computed (in units m)
    p = (q1 * q2)/(4 * numpy.pi * elecPermittivity * r)
    p *= 4184          # Convert from J to kcal/mol
    return p

# Use coord dcd file to determine LJ potential of an atom
def computeLJ():
    return 0

# Use force dcd file to examine total force interactions on particles
def computeTotalForce():
    return 0

# Pring debug info at example timestep
def exampleTimestepDebug(a, b):
    global plots
    if debug:
        print "\n-- Sample of Data"
        print "  Time Step: {}".format(exTs)
        print "  Collecting data on {} unique paring(s) of particles:".format(len(plots)/2)
        for i in range (0, len(plots)):
            if i % 2 == 0:
                print "\t{}".format(plots[i])
            else:
                print "\t Current size: {} interactions".format(len(plots[i][0]))
                print "\t Sample distance: {} A".format(plots[i][0][0])
                print "\t Sample Coulombic {} kcal/mol".format(plots[i][1][0])


# Draw subplots for each data set
def plotAllData(fitDegree=9, color='r'):
    global plots
    plt.close()
    for pair in range(0, len(plots)):
        if pair % 2 == 0:
            f, (ax0, ax1) = plt.subplots(2, sharex=True)
            ax0.scatter(plots[pair+1][0], plots[pair+1][1], c=color, s=15)
            ax0.axis([rMin, rMax, -1e-15, 1e-15])

            ax1.scatter(plots[pair+1][0], plots[pair+1][2], c=color, s=15)

            plt.xlabel("Intermolecular distance", fontsize=10)
            plt.suptitle("MD Analysis")
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
