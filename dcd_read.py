#!/Users/martinmccullagh/anaconda/bin/python


import scipy
import matplotlib.pyplot as plt
import numpy
import MDAnalysis

configFile = 'forcematchFiles.config'
psf = None
pdb = None
forceDcd = None
coordDcd = None
param = None
force = None
coord = None

debug = False

# default vars
rMin = 0
rMax = 25
rBinSize = 0.1
nBins = int((rMax - rMin) / rBinSize)
elecPermittivity = 710e-12      # electrical permittivity of water
epsilon = []        # list of epsilon values for all non-solvent atoms
sigma = []          # list of sigma values for all non-solvent atoms
magR = 0    # magnitude of distance between particles
q1 = 0      # charge of particle (C)
q2 = 0
rHat = [0, 0, 0]    # normal vector in R direction
coulombic = None
numPlots = 3

# average force array
fMagAvg = numpy.zeros(nBins, dtype=float)
# count array
fMagCount = numpy.zeros(nBins, dtype=int)
# coulombic potentials
coul = numpy.zeros(nBins, dtype=float)
# lennard jones potentials
lJ = numpy.zeros(nBins, dtype=float)
# total forces from force dcd
totForce = numpy.zeros(nBins, dtype=float)
totForceCount = numpy.zeros(nBins, dtype=float)



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
# The PDB file itself isn't yet utilized; this is here
#    for ease of adding pdb support later if needed.
def getPdb(cF):
    global pdb
    txt = open(cF, 'r')
    while psf is None:
        line = txt.next()
        if line == 'PDB FILE:\n':
            pdb = line.next()[:-1]
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

    global rHat, magR, q1, q2, param
    # choose a time step to be referenced as an example when debugging
    exampleTs = int(len(coord.trajectory) / 2)
    # junk variable for log purposes
    j = 0

    # Debugging reference time step
    if debug:
        print("Referencing time step {} for logging purposes.\n".format(exampleTs))

    # initialize distance vector
    r = [0., 0., 0.]
    # normal vector in distance direction
    rHat = 0
    # magnitude of distance vector
    magR = 0

    # get box dimensions (does not change between time steps)
    box = coord.trajectory.ts.dimensions[:3]

    # iterate over coordinate trajectory
    for ts in coord.trajectory:
        # iterate through all combinations of coord ion pairs
        for a in range(0, len(ionsCoord)):
            for b in range(a, len(ionsCoord)):
                # skip computing interactions of an atom with itself
                if a != b:
                    # initial data collecting step for each time step
                    debugAndComputeDistanceData(ts.frame, exampleTs, j, r, a, b, box)
                    # compute Coulombic interaction data, store in var "coul"
                    computeCoulombic(q1, q2, ts.frame, exampleTs)
                    # compute Lennard-Jones interaction data, store in var "lJ"
                    computeLennardJones(param, ts.frame, exampleTs, a, b)
                    # compute Total Force data, store in var "totForce"
                    computeTotalForce(ts.frame, exampleTs, j, ionsForce, a, b, rHat)

    # create bins to use for x-axis
    bins = numpy.arange(nBins)
    # plot all data (currently plots coulombic data)
    # if further plots are added, remember to change global var "numPlots"
    plotAllData(bins, coul, bins, lJ, bins, totForce)

def debugAndComputeDistanceData(frame, exampleTs, j, r, a, b, box):
    global debug, magR, rHat, q1, q2
    # print debug info for example timestep
    if (frame == exampleTs) & debug & (j == 0):
        j = 1
        print(" ** Coordinate computations for timestep {}:\n".format(exampleTs))
    # compute 3D distance vector, r
    for i in range(0, 3):
        r[i] = ionsCoord[a].pos[i] - ionsCoord[b].pos[i]

        # check which image is closest
        if r[i] > box[i] / 2:
            r[i] -= box[i]
        elif r[i] < -box[i] / 2:
            r[i] += box[i]
    # compute the magnitude of the distance vector (aka the distance)
    magR = numpy.linalg.norm(r)
    # compute the normal vector in the distance direction
    rHat = r / magR

    # obtain charges of the particles
    q1 = ionsCoord[a].charge
    q2 = ionsCoord[b].charge

    # print debug info for example timestep
    if (frame == exampleTs) & debug:
        print("A: {}\nB: {}".format(ionsCoord[a].name, ionsCoord[b].name))
        print("\tAtom A coords: {}    Charge: {}".format(ionsCoord[a].position, q1))
        print("\tAtom B coords: {}    Charge: {}".format(ionsCoord[b].position, q2))
        print("\tDist from A -> B: {} Angstroms".format(magR))

# compute coulombic potential
def computeCoulombic(q1, q2, frame, exampleTs):

    # access some global vars
    global elecPermittivity, coulombic, nBins, coul, rMin, rBinSize, magR

    q1 *= 1.60217646e-19
    q2 *= 1.60217646e-19
    magRMeters = magR * 1e-10     # convert angstroms to meters
    coulombic = (q1 * q2)/(4 * numpy.pi * elecPermittivity * magRMeters)
    coulombic /= 6.022e-23          # Multiply by Avogadro
    coulombic *= 0.000239005736     # Convert to J from kCal
    magRbin = int((magR - rMin) / rBinSize)
    # sort the data into bins
    if 0 <= magRbin < nBins:
        coul[magRbin] += coulombic

    # print debug info for example timestep
    if (frame == exampleTs) & debug:
        print("\n\t** Coulomb Potential Computation\n\tV = (q1 * q2) / (4 * pi * e0 * r)")
        print("\t\te0 = {} C^2 J^-1\tr = {} m\t\tq1 = {} C\tq2 = {} C".format(elecPermittivity, magRMeters, q1, q2))
        print("\n\t\tCoulombic potential between A--B = {} C    Sorted to bin: {}\n".format(coulombic, magRbin))

# compute Lennard-Jones potential
def computeLennardJones(paramFile, frame, exampleTs, a, b):
    global lJ, magR, epsilon, sigma, rMin, rBinSize, nBins
    if (len(epsilon) == 0) | (len(sigma) == 0):
        defineEpsilonSigma(paramFile)

    lJa = lJFunction(a)
    lJb = lJFunction(b)
    lJAvg = (lJa + lJb) / 2

    magRbin = int((magR - rMin) / rBinSize)
    if 0 <= magRbin < nBins:
        lJ[magRbin] += lJAvg

    if debug & (frame == exampleTs):
        print("\n\t**Lennard-Jones Potential Computation\n")
        for i in range(0, len(ionsCoord)):
            print("\t{}\t\tEpsilon: {}\t\tSigma: {}".format(ionsCoord[i].name, epsilon[i], sigma[i]))
            print("\tlJAvg: {}".format(lJAvg))

# LJ equation
def lJFunction(i):
    global epsilon, sigma, magR
    return (4 * epsilon[i] * (numpy.power((sigma[i] / magR), 12) - numpy.power((sigma[i] / magR), 6)))

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

# Use force dcd file to examine total force interactions on particles
def computeTotalForce(frame, exampleTs, j, ionsForce, a, b, rHat):
    global totForce, totForceCount, magR

    # Now work in Force Universe
    if (frame == exampleTs) & debug & (j == 1):
        print("\n ** Force computations for timestep {}:\n".format(exampleTs))
        j = 2
    ionAForce = ionsForce[a].pos
    if (frame == exampleTs) & debug:
        print("Atom pair {}, {}".format(a, b))
        print("\tAtom A force traj: {}".format(ionAForce))
    ionBForce = ionsForce[b].pos
    if (frame == exampleTs) & debug:
        print("\tAtom B force traj: {}".format(ionBForce))

    # compute the projection of the force in the distance direction
    AMagF = numpy.linalg.norm(numpy.dot(ionAForce, rHat))
    BMagF = numpy.linalg.norm(numpy.dot(ionBForce, rHat))

    if (frame == exampleTs) & debug:
        print("\t\tMagnitude of Force of atom A: {}".format(AMagF))
        print("\t\tMagnitude of Force of atom B: {}".format(BMagF))

    # First compute array index of magR
    magRbin = int((magR - rMin) / rBinSize)
    if (frame == exampleTs) & debug:
        print("\tMagnitude of R = {} Angstroms -> Bin: {}".format(magR, magRbin))

    # add to average force array
    if 0 <= magRbin < nBins:
        totForce[magRbin] += AMagF
        totForce[magRbin] += BMagF
        # add to counts
        totForceCount[magRbin] += 2

# Draw subplots for each data set
# Check "numPlots" var in global variables
def plotAllData(x0, y0, x1, y1, x2, y2, fitDegree=9, color='r'):

    # create figure with subplots
    global numPlots
    f, (ax0, ax1, ax2) = plt.subplots(numPlots)
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

    plt.show()

# main program
def main():
    # access global var for config file
    global configFile

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
    getCoordBounds(configFile)
    # Initialize MD Analysis
    initMDA()
    # Iterate over time steps, and perform MD calculations
    iterateAndCompute()


main()
