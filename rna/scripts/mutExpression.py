#!/usr/bin/python


MIN_N = 5

PROTEINS_OF_INTEREST = None
#PROTEINS_OF_INTEREST = ['SPI1', 'SPIB', 'SPIC']
#PROTEINS_OF_INTEREST = ['CEBPA', 'CEBPB', 'CEBPD', 'CEBPE', 'CEBPG', 'CEBPZ' ]

USE_CORRECTION = True
USE_SNAP_ALIGNED_DATA = False
SNAP_ALIGNED_RNA_PATH = '../data/rna-gene-expression'
TCGA_SUMMARY_RNA_PATH = '../data/tcga-laml.rnaseq.179_v1.0_gaf2.0_rpkm_matrix.txt.tcgaID.txt'

########################################################################

import sys
import os
import glob
import xml.dom.minidom
import csv
import numpy
import scipy
import scipy.stats
import itertools
import collections

if USE_CORRECTION:
    import rpy2.robjects as robjects

global outFile
def say(s):
    outFile.write(s)
    sys.stderr.write(s)

def patientIdFromSampleId(sampleId):
    if 'TCGA' in sampleId:
        return int(sampleId.split("-")[2])
    else:
        return int(sampleId.split("-")[0])

class Patients:
    def __init__(self):
        self.map = {}

        self.readMutationData(self)
        self.readExpressionData(self)

        # Print some statistics
        has_expressions = len([p for p in self.all() if p.getExpressionData()])
        has_muts = len([p for p in self.all() if p.getMutations().hasData()])
        has_both = len(self.valid())

        say("Got expression data for %d, mutations for %d, both for %d\n" % (
            has_expressions, has_muts, has_both))

    def get(self, patientId):
        if not patientId in self.map:
            self.map[patientId] = PatientData()

        return self.map[patientId]

    def all(self):
        return sorted(self.map.values())

    def valid(self):
        return [p for p in self.all() if p.isValid()]
        
    def size(self):
        return len(self.map.keys())


    # Populates each patient record with their expression ExpressionData 
    def readExpressionData(self, patients):
        if USE_SNAP_ALIGNED_DATA:
            self.readSnapAlignedData(patients)
        else:
            self.readTcgaSummaryData(patients)

    def readTcgaSummaryData(self, patients):
        if not os.path.exists(TCGA_SUMMARY_RNA_PATH):
            sys.stderr.write("Error: TCGA summary RNA data path '%s' does not exist" % (TCGA_SUMMARY_RNA_PATH))
            sys.exit(1)

        with open(TCGA_SUMMARY_RNA_PATH, 'r') as rnaMatrix:
            patientsInFile = []

            # Read the first (header) line and construct an array of the patient records we'll insert
            # each expression record into
            header = rnaMatrix.readline()
            sampleIds = header.split()

            seenFirst = False
            for sampleId in sampleIds:
                if not seenFirst:
                    if sampleId == 'GeneID':
                        seenFirst = True
                        continue
                    else:
                        sys.stderr.write("TCGA RNA summary data file doesn't have 'GeneID' where expected")
                        sys.exit(1)

                patientId = patientIdFromSampleId(sampleId)
                patient = patients.get(patientId)
                patientsInFile.append(patient)

            # Now read the data and insert
            for line in rnaMatrix:
                cols = line.split()
                geneName = cols[0].split("|")[0]
                expressions = cols[1:]

                if geneName == '?':
                    continue

                if len(expressions) != len(patientsInFile):
                    sys.stderr.write("Whoa. Header line doesn't seem to be the same length as data lines.")
                    sys.exit(1)

                for patient, expression in zip(patientsInFile, expressions):
                    patient.getExpressionData().set(geneName, expression)


    def readSnapAlignedExpressionData(self, patients):
        if not os.path.exists(SNAP_ALIGNED_RNA_PATH):
            sys.stderr.write("Error: RNA data path '%s' does not exist" % (SNAP_ALIGNED_RNA_PATH))
            sys.exit(1)

        if not os.path.isdir(SNAP_ALIGNED_RNA_PATH):
            sys.stderr.write("Error: RNA data path '%s' is not a directory" % (SNAP_ALIGNED_RNA_PATH))
            sys.exit(1)

        globList = glob.glob(os.path.join(SNAP_ALIGNED_RNA_PATH, "*"))

        if len(globList) == 0:
            sys.stderr.write("Error: RNA directory %s is empty" % (SNAP_ALIGNED_RNA_PATH))

        sys.stderr.write("Reading expression data from %s\n" % (SNAP_ALIGNED_RNA_PATH))

        for resultDir in globList:
            studyId = os.path.basename(resultDir)

            metadataFile = os.path.join(resultDir, studyId + ".cgquery.xml")
            countFile = os.path.join(resultDir, studyId + ".snap.gene_name.counts.txt")

            # Skip directories that only have metadata, no expression data
            if not os.path.exists(countFile):
                continue

            metadata = xml.dom.minidom.parse(metadataFile)

            sampleId = metadata.getElementsByTagName("ResultSet")[0] \
                               .getElementsByTagName("Result")[0] \
                               .getElementsByTagName("legacy_sample_id")[0] \
                               .childNodes[0] \
                               .nodeValue

            patientNumber = patientIdFromSampleId(sampleId)

            self.get(patientNumber).setExpressionData(ExpressionData.getFromCountFile(countFile))

        ############
        # Correct the expression data for coverage using the PoissonSeq library in R
        # http://cran.r-project.org/web/packages/PoissonSeq/index.html

        # First, generate the canonical list of proteins for which we have expression data
        # for every patient
        sys.stderr.write("Constructing protein list\n")

        expressionProteinSet = collections.Counter()
        numValidPatients = 0
        for patient in self.valid():
            numValidPatients += 1
            expressionProteinSet.update(patient.getExpressionData().getAllProteins())

        sys.stderr.write("Finding fully measured expressed proteins\n")
        proteinList = []
        for protein in sorted(expressionProteinSet.keys()):
            if expressionProteinSet[protein] == numValidPatients:
                proteinList.append(protein)

        sys.stderr.write("%d/%d proteins have expression data for all %d valid patients\n" % (
            len(proteinList), len(expressionProteinSet), numValidPatients))

        # Construct matrix to pass to PoissonSeq -- columns are samples, rows are genes.
        # R's matrix constructor takes items in column-major order by default
        sys.stderr.write("Constructing matrix for correction\n")
        expressionList = []
        for patient in self.valid():
            for protein in proteinList:
                expressionList.append(patient.getExpressionData().get(protein))

        if USE_CORRECTION:
            # Pass expression matrix to PoissonSeq for correction and get the corrections back!
            robjects.globalenv['eMat'] = robjects.r.matrix(robjects.IntVector(expressionList), nrow=len(proteinList))
            robjects.r('library("PoissonSeq")')
            scaleVector = robjects.r('scaleVector = PS.Est.Depth(eMat)')

            s = "Relative coverages per patient:\n%s\n" % (scaleVector)
            say(s)
            sys.stderr.write(s)

            i = 0
            for patient in self.valid():
                patient.getExpressionData().correct(scaleVector[i])
                i += 1

        
            
    def allValidPatientsHaveExpressionForProtein(self, protein):
        for patient in self.valid():
            if not patient.getExpressionData().has(protein):
                return False

        return True

    def readMutationData(self, patients):
        sys.stderr.write("Reading mutation data\n")

        with open('../data/tcga-laml-genes-frankannotations.csv', 'r') as mutationsCSV:
            mutations = csv.reader(mutationsCSV, delimiter=',')
            headerSkipped = False

            for mut in mutations:
                if not headerSkipped:
                    headerSkipped = True
                    continue

                # original tcga-laml
#                gene = mut[0]
#                mutType = mut[8]
#                sampleId = mut[15]

                # frank's annotations
                gene = mut[0]
                mutType = mut[1]
                sampleId = mut[2]

                patientNumber = patientIdFromSampleId(sampleId)

                # Skip "harmless" mutations
                if mutType in ['Silent', 'Intron']:
                    continue

                self.get(patientNumber).getMutations().add(gene)

class PatientData:
    def __init__(self):
        self.mutations = Mutations()
        self.expressionData = ExpressionData()

    def isValid(self):
        return self.getMutations().hasData() and self.getExpressionData()

    def getMutations(self):
        return self.mutations

    def getExpressionData(self):
        return self.expressionData

    def setExpressionData(self, expressionData):
        self.expressionData = expressionData

class ExpressionData:
    def __init__(self):
        self.proteins = {}

    @staticmethod
    def getFromCountFile(countFilename):
        retval = ExpressionData()

        for line in open(countFilename):
            cols = line.split()
            retval.set(cols[0], cols[1])

        return retval

    def getAllProteins(self):
        return self.proteins.keys()

    def has(self, proteinName):
        return proteinName in self.proteins

    def get(self, proteinName):
        return self.proteins[proteinName]

    def set(self, proteinName, level):
        self.proteins[proteinName] = float(level)

    def correct(self, scaleFactor):
        for proteinName in self.proteins:
            self.proteins[proteinName] /= scaleFactor

class Mutations:
    def __init__(self):
        self.mutations = set()

    def add(self, gene):
        self.mutations.add(gene)

    def has(self, mutation):
        return mutation in self.mutations

    def all(self):
        return self.mutations

    def hasData(self):
        return len(self.mutations) >= 1

def meanAndStdDev(d):
    return "%8.2f +/- %8.2f" % (numpy.mean(d), numpy.std(d))

class Result:
    def __init__(self, mutationList, protein, significance, normalSet, mutatedSet):
        self.mutationList = mutationList
        self.protein = protein
        self.significance = significance
        self.normalSet = normalSet
        self.mutatedSet = mutatedSet

    def getMutationList(self):
        return "+".join(self.mutationList)

    def csv(self):
        return "%s,%s,%f,%f,%f,%f,%d,%f" % (
            self.getMutationList(),
            self.protein,
            numpy.mean(self.normalSet),
            numpy.std(self.normalSet),
            numpy.mean(self.mutatedSet),
            numpy.std(self.mutatedSet),
            len(self.mutatedSet),
            self.significance)
            
    def __repr__(self):
        return "%15s: %15s %s  -->  %s, n=%2d, p=%f" % (
            self.getMutationList(),
            self.protein,
            meanAndStdDev(self.normalSet),
            meanAndStdDev(self.mutatedSet),
            len(self.mutatedSet),
            self.significance)

def getExpressionDataForProtein(patientList, protein):
    expressionDataList = []

    hasData = False

    for patient in patientList:
        expressionData = patient.getExpressionData()

        if expressionData.has(protein):
            val = expressionData.get(protein)
            expressionDataList.append(val)

            if val > 0:
                hasData = True

    return (expressionDataList, hasData)

def classifyByMutation(patients, proteinSet, mutationList):
    retval = []

    # Split the patients by whether or not they've got this mutation
    normalPatients = []
    mutatedPatients = []

    for patient in patients.valid():
        hasAllMutations = True

        for mutation in mutationList:
            if not patient.getMutations().has(mutation):
                hasAllMutations = False
                break

        if hasAllMutations:
            mutatedPatients.append(patient)
        else:
            normalPatients.append(patient)

    # Ignore very small sample sizes
    if len(normalPatients) < MIN_N or len(mutatedPatients) < MIN_N:
        say("  %s: mutated in only %d patients\n" % (mutationList, len(mutatedPatients)))
        return []

    sys.stderr.write("Testing %s (n=%d)...\n" % ("+".join(mutationList), len(mutatedPatients)))

    # Find significance of this split for each protein's expression data
    for protein in proteinSet:
        if PROTEINS_OF_INTEREST:
            if not protein in PROTEINS_OF_INTEREST:
                continue

        (normalExpressions, normalHasData) = getExpressionDataForProtein(normalPatients, protein)
        (mutatedExpressions, mutatedHasData) = getExpressionDataForProtein(mutatedPatients, protein)

        # If all expressions are 0 in both classes, ignore. (Invalid for mann-whitney u test)
        if not normalHasData and not mutatedHasData:
            continue

        # Check again for small sample sizes, in case not all patients have
        # expression data for all proteins
        if len(normalExpressions) < MIN_N or len(mutatedExpressions) < MIN_N:
            continue

        try:
            (u, p) = scipy.stats.mannwhitneyu(normalExpressions, mutatedExpressions)
        except Exception, e:
            print e
            print normalExpressions
            print mutatedExpressions
            continue

        retval.append(Result(mutationList, protein, p, normalExpressions, mutatedExpressions))

    return retval

# Generate statistics about number of patients who share mutations
def getCommonMutationsAndPrintMutationStatistics(patients, mutationSet):
    patientsPerMutation = collections.Counter()
    for patient in patients.valid():
        patientsPerMutation.update(patient.getMutations().all())

    # Get common mutations
    retval = [item[0] for item in patientsPerMutation.items() if item[1] >= MIN_N]

    # Sort first by number of mutations (descending), then by mutation name.
    # We do this by constructing a string as the key that contains the (negated) number of mutations
    # followed by the mutation name, to break ties.
    patientsPerMutationSorted = sorted(patientsPerMutation.items(),
                                       key=lambda item: "%05d%s" % (len(patients.valid()) - item[1], item[0]))

    totalP = len(patients.valid())
    cumulativeSet = set()

    say("\nMutation prevalence\n")
    cutoffReached = False

    for mutationStat in patientsPerMutationSorted:
        if mutationStat[1] <= 1:
            break

        if not cutoffReached and mutationStat[1] < MIN_N:
            say("---cutoff--- (mutations below are not considered)\n")
            cutoffReached = True
            
        for patient in patients.valid():
            if patient.getMutations().has(mutationStat[0]):
                cumulativeSet.add(patient)

        say("%-10s %2d patients, cumulative %2d/%2d (%.2f%%)\n" % (
            mutationStat[0],
            mutationStat[1],
            len(cumulativeSet),
            totalP,
            100.0 * len(cumulativeSet) / totalP
            ))

    return retval


def classifyMutationSignificance(patients, expressionProteinSet, mutationSet, csvOutputFile):
    # Test each mutation for significance
    results = []
    i = 0

    csvOutputFile.write(",".join(["Mutation","ProteinExpressed","UnmutatedMean","UnmutatedStdDev","MutatedMean","MutatedStdDev","n","p_post_bonferroni"]))
    csvOutputFile.write("\n")

    for setSize in [1, 2]:
        for mutationList in itertools.combinations(mutationSet, setSize):
            i += 1
            retval = classifyByMutation(patients, expressionProteinSet, mutationList)

            if len(retval) > 0:
                sys.stderr.write("%3d/%3d [%-10s] %d candidate proteins\n" % (
                    i, len(mutationSet),
                    "+".join(mutationList),
                    len(retval),
                ))
                results.extend(retval)

    # Apply bonferonni correction
    correctionFactor = len(results)
    for result in results:
        result.significance *= correctionFactor

    # Sort 
    results.sort(key=lambda result: result.significance)
    correlationCounter = collections.Counter()

    for result in results:
        if not PROTEINS_OF_INTEREST and result.significance > 0.01:
            break;

        say("%s\n" % (result))
        csvOutputFile.write("%s\n" % (result.csv()))

        correlationCounter[result.getMutationList()] += 1

    say("\nSignificant effects:\n")

    for mut, count in correlationCounter.most_common():
        say("%20s %d\n" % (mut, count))


def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: %s <output-filename-prefix>\n" % (sys.argv[0]))
        sys.exit(1)

    global outFile
    outFile = open(sys.argv[1] + ".txt", "w")
    csvOutputFile = open(sys.argv[1] + ".csv", "w")

    patients = Patients()

    # Generate the set of proteins for which at least one valid patient has expression data
    expressionProteinSet = set()
    for patient in patients.valid():
        expressionProteinSet.update(patient.getExpressionData().getAllProteins())

    # Generate the set of mutations in at least one valid patient
    mutationSet = set()
    for patient in patients.valid():
        mutationSet.update(patient.getMutations().all())

    say("Valid patients: %d mutations, expression data for %d proteins\n" % (
        len(mutationSet), len(expressionProteinSet)))

    commonMutations = getCommonMutationsAndPrintMutationStatistics(patients, mutationSet)

    sys.stderr.write("Common mutations: %s\n" % (commonMutations))

    classifyMutationSignificance(patients, expressionProteinSet, commonMutations, csvOutputFile)

main()


        
