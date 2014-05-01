#!/usr/bin/python

import sys
import os
import glob
import xml.dom.minidom
import csv
import numpy
import scipy
import scipy.stats

MIN_N = 5

def patientIdFromSampleId(sampleId):
    return int(sampleId.split("-")[2])

class Patients:
    def __init__(self):
        self.map = {}

        self.getMutationData()
        self.getFunctionalData()

        # Print some statistics
        has_funcs = len([p for p in self.all() if p.getFunctionalData()])
        has_muts = len([p for p in self.all() if p.getMutations().hasData()])
        has_both = len(self.valid())

        sys.stdout.write("Got functional data for %d, mutations for %d, both for %d\n" % (
            has_funcs, has_muts, has_both))

    def get(self, patientId):
        if not patientId in self.map:
            self.map[patientId] = PatientData()

        return self.map[patientId]

    def all(self):
        return self.map.values()

    def valid(self):
        return [p for p in self.all() if p.isValid()]
        
    def size(self):
        return len(self.map.keys())


    # Populates each patient record with their expression FunctionalData 
    def getFunctionalData(self):
        for resultDir in glob.glob('rna-gene-expression/*'):
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

            self.get(patientNumber).setFunctionalData(FunctionalData(countFile))

    def getMutationData(self):
        with open('tcga-laml.csv', 'r') as mutationsCSV:
            mutations = csv.reader(mutationsCSV, delimiter=',')
            headerSkipped = False

            for mut in mutations:
                if not headerSkipped:
                    headerSkipped = True
                    continue

                gene = mut[0]
                mutType = mut[8]
                sampleId = mut[15]
                patientNumber = patientIdFromSampleId(sampleId)

                # Skip "harmless" mutations
                if mutType in ['Silent', 'Intron']:
                    continue

                self.get(patientNumber).getMutations().add(gene)

class PatientData:
    def __init__(self):
        self.mutations = Mutations()
        self.functionalData = None

    def isValid(self):
        return self.getMutations().hasData() and self.getFunctionalData()

    def getMutations(self):
        return self.mutations

    def getFunctionalData(self):
        return self.functionalData

    def setFunctionalData(self, functionalData):
        self.functionalData = functionalData

class FunctionalData:
    def __init__(self, filename):
        self.protein = {}

        for line in open(filename):
            cols = line.split()
            self.protein[cols[0]] = int(cols[1])

    def proteins(self):
        return self.protein.keys()

    def has(self, protein):
        return protein in self.protein

    def get(self, protein):
        return self.protein[protein]


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
    def __init__(self, mutation, protein, significance, normalSet, mutatedSet):
        self.mutation = mutation
        self.protein = protein
        self.significance = significance
        self.normalSet = normalSet
        self.mutatedSet = mutatedSet

    def __repr__(self):
        return "%8s: %15s %s  -->  %s, n=%d, p=%f" % (
            self.mutation,
            self.protein,
            meanAndStdDev(self.normalSet),
            meanAndStdDev(self.mutatedSet),
            len(self.mutatedSet),
            self.significance)

def getFunctionalData(patientList, protein):
    functionalDataList = []

    for patient in patientList:
        funcData = patient.getFunctionalData()

        if funcData.has(protein):
            functionalDataList.append(funcData.get(protein))

    return functionalDataList

def classifyByMutation(patients, proteinSet, mutation):
    retval = []

    # Split the patients by whether or not they've got this mutation
    normalPatients = []
    mutatedPatients = []

    for patient in patients.valid():
        if patient.getMutations().has(mutation):
            mutatedPatients.append(patient)
        else:
            normalPatients.append(patient)

    # Ignore very small sample sizes
    if len(normalPatients) < MIN_N or len(mutatedPatients) < MIN_N:
        #sys.stdout.write("  %s: mutated in only %d patients\n" % (mutation, len(mutatedPatients)))
        return []

    # Find significance of this split for each protein's functional data
    for protein in proteinSet:
        normalFunctions = getFunctionalData(normalPatients, protein)
        mutatedFunctions = getFunctionalData(mutatedPatients, protein)

        # Check again for small sample sizes, in case not all patients have
        # functional data for all proteins
        if len(normalFunctions) < MIN_N or len(mutatedFunctions) < MIN_N:
            continue

        try:
            (u, p) = scipy.stats.mannwhitneyu(normalFunctions, mutatedFunctions)
        except:
            continue

        retval.append(Result(mutation, protein, p, normalFunctions, mutatedFunctions))

    return retval

# Generate statistics about number of patients who share mutations
def printMutationStatistics(patients, mutationSet):
    patientsPerMutation = {}
    for mutation in mutationSet:
        patientsPerMutation[mutation] = len([1 for p in patients.valid() if p.getMutations().has(mutation)])
    patientsPerMutationSorted = sorted(patientsPerMutation.items(), key=lambda item: -item[1])

    totalP = len(patients.valid())
    cumulativeSet = set()

    sys.stdout.write("\nMutation prevalance\n")

    for mutationStat in patientsPerMutationSorted:
        if mutationStat[1] <= 1:
            break

        for patient in patients.valid():
            if patient.getMutations().has(mutationStat[0]):
                cumulativeSet.add(patient)

        sys.stdout.write("%-10s %2d patients, cumulative %2d/%2d (%.2f%%)\n" % (
            mutationStat[0],
            mutationStat[1],
            len(cumulativeSet),
            totalP,
            100.0 * len(cumulativeSet) / totalP
            ))
            

def classifyMutationSignificance(patients, functionalProteinSet, mutationSet):
    # Test each mutation for significance
    results = []
    i = 0
    for mutation in mutationSet:
        i += 1
        retval = classifyByMutation(patients, functionalProteinSet, mutation)

        if len(retval) > 0:
            sys.stderr.write("%3d/%3d [%-10s] %d candidate proteins\n" % (
                i, len(mutationSet),
                mutation,
                len(retval),
            ))
            results.extend(retval)

    # Apply bonferonni correction
    correctionFactor = len(results)
    for result in results:
        result.significance *= correctionFactor

    # Sort 
    results.sort(key=lambda result: result.significance)

    for result in results:
        if result.significance > 0.01:
            break;

        sys.stdout.write("%s\n" % (result))


def main():
    patients = Patients()

    # Generate the set of proteins for which at least one valid patient has functional data
    functionalProteinSet = set()
    for patient in patients.valid():
        functionalProteinSet.update(patient.getFunctionalData().proteins())

    # Generate the set of mutations in at least one valid patient
    mutationSet = set()
    for patient in patients.valid():
        mutationSet.update(patient.getMutations().all())

    sys.stdout.write("Valid patients: %d mutations, functional data for %d proteins\n" % (
        len(mutationSet), len(functionalProteinSet)))

    printMutationStatistics(patients, mutationSet)

    classifyMutationSignificance(patients, functionalProteinSet, mutationSet)

main()


        
