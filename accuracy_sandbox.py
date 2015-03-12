#!/usr/bin/env python3

"""
"""
__author__ = 'Nick Frosst (nickfrosst@gmail.com)'

import sys
import csv
import logging
import patient_similarity
import is_same_cohort
import bisect
from random import choice, uniform 
from numpy import array, dot, random
from sklearn import tree, svm


from sklearn.externals.six import StringIO  

logger = logging.getLogger(__name__)



def script (manuscript_file_name):

    entries = []
    with open(manuscript_file_name) as ifp:
        header = ifp.readline().strip().split('\t')
        for line in ifp:
            entries += [line[:-1].split('\t')]

    entries = is_same_cohort.direct_script(entries,'data/cohorts.txt')
    
    entries = normalize(entries)

    #print(sum([entry[-1] for entry in entries])/2)
    
    #train_svm(entries)
    print (header)
    #w=perceptron(entries)
    #print("perceptron")
    #print(perceptron_accuracy(entries,0,5))
    print("simgic")
    print(index_accuracy(entries,header.index('simgic'),5))

    print("icca")
    print(index_accuracy(entries,header.index('icca'),5))
    #for i in range(len(entries[0]) - 3):
    #    index_accuracy(entries,i,5)


def similarities_accuracy(similarities,n):

    similarities = sorted(similarities, key=lambda entry: entry[2],reverse=True)
    # now check accuracy
    patients = {}
    potential = 0;
    potential_matches = 0
    for sim in similarities: 
        if sim[0] not in patients:

            patients[sim[0]],number_of_matches = top_n_similar_patients(sim[0],similarities,n) #top_n_similar_patients_get_number(sim[0],similarities,n)  = 168
            if (number_of_matches > 0):
                potential_matches += 1

    positives = 0
    for patient in patients:
        for matches in patients[patient]:
            if matches[2] == True:
                positives +=1
                break

    return (positives,potential_matches)




def index_accuracy(entries,i,n):
    similarities = []
    for entry in entries:
        similarities += [entry[:2] + [entry[i],entry[-1]]]

    return similarities_accuracy(similarities,n)


def testTrainSplit(entries,p):

    trainX = []
    trainY = []
    testX = []
    testY = []

    for entry in entries:
        if (uniform(0,1) > p):
            testX += [entry[2:-1]]
            testY += [entry[-1]]
        else:
            trainX += [entry[2:-1]]
            trainY += [entry[-1]]

    return [trainX,trainY,testX,testY]

def train_tree(entries):

    [trainX,trainY,testX,testY] = testTrainSplit(entries,0.7)

    clf = tree.DecisionTreeClassifier()
    clf = clf.fit(trainX, trainY)

    print(sum(map (lambda x,y: x !=y & y == 1,clf.predict(testX),testY)))
    print(sum(testY))

    #with open("sims.dot", 'w') as f:
    #    f = tree.export_graphviz(clf, out_file=f)

def train_svm(entries):

    [trainX,trainY,testX,testY] = testTrainSplit(entries,0.7)

    clf = svm.SVC()
    clf.fit(trainX, trainY)

    print(sum(map (lambda x,y: x ==y & y == 1,clf.predict(testX),testY)))
    print(sum(testY))

def perceptron_accuracy (entries,weights,n):
    #weights = [-0.43273271,-1.04473798,-0.33381223,-0.12731128,1.90711059,-0.78914891,-3.41655694,0.70454378,-0.96846391,2.79248341,0.53295195,-0.15101959,-0.56863352,0.70496117,3.15530632,1.51388978,0.49137292,1.14738013,-1.03259826] #perceptron learned
    #weights =  [-0.98548433,-1.13133962,-0.15100226,-0.36933738,1.86683414,-0.24893899,-3.50280371,0.20754024,-4.79387066,6.76799445,0.42069015,0.62662121,-0.67272341,0.02340744,3.30635255,0.7218613,0.9191414,0.10088303,-0.66337869]
    #weights = [0,0,0,0,1.90711059,0,0,0,0,2.79248341,0,0,0,0,3.15530632,0,0,0,0] #perceptron learned and prunned
    #weights = [108/317,83/317,102/317,105/317,106/317,78/317,98/317,121/317,116/317,116/317,117/317,110/317,95/317,82/317,100/317,120/317,115/317,120/317,121/317] #post-hoc accuracy
    #weights = [108/317,0,102/317,105/317,106/317,0,0,121/317,116/317,116/317,117/317,110/317,0,0,100/317,120/317,115/317,120/317,121/317]#post-hoc accuracy prunned

    similarities = []
    for entry in entries:
        similarities += [entry[:2] + [dot(entry[2:-1],weights),entry[-1]]]

    return similarities_accuracy(similarities,n)    

def perceptron (entries):


    w = random.rand(len(entries[0])-3)
    errors = [] 
    eta = 0.2 
    n = 10000


    unit_step = lambda x: 0 if x < 0 else 1
    percept = lambda i: unit_step(dot(w, x))

    for i in range(1000):
        for j in range(n):
            entry = choice(entries) 
            x = array(entry[2:-1])
            expected = entry[-1]
            error = expected - percept(x)
            errors.append(error) 
            w += eta * error * x
        print (n*i)
        accuracy = perceptron_accuracy(entries,w,5)
        print (accuracy)
        if (accuracy > 200):
            print (w)

    print(w)
    return(w)

def normalize (entries):
    m=entries[0]
    for i in range (2, len (entries[0]) - 1):
        m[i ] = max([el[i] for el in entries])

    normalized_entries = []
    for entry in entries:
        normalized_entries += [entry[:2]+ list(map(lambda a,b: a/b, entry[2:-1], m[2:-1])) + [entry[-1]]] 
    return normalized_entries

def top_n_similar_patients_get_number (patient, sims,n):
    top_n_sims = []
    number_of_matches = 0
    
    for sim in sims:
        patient2 = [x for x in sim[:2] if x != patient]
        if (len(patient2) == 1):
            if len(top_n_sims) < n:
                top_n_sims.append([patient2[0],[2],sim[3]])
            if sim[-1]:
                number_of_matches += 1


    return [top_n_sims, number_of_matches]

def top_n_similar_patients (patient, sims,n):
    top_n_sims = []
    number_of_matches = 0
    
    for sim in sims:
        if len(top_n_sims) >= n:
            break
        
        patient2 = [x for x in sim[:2] if x != patient]
        
        if (len(patient2) == 1):
            top_n_sims.append([patient2[0],[2],sim[3]])


    return [top_n_sims,0]



def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()
    
    parser = ArgumentParser(description=description)
    parser.add_argument('manuscript_file_name', metavar='manuscript.file')

    return parser.parse_args(args)


def main(args=sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
