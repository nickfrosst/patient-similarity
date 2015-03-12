#!/usr/bin/env python3

"""
out put the accuracy of patient-similarity 
"""

__author__ = 'Nick Frosst (nickfrosst@gmail.com)'

import sys
import csv
import logging
import patient_similarity
import is_same_cohort
import bisect

logger = logging.getLogger(__name__)



def script (patient_hpo_filename, hpo_filename, disease_phenotype_filename, 
           orphanet_lookup_filename=None, orphanet_prevalence_filename=None, proto=None, 
           use_disease_prevalence=False, use_phenotype_frequency=False, 
           use_patient_phenotypes=False, distribute_ic_to_leaves=False,
           use_aoo=False, scores=None):
    similarities = patient_similarity.script(patient_hpo_filename, hpo_filename, disease_phenotype_filename, orphanet_lookup_filename, orphanet_prevalence_filename, proto, use_disease_prevalence, use_phenotype_frequency, use_patient_phenotypes, distribute_ic_to_leaves, use_aoo, scores)
    similarities = is_same_cohort.script(similarities,'data/pc.manuscript.hpo')

    f1=open('./temp', 'w+')
    for sim in similarities: 
        new_sim = [sim[3]] + sim[:3]
        if new_sim[0]:
            new_sim[0] = 1
        else:
            new_sim[0] = 0
        f1.write ("\t".join(map(str,new_sim)) + "\n")

    # sort by similarity 
    similarities = sorted(similarities, key=lambda entry: entry[2],reverse=True)
    # now check accuracy
    patients = {}
    for sim in similarities: 
        if sim[0] not in patients:
            patients[sim[0]] = top_n_similar_patients(sim[0],similarities,5)

    positives =0
    for patient in patients:
        for matches in patients[patient]:
            if matches[2] == True:
                positives +=1
                break

    print(positives)
    print(len( [entry for entry in similarities if entry[3]]))


#def top_n_similar_patients (patient, sims,n):
#    top_n_sims = [0] * n
#    top_n_patients = ['a'] * n
#    top_n_is_same_co = [False] * n
#    for sim in sims:
#        if sim[0] == patient or sim[1] == patient:
#            patient2 = [x for x in sim[:2] if x not in [patient]][0] 
#            new_top_n_sims = top_n_sims
#           bisect.insort_left(new_top_n_sims,sim[2]) 
#            cut_top_n_sims = new_top_n_sims[1:]
#            if (cut_top_n_sims != top_n_sims):
#                top_n_patients.insert(new_top_n_sims.index(sim[2]),patient2)
#                top_n_is_same_co.insert(new_top_n_sims.index(sim[2]),sim[3])
#            top_n_patients = top_n_patients[1:]
#            top_n_sims = top_n_sims[1:]
#            top_n_is_same_co = top_n_is_same_co[1:]
#    return list(map((lambda x,p,c: [x,p,c]),top_n_sims,top_n_patients,top_n_is_same_co))


def top_n_similar_patients (patient, sims,n):
    top_n_sims = []
    i = 0;
    while len(top_n_sims) < n:
        patient2 = [x for x in sims[i][:2] if x not in [patient]]
        if (len(patient2) == 1):
            top_n_sims.append([patient2[0],sims[i][2],sims[i][3]])
        i+=1
    return top_n_sims



def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()
    
    parser = ArgumentParser(description=description)
    parser.add_argument('patient_hpo_filename', metavar='patients.hpo')
    parser.add_argument('hpo_filename', metavar='hp.obo')
    parser.add_argument('disease_phenotype_filename', metavar='phenotype_annotations.tab')
    parser.add_argument('--orphanet-lookup', metavar='en_product1.xml', 
                        dest='orphanet_lookup_filename', default=None)
    parser.add_argument('--orphanet-prevalence', metavar='en_product2.xml',
                        dest='orphanet_prevalence_filename', default=None)
    parser.add_argument('--use-disease-prevalence', default=False, action='store_true')
    parser.add_argument('--use-phenotype-frequency', default=False, action='store_true')
    parser.add_argument('--use-patient-phenotypes', default=False, action='store_true')
    parser.add_argument('--distribute-ic-to-leaves', default=False, action='store_true')
    parser.add_argument('--use-aoo', default=False, action='store_true')
    parser.add_argument('--log', dest='loglevel', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default='WARNING')
    parser.add_argument('--proto', metavar="file", help="HPO file of disease prototypes to compare against as well")
    parser.add_argument('-s', '--score', dest='scores', action='append', default=[],
                        choices=['jaccard', 'resnik', 'lin', 'jc', 'owlsim', 'ob', 'jz', 'ui', 'simgic', 'nsimgic', 'icca'],
                        help='Include this score in the output for each pair of patients')

    return parser.parse_args(args)


def main(args=sys.argv[1:]):
    args = parse_args(args)
    logging.basicConfig(level=args.__dict__.pop('loglevel'))
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
