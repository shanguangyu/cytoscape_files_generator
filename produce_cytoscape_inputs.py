#! /usr/bin/env python
# -*- coding:utf-8 -*-
# <Author: work@shanguangyu.com>

import os
import sys
import itertools
import json
import scipy.stats as sci
import shutil
import time


"""Produce files required for cytoscape."""


def parseArgs():
    usage = "Usage:\n%s differential_expression_file, molecular_name_file" % sys.argv[
        0]
    if len(sys.argv) != 3:
        print "Error found: \n\n", usage
    differential_expression_file = sys.argv[1]
    molecular_name_file = sys.argv[2]
    return differential_expression_file, molecular_name_file


def readData():
    differential_expression_dict = {}
    mrna_dict, lncrna_dict, mirna_dict = {}, {}, {}

    # produce differential expression value dict
    with open("./1.OriginInputFiles/DE.csv", 'rb') as fh:
        handle = fh.readlines()
    for line in handle[1:]:
        line_lst = line.strip().split(',')
        name = line_lst[0]
        value = [float(x) for x in line_lst[1:]]
        differential_expression_dict[name] = value

    # produce mrna, lncrna, mirna dict
    with open("./1.OriginInputFiles/NAME.csv", 'rb') as fh:
        handle = fh.readlines()
    for line in handle[1:]:
        line_lst = line.strip().split(',')
        # eg. PH_mr/lnc/hs_0009102
        molecular_id = line_lst[0]
        if molecular_id.startswith('PH_mr'):
            mirna_name = '-'.join(line_lst[1].split('-')[0:3])
            mirna_dict[mirna_name] = molecular_id
        # eg. {'hsa-miR-6751': 'PH_mr_0009102', ...}
        if molecular_id.startswith('PH_lnc'):
            lncrna_name = line_lst[1]
            lncrna_dict[lncrna_name] = molecular_id
        # eg. {'LINC00317': 'PH_lnc_0001282', ...}
        if molecular_id.startswith('PH_hs'):
            mrna_name = line_lst[1]
            entrez_gene = line_lst[2]
            mrna_dict[entrez_gene] = (mrna_name, molecular_id)
        # eg. {'28972': ('SPCS1', 'PH_hs_0027073'), ...}

    return differential_expression_dict, mrna_dict, lncrna_dict, mirna_dict


class ProduceCytoscapeInputs(object):
    """A class for produce cytoscape input files

    This class maintains two main structure:
    mirna-mrna network and lncrna-others network.
    """

    def __init__(self, four_dicts):
        self.de = four_dicts[0]
        self.mrna, self.lncrna, self.mirna = four_dicts[1:]
        self.interfolder = "./3.Intermediate"
        self.finalfolder = "./4.Result"

    def constructLocationDict(self):
        gtf_database = "./2.OriginDatabases/gencode_v21_annotation.gtf"
        loc_dict = {}
        if not os.path.exists(self.interfolder):
            os.makedirs(self.interfolder)
        lncrna_names = self.lncrna.keys()  # ['LINC00312', ...]
        mirna_names = [''.join([i.upper() for i in x.split('-')[1:3]])
                       for x in self.mirna]  # ['MIR27A', ...]
        mrna_names = [a for a, b in self.mrna.values()]  # ['ZNF277', ...]
        names = list(set(lncrna_names + mirna_names + mrna_names))
        names.pop(names.index('NA'))

        with open(gtf_database, 'r') as fh:
            handle = fh.readlines()

        for line in handle[5:]:  # pass comments
            line_lst = line.strip().split('\t')
            if line_lst[2] == "gene":
                gene_name = line_lst[-1].split(';')[4].split(' ')[
                    2].replace("\"", '')
                loc = line_lst[3:5]
                loc.append(line_lst[0])  # ['166903698', '166904835', 'chr6']
                if gene_name in names:
                    loc_dict[gene_name] = loc
        with open('%s/loc_db.json' % self.interfolder, 'w') as fp:
            json.dump(loc_dict, fp)

    def mirnaTargetFile(self):
        target_database = "./2.OriginDatabases/microRNA-target_hsa.txt"
        # [('hsa-miR-518a', '1826'), ...]
        mi_m_pairs = list(itertools.product(self.mirna.keys(),
                                            self.mrna.keys()))

        with open(target_database, 'r') as fh:
            handle = fh.readlines()
        target_pairs = [('-'.join(line.strip().split('\t')[0].split('-')[0:3]),
                         line.strip().split('\t')[1]) for line in handle]
        hits_pairs = list(set(mi_m_pairs) & set(target_pairs))

        with open('%s/mi_m_target.txt' % self.interfolder, 'w') as f:
            sys.stdout = f
            for hit in hits_pairs:
                print hit[0], hit[1]

    def mirnaCorrelationFile(self):
        with open('%s/mi_m_target.txt' % self.interfolder, 'r') as fo, open('%s/mi_m_correlation.txt' % self.interfolder, 'w') as fw:
            sys.stdout = fw
            for line in fo:
                line_lst = line.strip().split(' ')
                mirna, mrna = line_lst[:]
                mirna_id = self.mirna[mirna]
                mrna_id = self.mrna[mrna][1]
                vector_mirna = self.de[mirna_id]
                vector_mrna = self.de[mrna_id]
                spearman_p = sci.spearmanr(vector_mirna, vector_mrna)
                if (spearman_p[0] < 0) & (spearman_p[1] < 0.05):
                    print mirna + ',' + self.mrna[mrna][0] + ',' + ','.join([str(x) for x in list(spearman_p)])

    def lncrnaTargetFile(self, distance=10000):
        with open('%s/loc_db.json' % self.interfolder, 'r') as file:
            loc_dict = json.load(file)

        lnc_mi_pairs = list(itertools.product(
            self.lncrna.keys(), self.mirna.keys()))
        lnc_m_pairs = list(itertools.product(
            self.lncrna.keys(), self.mrna.keys()))
        with open("%s/lnc_other_target.txt" % self.interfolder, 'w') as file:
            sys.stdout = file
            for pair in lnc_m_pairs + lnc_mi_pairs:
                lncrna = pair[0]
                if 'hsa-miR' in pair[1]:
                    # Convert to std MIRNAXXX to make comparison with DB
                    agent = ''.join([x.upper()
                                     for x in pair[1].split('-')[1:3]])
                else:
                    agent = self.mrna[pair[1]][0]
                try:
                    lnc_loc = loc_dict[lncrna]
                    agent_loc = loc_dict[agent]
                    merge = zip(lnc_loc, agent_loc)
                    # merge format: [(u'157081667', u'101030041'), (u'157123004', u'101030153'), (u'chr3', u'chr14')]
                    # (A,B) (C,D) -> dist = abs(C-B)
                    dist = abs(int(merge[0][1]) - int(merge[1][0]))
                    if (len(set(merge[-1])) == 1) & (dist < distance):
                        if lncrna != agent:
                            print lncrna, agent, dist
                except:
                    pass

    def lncrnaCorrelationFile(self):
        m_dict = {v[0]: v[1] for k, v in self.mrna.iteritems()}
        # {'GGACT': 'PH_hs_0022883', ...}
        mi_dict = {''.join([x.upper() for x in k.split('-')[1:3]]): v for k, v in self.mirna.iteritems()}
        # {'MIR3665': 'PH_mr_0004512', ...}
        with open('%s/lnc_other_target.txt' % self.interfolder, 'r') as fo, open('%s/lnc_other_correlation.txt' % self.interfolder, 'w') as fw:
            sys.stdout = fw
            for line in fo:
                line_lst = line.strip().split(' ')
                lnc = self.lncrna[line_lst[0]]
                other_name = line_lst[1]
                if other_name.startswith("MIR"):
                    other = '-'.join(mi_dict[other_name].split('-')[0:3])
                else:
                    other = m_dict[other_name]
                vector_lnc = self.de[lnc]
                vector_other = self.de[other]
                spearman_p = sci.spearmanr(vector_lnc, vector_other)
                if (spearman_p[0] > 0) & (spearman_p[1] < 0.05):
                    print ','.join(line_lst) + ',' + ','.join([str(x) for x in list(spearman_p)])

    def mergeFiles(self):
        if not os.path.exists(self.finalfolder):
            os.makedirs(self.finalfolder)
        Files = ['/lnc_other_correlation.txt', '/mi_m_correlation.txt']
        with open('%s/combined_file.txt' % self.finalfolder, 'wb') as wfd:
            for f in Files:
                with open(self.interfolder + f, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd, 1024 * 1024 * 10)
                    # 10MB per writing chunk to avoid reading big file into
                    # memory.

    def produceResults(self):
        with open("%s/combined_file.txt" % self.finalfolder, 'r') as fh:
            handle = fh.readlines()

        # eg['SLC7A11': 'PH_hs_0048259', ...]
        mrna_dict = {v[0]: v[1] for k, v in self.mrna.iteritems()}
        names = []

        with open("%s/relations.txt" %self.finalfolder, 'w') as file:
            sys.stdout = file
            print "R1", "R2"
            for line in handle:
                line_lst = line.strip().split(',')
                print line_lst[0],line_lst[1]
                names.append(line_lst[0])
                names.append(line_lst[1])


        with open("%s/type.txt" % self.finalfolder, 'w') as file:
            sys.stdout = file
            print "key","attr"
            for name in set(names):
                value = 0
                try:
                    if name in self.lncrna.keys():
                        agent = lncrna[name]
                        attr = "lnc"
                    if name in mrna_dict.keys():
                        agent = mrna_dict[name]
                        attr = "mrna"
                    else:
                        agent = '-'.join(self.mirna(name).split('-')[0:3])
                        attr = "mirna"
                except:
                    pass
                print name, attr

    def workflowController(self):
        self.constructLocationDict()
        self.mirnaTargetFile()
        self.mirnaCorrelationFile()
        self.lncrnaTargetFile()
        self.lncrnaCorrelationFile()
        self.mergeFiles()
        self.produceResults()


def main():
    four_dicts = readData()
    rna = ProduceCytoscapeInputs(four_dicts)
    rna.workflowController()


if __name__ == '__main__':
    main()
