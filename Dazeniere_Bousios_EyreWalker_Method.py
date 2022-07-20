#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 15:25:40 2020

@author: Adam Eyre-Walker

Programme to perform the analysis in Dazeniere et al. (2022) G3
The programe tabulates variants at synonymous and non-synonymous sites, assuming the universal genetic code, at frequencies in the alignment; note the programme refers to this as the SFS, but it is not a site frequency spectrum.
The programme expectes data in FASTA format
"""

import re
import math

aminoAcid = {'nnn':'nothing','---':'gap','atg':'methionine', 'att': 'isoleucine', 'atc': 'isoleucine', 'ata': 'isoleucine', 'ctt': 'leucine', 'ctc': 'leucine', 'cta': 'leucine', 'ctg': 'leucine', 'tta': 'leucine', 'ttg': 'leucine', 'gtt': 'valine', 'gtc': 'valine', 'gta': 'valine', 'gtg': 'valine', 'ttt': 'phenylalanine', 'ttc': 'phenylalanine', 'atg': 'methionine', 'tgt': 'cysteine', 'tgc': 'cysteine', 'gct': 'alanine', 'gcc': 'alanine', 'gca': 'alanine', 'gcg': 'alanine', 'ggt': 'glycine', 'ggc': 'glycine', 'gga': 'glycine', 'ggg': 'glycine', 'cct': 'proline', 'ccc': 'proline', 'cca': 'proline', 'ccg': 'proline', 'act': 'threonine', 'acc': 'threonine', 'aca': 'threonine', 'acg': 'threonine', 'tct': 'serine', 'tcc': 'serine', 'tca': 'serine', 'tcg': 'serine', 'agt': 'serine', 'agc': 'serine', 'tat': 'tyrosine', 'tac': 'tyrosine', 'tgg': 'tryptophan', 'caa': 'glutamine', 'cag': 'glutamine', 'aat': 'asparagine', 'aac': 'asparagine', 'cat': 'histidine', 'cac': 'histidine', 'gaa': 'glutamic acid', 'gag': 'glutamic acid', 'gat': 'aspartic acid', 'gac': 'aspartic acid', 'aaa': 'lysine', 'aag': 'lysine', 'cgt': 'arginine', 'cgc': 'arginine', 'cga': 'arginine', 'cgg': 'arginine', 'aga': 'arginine', 'agg': 'arginine', 'taa': 'stop', 'tag': 'stop', 'tga': 'stop'}

#fileName = 'Zmay_v4_chr_ALL.orfis.multilevel.ALL.list_ids.nt.fasta.gag.opie.selected_len.TranslatorX_muscle.nt_ali.fasta'
fileName = 'Zmay_v4_chr_ALL.orfis.multilevel.ALL.list_ids.fasta.gag.opie.mafft.iter5.aln'

#need to supply the filename of the alignment in this line
infile = open('/Users/bafg9/Dropbox/Adaptive TEs/alignments/'+fileName,'r')
allText = infile.read()

allText = allText.lower()

seqs = re.split('>.+',allText)

for i in range(len(seqs)):
    seqs[i] = re.sub('\s','',seqs[i])

del(seqs[0])
    
#noOfGapCodons = 0
#noOfStops = 0
#noOfFrameShifts = 0
#noOfTriAllelics = 0

stopCodons = ['taa','tga','tag']

codonSets = [['Ser/Pro/Thr/Ala',1,3,'NC',['tct','tcc','tca','tcg','cct','ccc','cca','ccg','act','acc','aca','acg','gct','gcc','gca','gcg']],
             ['Phe/Leu',1,3,'NC',['ttt','ttc','ctt','ctc']],
             ['Tyr/His',1,3,'NC',['tat','tac','cat','cac']],
             ['Lys/Glu',1,3,'NC',['aaa','aag','gaa','gag']],
             ['Cys/Arg',1,3,'NC',['tgt','tgc','cgt','cgc']],
             ['Arg/Gly',1,3,'NC',['aga','agg','gga','ggg']],
             ['Iso/Val',1,3,'C',['att','atc','gtt','gtc']],
             ['Asn/Asp',1,3,'C',['aat','aac','gat','gac']],
             ['Ser/Gly',1,3,'C',['agt','agc','ggt','ggc']]]

"""
These are not independent of the above but can be used if data is limited
             ['Phe/Ser',2,3,'NC',['ttt','ttc','tct','tcc']],
             ['Leu/Pro',2,3,'NC',['ctt','ctc','cct','ccc']],
             ['Iso/Thr',2,3,'NC',['att','atc','act','acc']],
             ['Val/Ala',2,3,'NC',['gtt','gtc','gct','gcc']],
             ['Glu/Arg',2,3,'NC',['caa','cag','cga','cgg']],
             ['Lys/Arg',2,3,'NC',['aaa','aag','aga','agg']],
             ['Glu/Gly',2,3,'NC',['gaa','gag','gga','ggg']],
"""

SFSlength = len(seqs) // 2 + 1

SFSn = []
SFSs = []
for i in range(len(codonSets)):    
    SFSn.append([0]*SFSlength)   
    SFSs.append([0]*SFSlength)
    
noOfCodons = 0

for cs in range(len(codonSets)):
    allowableCodons = codonSets[cs][4]
    for i in range(0,len(seqs[0]),3):
        codons = {}
        codonTotal = 0
        for j in range(len(seqs)):
            codon = seqs[j][i:i+3]
            if codon in allowableCodons:  
                if codon not in codons:
                    codons[codon] = 0
                codons[codon] += 1
                codonTotal += 1
        if len(codons) > 0 and codonTotal > 10:   #the site has some valid codons
            SFSn[cs][0] += 1    #the zero element of the SFS is the number of sites
            SFSs[cs][0] += 1
            if len(codons) > 1:
                nonPos = codonSets[cs][1]-1
                synPos = codonSets[cs][2]-1
                
                codonList = list(codons.keys())
                nonPosition = {}
                synPosition = {}
                for codon in codons:
                    nonNucleotide = codon[nonPos]
                    if nonNucleotide not in nonPosition:
                        nonPosition[nonNucleotide] = 0
                    nonPosition[nonNucleotide] += codons[codon]
                    synNucleotide = codon[synPos]
                    if synNucleotide not in synPosition:
                        synPosition[synNucleotide] = 0
                    synPosition[synNucleotide] += codons[codon]
                
                nonFreqs = list(nonPosition.values())
                nonFreqs.sort()
                synFreqs = list(synPosition.values())
                synFreqs.sort()
                if len(nonPosition) > 1:     #analyses tri-allelics taking the rarest allele
                    SFSn[cs][nonFreqs[0]] += 1
                if len(synPosition) > 1:
                    SFSs[cs][synFreqs[0]] += 1


#need to supply a file name for the results
outfile = open('/Users/bafg9/Dropbox/Adaptive TEs/alignments/'+fileName+'_results.txt','w')

outfile.write('Number of sequences = '+str(len(seqs))+'\n\n')

summedSFSn = [0]*int(math.log2(len(seqs)//2)+2)
summedSFSs = [0]*int(math.log2(len(seqs)//2)+2)

for cs in range(len(codonSets)):
    name = codonSets[cs][0]
    outfile.write(name+'\n')
    outfile.write('Site type\tNo of sites')
    for i in range(0,int(math.log2(len(seqs)))):
        outfile.write('\tP'+str(2**i)+'-'+str(2**(i+1)-1))
    outfile.write('\n')

#output the summarised SFS for non-synonymous
    outfile.write('Non-synon\t'+str(SFSn[cs][0]))
    summedSFSn[0] += SFSn[cs][0]
    
    subSummedSFSn = [0]*int(math.log2(len(seqs)//2)+2)
    for i in range(0,int(math.log2(len(seqs)))):
        s = 0
        for j in range(2**i,2**(i+1)):
            if j < len(SFSn[cs]):
                s += SFSn[cs][j]
            else:
                break
        outfile.write('\t'+str(s))
        summedSFSn[i+1] += s
        subSummedSFSn[i+1] += s
    outfile.write('\n')
    
#output the summarised SFS for synonymous
    outfile.write('Synon\t'+str(SFSs[cs][0]))
    summedSFSs[0] += SFSs[cs][0]
    
    subSummedSFSs = [0]*int(math.log2(len(seqs)//2)+2)
    for i in range(0,int(math.log2(len(seqs)))):
        s = 0
        for j in range(2**i,2**(i+1)):
            if j < len(SFSs[cs]):
                s += SFSs[cs][j]
            else:
                break
        outfile.write('\t'+str(s))
        summedSFSs[i+1] += s
        subSummedSFSs[i+1] += s
    outfile.write('\n')
    
    outfile.write('pN/pS\t')
    for i in range(1,int(math.log2(len(seqs)))+1):
        if subSummedSFSs[i] > 0:
            outfile.write('\t'+str(subSummedSFSn[i]/subSummedSFSs[i]))
        else:
            outfile.write('\t')
    outfile.write('\n\n')
            
    
outfile.write('Summed results counts\n')

outfile.write('Site type\tNo of sites')
for i in range(0,int(math.log2(len(seqs)))):
    outfile.write('\tP'+str(2**i)+'-'+str(2**(i+1)-1))
outfile.write('\n')

outfile.write('Non-synon\t'+str(summedSFSn[0]))
for i in range(1,len(summedSFSn)):
    outfile.write('\t'+str(summedSFSn[i]))
outfile.write('\n')

outfile.write('Synon\t'+str(summedSFSs[0]))
for i in range(1,len(summedSFSs)):
    outfile.write('\t'+str(summedSFSs[i]))
outfile.write('\n')

outfile.write('pN/pS\t')
for i in range(1,len(summedSFSs)):
    if summedSFSs[i] > 0:
        outfile.write('\t'+str(summedSFSn[i]/summedSFSs[i]))
    else:
        outfile.write('\t')
outfile.write('\n')



noFreqCats = 6
def freqCat(x,n):
    y = -int(math.log2(x/n))
    if y > noFreqCats:
        z = 1
    else:
        z = int(noFreqCats-y+1)
    
    return z
    
    

outfile.write('\n\nSummed results freq\n')

outfile.write('\tSite type\tNo of sites')
for i in range(noFreqCats,0,-1):
    outfile.write('\tP'+str(1/2**i))
outfile.write('\n')

summedSFSn = [0] * (noFreqCats + 1)
summedSFSs = [0] * (noFreqCats + 1)
for cs in range(len(codonSets)):
    n = len(seqs)
    for i in range(SFSlength):
        if i == 0:
            fc = 0
        else:
            fc = freqCat(i,n)
        try:
            summedSFSn[fc] += SFSn[cs][i]
            summedSFSs[fc] += SFSs[cs][i]
        except:
            x = 0
            


outfile.write('Non-synon') 

for i in range(noFreqCats+1):
    outfile.write('\t'+str(summedSFSn[i]))
outfile.write('\n')

outfile.write('Synon')
for i in range(noFreqCats+1):
    outfile.write('\t'+str(summedSFSs[i]))
outfile.write('\n')

outfile.write('pN/pS\t')
for i in range(1,noFreqCats+1):
    if summedSFSs[i] > 0:
        outfile.write('\t'+str(summedSFSn[i]/summedSFSs[i]))
    else:
        outfile.write('\t')
outfile.write('\n')

   
outfile.close()
