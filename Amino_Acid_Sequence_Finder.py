# -*- coding: utf-8 -*-
"""
Created on Fri May 20 07:49:34 2022

@author: Mitch
"""
import csv 
import matplotlib.pyplot as plt
import numpy as np
import sys
import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis

name = list()
sequence = list()
with open('Pacific_Oyster_Proteome_NameSequence.csv', newline= '') as csvfile:
    bigfile = csv.reader(csvfile, delimiter=',')
    for row in bigfile:
      name.append(row[0])  
      sequence.append(row[1])
      
##chunking a list by n parts

def chunk_list_by_n (listx, n):
    newlist = list()
    for i in range(0, len(listx), n):
        newlist.append(listx[i:i+n])
    return (newlist)

##counting amount of amino acids in a given sequence 
def frequency_of_aminoaids (listx):
    listofaminoacids = ('A','R','D','N','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V')
    aafreq = list()
    for aa in listx:
        for saa in listofaminoacids:
            aafreq.append(aa.count(saa))
    return(aafreq)

def MW_of_protein (selectedsequences):
    aaweights = {'A': 89.1, 'C': 121.2, 'D' : 133.1, 'E': 147.1, 'F': 165.2,
           'G': 75.12, 'H': 155.2, 'I': 131.2, 'K': 146.2, 'L': 131.2,
           'M': 149.2, 'N': 132.1, 'P': 115.1, 'Q': 146.2, 'R': 174.2,
           'S': 105.1, 'T': 119.1, 'V': 117.1, 'W': 204.2, 'Y': 181.2, 
           'U':0.0, 'X': 0.0}
    mass = ''
    masslist = list()
    for seq in selectedsequences:
        for letter in seq:
            mass = sum(aaweights[a] for a in seq)
        masslist.append(mass)
    return (masslist)


def make_dictionary (listx):
    numberkey = enumerate(listx)
    dkey = dict((key,value) for key,value in numberkey)
    return(dkey)


def find_motif (motif, dkey, dvalue,):
    
    requestedsequence = str(motif)
    pronamenumbers = list()
    for pronamenumber, sequence in dvalue.items(): # iterating with varibles pronamenumber(key) and sequence(value) in the dvalues dictionar with sequences
        if re.search(requestedsequence, sequence): #if the search for "requestedsequence" is in sequence
            pronamenumbers.append(pronamenumber)  #append the number of that sequence to pronamenumbers 
    return(pronamenumbers)


def call_selected_keyandvalue (pronamenumbers, dkey, dvalue):
    
    selectedpronames = list()
    selectedsequences = list()
    for i in pronamenumbers:
        selectedpronames.append(dkey[i])
        selectedsequences.append(dvalue[i])
    return (selectedpronames, selectedsequences)


def length_of_sequence (selectedsequences):
    
    lengthofseq = list()
    for seq in selectedsequences:
        lengthofseq.append(int(len(seq)))
    return(lengthofseq)


def divide_array (aakeydict, seqlenkey):
    
    percentaalist = list()
    for number in range(0, len(aakeydict)):
        n = np.array(aakeydict[number])
        m = np.array(seqlenkey[number])
        di = np.divide(n,m)*100
        percentaalist.append(di)
    return (percentaalist)


def assign_amino_acid_to_frequency (listx):
    listofaminoacids = ('A','R','D','N','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V')

    arraylist = list()
    for freq_array in listx:
        for val in freq_array:
            arraylist.append(val)
    arraychunck = chunk_list_by_n(arraylist, 20)
    for lista in arraychunck:
        return(dict(zip(listofaminoacids, lista)))
        
#print(pronamenumbers)
#################################################################################### 
##Work Book


##creating a dictionary for the FASTA title (key) and sequences (value)
dkey = make_dictionary(name)
dvalue = make_dictionary(sequence)


##selecting sequences wit the respective motif and assinging the titles to pronames and values to sequences 
pronamenumbers = find_motif('AAAAAAAA', dkey, dvalue) 
selectedkeyvalues = call_selected_keyandvalue(pronamenumbers, dkey, dvalue)
selectedpronames = (selectedkeyvalues[0])
selectedsequences = (selectedkeyvalues[1])

##Length of selected amino acids and their respective amino acid frequencies
lengthofsequences = length_of_sequence(selectedsequences)
aacount = frequency_of_aminoaids(selectedsequences)

##Chunking aacouts by 20 (there are 20 major amino acids) and length of sequences by 1 (1 total sequence for protein)
blockedaacount = chunk_list_by_n(aacount, 20)
blockedseq = chunk_list_by_n(lengthofsequences, 1)

##making blocked sections into a dictionary 
aakeydict = make_dictionary(blockedaacount)
seqlenkey = make_dictionary(blockedseq)

##getting percent frequency of aa 
percentaalist = divide_array(aakeydict, seqlenkey)
#print(percentaalist)

##getting mass of proteins 
masslist = MW_of_protein(selectedsequences)


## associating frequency with amino acids 
listofaminoacids1 = ('A','R','D','N','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V')

##putting amino acid count into a library with the letter as the key and the frequency percent (mole percent) as the value

arraylist = list()
aafreqlist = list()
for freq_array in percentaalist:
    for val in freq_array:
        arraylist.append(val)
arraychunck = chunk_list_by_n(arraylist, 20)
for lista in arraychunck:
    aafreqlist. append(dict(zip(listofaminoacids1, lista)))

##printing the protein name, sequence, protein mass and amino acid composition. With number of hits
zipfile = zip(selectedpronames, selectedsequences, masslist,  aafreqlist)
listzipfile = list()
for j in zipfile:
    listzipfile.append(j)
print(listzipfile)
print( "\nNumber of Proteins Found:", len(selectedpronames))
        
##Create a folder where each element in the list is in a coloumn
#with open ('Oyster_Sequence_.SKG.csv', 'w', newline='')as f:
    #wr = csv.writer(f)
    #for row in listzipfile:
        #wr.writerow(row)

