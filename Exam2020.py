#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 17:03:30 2020

@author: mario
"""

#%%
#Exercise 1a

def iruns(seq): 
    list_iruns = []
    #Initialize first irun with first base
    irun = seq[0]
    #Iterate over sequence starting by base 2
    for i in range(1, len(seq)):
        #If base is equal or "bigger" than previous one, add it to the current irun
        if seq[i] >= seq[i-1]:
            irun += seq[i]
        #Else, the irun ios over, append and start a new one with current value as first value of the new one
        else:
            list_iruns.append(irun)
            irun = seq[i]
    
    #Append last irun (it wont otherwise)
    list_iruns.append(irun)
        
    #Return iruns list
    return list_iruns

print(iruns("AAACCGGTTTAAAATTTT")) # Returns ["AAACCGGTTT","AAAATTTT"]

print(iruns("ACGTACGTACGT")) # Returns ["ACGT","ACGT","ACGT"]

print(iruns("ACACTTTTTT"))

#%%
#Exercise 1b

def riruns(seq):
    if len(seq) == 1:
        return seq[0]
    
    return seq[0] + riruns(seq[1:])

print(riruns("AAACCGGTTTAAAATTTT"))
    
    
#%%
#Exercise 1c
def dastq_fasta_iruns(fastq):
    out = open('iruns.fasta', 'wt')
    
    with open(fastq, 'rt') as f:
        while True:
            des = f.readline().rstrip() #Header line
            seq = f.readline().rstrip() #Sequence line            
            f.readline.rstrip() # + line, ignore
            f.readline.rstrip() #Quality line, ignore 
            
            list_iruns = iruns(seq)
            
            for i in range(len(list_iruns)):
                out.write(des + '.' + str((i+1)) + '\n')
                out.write(list_iruns[i])
    out.close()
    
#%%
#Exercise 1d


#%%
#Exercise 1e - Works for not too long l and n, but there is a more efficient way  for sure, this is way too complicated and unefficient
def iruns_gen(n,l, seq=''):   
    def is_irun(seq):
        for i in range(1, len(seq)):
            if seq[i] < seq[i-1]:
                return False
        return True
     
    kmer_list = []
    def kmers(k, kmer=''):
        bases = 'ACGT'
        if k == 0:
                return kmer
        else:
            for b in bases:
                kmer_list.append(kmers(k-1, b+kmer))
        return 'XAX'
    
    kmers(l)
    new_kmer_list = [kmer for kmer in kmer_list if is_irun(kmer)]
    
    def gen_iruns(k, list_iruns, seq=''):
        if k == 0:
            if len(iruns(seq)) == n:
                print(seq)
        else: 
            for kmer in list_iruns:
                gen_iruns(k-1, list_iruns, seq=kmer+seq)
                
    return gen_iruns(n, new_kmer_list)
        


print(iruns_gen(4,4))

#%%
#Exercise 1e v2 - Checking each seq with first function to confirm they are made of n  iruns
def iruns_gen2(n,l, seq=''):
    
    def is_irun(seq):
        for i in range(1, len(seq)):
            if seq[i] < seq[i-1]:
                return False
        return True
     
    kmer_list = []
    def kmers(k, kmer=''):
        bases = 'ACGT'
        if k == 0:
                return kmer
        else:
            for b in bases:
                kmer_list.append(kmers(k-1, b+kmer))
        return 'XAX'
    
    kmers(l)
    new_kmer_list = [kmer for kmer in kmer_list if is_irun(kmer)]
    
    def gen_iruns(k, list_iruns, seq=''):
        if k == 0:
            if len(iruns(seq)) == n:
                print(iruns(seq))
        else: 
            for kmer in list_iruns:
                gen_iruns(k-1, list_iruns, seq=kmer+seq)
                
    return gen_iruns(n, new_kmer_list)
 

iruns_gen2(4,4)
 

#%%
#Exercise 1e



   
#%%
#Extra: function that checks if sequence is a irun
def is_irun(seq):
        for i in range(1, len(seq)):
            if seq[i] < seq[i-1]:
                return False
        return True
    
print(is_irun('AAAAAAAAAATTCCT')) #Returns False
print(is_irun('AAAAAAAAAACCCGTTTT')) #Returns True



    
    
    

    
    
    
    
    
