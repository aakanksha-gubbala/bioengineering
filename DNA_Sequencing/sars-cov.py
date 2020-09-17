# Severe acute respiratory syndrome coronavirus 2 isolate SARS-CoV-2/human/IND/GBRC199/2020, complete genome
# LOCUS       MT635858               29800 bp    RNA     linear   VRL 18-JUN-2020
# DEFINITION  Severe acute respiratory syndrome coronavirus 2 isolate
#             SARS-CoV-2/human/IND/GBRC199/2020, complete genome.
# ACCESSION   MT635858

from collections import Counter

sequence = open("sars-cov.txt").read().replace("\n", '')

print("\nFrequency of bases in the sequence: ", Counter(sequence))

dict = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G'}


def GC_perecent(sequence):
    gc, total = 0, 0
    for x in sequence:
        if x == 'C' or x == 'G':
            gc += 1
        total += 1
    return gc * 100 / total


print("\nPerecent GC : %0.3f" % GC_perecent(sequence))


def direct_complement(sequence):
    complement = ''
    for x in sequence:
        complement += dict[x]
    return complement


def reverse_complement(sequence):
    complement = ''
    for x in sequence:
        complement = dict[x] + complement
    return complement


pre_mRNA = reverse_complement(sequence)


def pre_mRNA_to_mRNA(mRNA):
    rna = []
    for x in range(0, len(mRNA) - len(mRNA) % 3, 3):
        rna += [mRNA[x:x + 3]]

    rna_proteins = []
    s = rna.index('AUG')
    while s > 0:
        s = rna.index('AUG')
        rna = rna[s:]
        n = min(rna.index('UAA'), rna.index('UAG'), rna.index('UGA'))
        rna_proteins += rna[:n + 1]
        rna = rna[n + 1:]
    return rna_proteins


mRNA = pre_mRNA_to_mRNA(pre_mRNA)


def RNA_to_protein(mRNA):
    RNA_codon = {'UUU': 'Phe', 'UUC': 'Phe',
                 'UUA': 'Leu', 'UUG': 'Leu', 'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
                 'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
                 'AUG': 'Met',
                 'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
                 'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser', 'AGU': 'Ser', 'AGC': 'Ser',
                 'UAU': 'Tyr', 'UAC': 'Tyr',
                 'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
                 'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
                 'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
                 'CAU': 'His', 'CAC': 'His',
                 'CAA': 'Gln', 'CAG': 'Gln',
                 'AAU': 'Asn', 'AAC': 'Asn',
                 'AAA': 'Lys', 'AAG': 'Lys',
                 'GAU': 'Asp', 'GAC': 'Asp',
                 'GAA': 'Glu', 'GAG': 'Glu',
                 'UGU': 'Cys', 'UGC': 'Cys',
                 'UGG': 'Trp',
                 'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'AGA': 'Arg', 'AGG': 'Arg',
                 'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
                 'UAA': 'STOP', 'UAG': 'STOP', 'UGA': 'STOP'}
    protein = ''
    for x in range(0, len(mRNA)):
        protein += RNA_codon[mRNA[x]] + '-'
    return protein.split('STOP')


proteins = RNA_to_protein(mRNA)

print("\nDNA :", sequence)
print("RNA :", ''.join(mRNA))
print("AA sequences (from exons) :", proteins)


def DNA_RNA_Protein(DNA):
    """ Input DNA is in 3' to 5' direction and the output is the corresponding RNA strand and amino acid sequence. """
    RNA = direct_complement(DNA)
    rna = []
    for x in range(0, len(RNA) - len(RNA) % 3, 3):
        rna += [RNA[x:x + 3]]
    protein = ''.join(RNA_to_protein(rna))
    return print("\nDNA: 5'-%s-3'" % DNA, "\nRNA: 5'-%s-3'" % RNA, "\nProtein: -%s" % protein)


DNA_RNA_Protein('GGCCCTGCT')
