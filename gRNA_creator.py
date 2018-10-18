#Python script for automated gRNA cretaion

#Asking the author to input the mRNA sequence
print('Insert targeting sequence:')    
input1=input()

input1= input1.upper()

print('Length of Sequence:',len(input1))

#Predefined 5' and 3' extra sequences
hairpin_domain = 'GUUGUGGAAGGUCCAGUUUUGGGGGCUAUUACAACA'
u6_terminator = 'UUUUUU'
grna_3prime = hairpin_domain + u6_terminator
flanking_3prime = 'GTTTCGGTCTTCGCGGCG'
flanking_5prime = 'TATTATGAAGACTACACCG'


before_seq ='TATTATGAAGACTACACCG'
after_seq = 'GTTGTGGAAGGTCCAGTTTTGGGGGCTATTACAACATTTTTTGTTTCGGTCTTCGCGGCG'
after_seq2 = 'GTTGTGGAAGGTCCAGTTTTGGGGGCTATTACAACA'


kozak_seq = 'GCCACCATGG'

mut_site1 = kozak_seq.find('T')

def convert_dna_to_mrna(seq):
    new_seq= seq
    for i in range(len(seq)):
        if seq[i] == 'T':
            new_seq = new_seq[:i] + 'U' + new_seq[i+1:]
    return new_seq

def convert_mrna_to_dna(seq):
    new_seq= seq
    for i in range(len(seq)):
        if seq[i] == 'U':
            new_seq = new_seq[:i] + 'T' + new_seq[i+1:]
    return new_seq

def complementary_seq_DNA(seq):
    new_seq = seq
    for i in range(len(seq)):
        if seq[i] == 'A':
            new_seq = new_seq[:i] + 'T' + new_seq[i+1:]
        elif seq[i] == 'T':
            new_seq = new_seq[:i] + 'A' + new_seq[i+1:]
        elif seq[i] == 'C':
            new_seq = new_seq[:i] + 'G' + new_seq[i+1:]
        elif seq[i] == 'G':
            new_seq = new_seq[:i] + 'C' + new_seq[i+1:]
    return new_seq
    
def reverse_DNA(seq):
    new_seq = seq[::-1]
    return new_seq

def reverse_comp_DNA(seq):
    return reverse_DNA(complementary_seq_DNA(seq))  
    
def return_mutation_site(seq,check_seq,mut_site): #(mut_site ==7)
    index_pos = seq.find(check_seq)
    if index_pos >= 0:
        print('Recognised mutation site is at: ',index_pos+1)
        return index_pos + mut_site
    else:
        print('The sequence does not contains any of the required site.')

def U_to_C_mutation(seq,index):
    new_seq = seq
    print('New_seq = ',new_seq)
    print('The index for mutation: ',index)
    if seq[index] == 'T':
        new_seq = new_seq[:index] + 'T' + new_seq[index+1:]
    return new_seq

index1 = return_mutation_site(input1, kozak_seq, mut_site1)
seq_updated = U_to_C_mutation(input1,index1)
seq_updated = convert_dna_to_mrna(reverse_comp_DNA(seq_updated))
    
print('spacer gRNA: ',seq_updated)
print('Spacer cDNA of gRNA: ',convert_mrna_to_dna(seq_updated))

print('Simply gRNA: ', convert_dna_to_mrna(convert_mrna_to_dna(seq_updated) + after_seq2))
print('Final gRNA: ', before_seq + convert_mrna_to_dna(seq_updated) + after_seq )
