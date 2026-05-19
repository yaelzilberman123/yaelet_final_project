RNA_codon_table=dict() #prot-codons
codon_dict=dict() #codon-how many
amino_dict=dict() # amino-how many

# file_path -> gets the path of the DNA file (for example: "DNA.fa")
def slice_gen(file_path):
    fasta_file=open(file_path ,'r')
    for line in fasta_file:
        gene=""
        ll=str(line.rstrip("\n"))
        for line, i in enumerate(faste_file):
            if (line[0]='A')|(line[0]='C')|(line[0]='C')|(line[0]='T'):
                gene+=line[i]
            else: break
        codon_counter(gene)

    

def Read_dict(): 
  with open("data/AA_codons.txt","r") as c: # open file
    for line in c: # run on every line in the file
        ll=str(line.rstrip("\n"))
        name=ll[:2]
        amino_dict[name]=0
        tmp_codon=ll[2:]
        list_codon=tmp_codon.split(";")
        RNA_codon_table[name]= list_codon #the key is the prot and the value is the codons        
        for i in list_codon:
          codon_dict[i]=0 #the key is the codon and the value is the amount of times he is in the dna           
       
            
def codon_counter(gene):
    for i in range(0, (len(gene)/3),3):  #האם זה אמור להיות מינוס 3?
        codon=gene[i]+gene[i+1]+gene[i+2]
        codon_dict[codon]+=1

def amino_counter():
    for amino in RNA_codon_table:
        sum_amino=0
        for codon in RNA_codon_table[amino]:
            sum_amino+=codonn_dict[codon]
        
        ammino_dict[amino]+= sum_amino


def create_profile():
    profil_dict=dict()
    for amino in RNA_codon_table:
        for codon in amino:
            pre_codon= ((codon_dict[codon]*100)/ amino_dict[amino])
            profil_dict[codon]= pre_codon


            

            
        

        



    





