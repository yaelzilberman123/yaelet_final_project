RNA_codon_table=dict() #prot-codons
codon_dict=dict() #codon-how many
amino_dict=dict() # amino-how many



#פונקציה שמקבלת קובץ עם גנים של יצור, בודקת שהגנים טובים ומחזירה רשימה של כל הגנים מעל 300 תווים
def sort_genes(file):
    gene = ""
    gene300 = list()
    gene_counter = 0
    for line in file:
        if line[0] == ">":
            gene_length = len(gene)
            if len(gene) >= 300 and len(gene) % 3 == 0:
                gene_counter += 1
                gene300.append(gene)
            gene = ""
            continue
        else:
            line = line.rstrip("\n\r")
            line = line.upper()
            gene = gene + line
    gene_length = len(gene)
    if len(gene) >= 300 and len(gene) % 3 == 0:
        gene_counter += 1
        gene300.append(gene)
    return gene300, gene_counter


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
    for i in range(0, (len(gene)/3),3):  
        codon=gene[i:i+2]
        codon_dict[codon]+=1

def amino_counter():
    for amino in RNA_codon_table:
        sum_amino=0
        for codon in RNA_codon_table[amino]:
            sum_amino+=codon_dict[codon]
        
        amino_dict[amino]+= sum_amino


def create_profile():
    profil_dict=dict()
    for amino in RNA_codon_table:
        for codon in amino:
            pre_codon= ((codon_dict[codon]*100)/ amino_dict[amino])
            profil_dict[codon]= pre_codon


### main program ###
