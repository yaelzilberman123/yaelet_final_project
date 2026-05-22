import os

def sort_genes(file):
    gene = ""
    gene300 = list()
    for line in file:
        if line[0] == ">":
            gene_length = len(gene)
            if len(gene) >= 300 and len(gene) % 3 == 0:
                gene300.append(gene)
            gene = ""
            continue
        else:
            line = line.rstrip("\n\r")
            line = line.upper()
            gene = gene + line
    gene_length = len(gene)
    if len(gene) >= 300 and len(gene) % 3 == 0:
        gene300.append(gene)
    return gene300

#פונקציה המקבלת קובץ דנ"א והופכת אותו לרצף רנ"א
def DNA_RNA_Cod(DNA):
  delimiter = ""
  line_list = []
  DNA = DNA.upper()
  for ch in DNA:
    if ch == "T":
      line_list.append("U")
    else:
      line_list.append(ch)
  RNA = delimiter.join(line_list)
  return RNA     

def Read_dict(): 
  with open("data/AA_codons.txt","r") as c: # open file
    for line in c: # run on every line in th file
        ll=str(line.rstrip("\n"))
        name=ll[:2]
        amino_dict[name]=0
        tmp_codon=ll[2:]
        list_codon=tmp_codon.split(";")
        RNA_codon_table[name]= list_codon #the key is the prot and the value is the codons        
        for i in list_codon:
          codon_dict[i]=0 #the key is the codon and the value is the amount of times he is in the dna           
       
            
def codon_counter(gene):
    for i in range(0, len(gene),3):  
        codon=gene[i:i+3]
        if codon == "UAA" or "UGA" or "UAG":
            continue
        else:
            codon_dict[codon]+=1

def amino_counter():
    for amino in RNA_codon_table:
        sum_amino=0
        for codon in RNA_codon_table[amino]:
            sum_amino+=codon_dict[codon]
        
        amino_dict[amino]+= sum_amino


def create_profile():
    profile_dict=dict()
    for amino in RNA_codon_table:
        for codon in RNA_codon_table[amino]:
            pre_codon= ((codon_dict[codon]*100)/ amino_dict[amino])
            profile_dict[codon]= pre_codon
            return 






#main

if __name__ == "__main__":
    RNA_codon_table=dict() #prot-codons
    codon_dict=dict() #codon-how many
    amino_dict=dict() # amino-how many


    #open file
    
folder_path = "data/gene"

for file_name in os.listdir(folder_path):
    full_path = os.path.join(folder_path, file_name)
    with open (file_name ,"r") as file:
        with open ("result_" + file_name, "w") as out_file:
            temp_list=[]
            list_gene = sort_genes(file_path)
            count_gene=len(list_gene)
            while len(test_list) > 0.9*length:
                temp_list.append(list_gene.pop(i))
                for seq in temp_list:
                    gene= DNA_RNA_Cod(seq)
                    codon_counter(gene)
                    amino_counter()
                    file.write(create_profile())
                    list_gene.append(temp_list.pop(i))
            file.write(list_gene)