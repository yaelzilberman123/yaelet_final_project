import os
#פונקציה המקבלת קובץ עם גנים מאורגניזם, מעבירה כל גן לשורת סטריטנג אחת ומחזירה רשימה של כל הגנים עם מעל ל300 קודונים
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

#פונקציה שמקבלת קובץ חומצות אמיניות ואת הקיצורים שלהן ומכניסה אותן למילון גלובלי + פותחות מילונים לספירת החומצות האמיניות והקודונים
def Read_dict():
  global RNA_codon_table
  file = open("data/codon_AA.txt")
  for line in file:
    line = line.rstrip("\r\n")
    (codon,dev,AA) = line.partition("\t")
    codon_dict[codon] = 0
    amino_dict[AA] = 0
    RNA_codon_table[codon] = AA
  file.close()       
       
            
def codon_counter(gene):
    for i in range(0, len(gene),3):  
        codon=gene[i:i+3]
        if codon == "UAA" or "UGA" or "UAG":
            continue
        else:
            codon_dict[codon]+=1
    return codon_dict

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

Read_dict()
print(RNA_codon_table)

    #open file



folder_path = "data/gene"
for file_name in os.listdir(folder_path):
    full_path = os.path.join(folder_path, file_name)
    with open (full_path ,"r") as file:
        with open ("results/result_" + file_name, "w") as out_file:
            temp_list=[]
            list_gene = sort_genes(file)
            count_gene=len(list_gene)
            while len(list_gene) > 0.9*count_gene:
                temp_list.append(list_gene.pop(0))
                for seq in temp_list:
                    gene= DNA_RNA_Cod(seq)
                    codon_counter(gene)
                    amino_counter()
                    out_file.write(create_profile())
            out_file.write(list_gene)