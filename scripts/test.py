#functions!!!


def Read_dict():
  global RNA_codon_table
  global codon_dict
  global amino_dict
  file = open("data/codon_AA.txt")
  for line in file:
    line = line.rstrip("\r\n")
    parts = line.split()
    codon, AA = parts[0], parts[1]
    if AA not in RNA_codon_table:
        RNA_codon_table[AA] = [codon]
    else:
        RNA_codon_table[AA].append(codon)
    codon_dict[codon] = 0
  file.close()       
#i think is ok but not sure

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
#works!!!!

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
#works!!!!

def codon_counter(gene):
    for i in range(0, len(gene),3):  
        codon=gene[i:i+3]
        if codon in ["UAA", "UGA", "UAG"]:
            continue
        else:
            codon_dict[codon]+=1
#fixed & workssss!!! thank goddddd

def amino_counter():
    for amino in RNA_codon_table:
        sum_amino=0
        for codon in RNA_codon_table[amino]:
            sum_amino += codon_dict[codon]
        amino_dict[amino] = sum_amino
    print(amino_dict)

def create_profile():
    profile_dict=dict()
    for amino in RNA_codon_table:
        for codon in RNA_codon_table[amino]:
            if amino_dict[amino] == 0:
               profile_dict[codon] = 0
            else:
                pre_codon= ((codon_dict[codon]*100)/ amino_dict[amino])
                profile_dict[codon]= pre_codon
    return profile_dict

# main program

RNA_codon_table = {}
codon_dict = {}
amino_dict = {}
Read_dict()

test_file = open("data/test_file")
gene300 = sort_genes(test_file)
for gene in gene300:
   rna_gene = DNA_RNA_Cod(gene)
   codon_counter(rna_gene)

amino_counter()
test_profile = create_profile()

print(test_profile)


