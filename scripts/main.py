import os
#functions!!!

#פונקציה הפותחת מילון שהמפתחות הם החומצות האמיניות והערכים הם רשימה של הקודונים המקודדים אליהן, ופותחת גם מילונים לספירת כלל הקודונים והחומצות האמיניות
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

#פונקציה המקבלת קובץ עם גנים של אורגניזם ומחזירה רשימה של הגנים בעלי יותר מ300 נוקלאוטידים
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
#works!!!!
#פונקציה המקבלת רצף של גן, סופרת כמה פעמים מופיע כל קודון ומחזירה מילון עם הכמויות של כל קודון
def codon_counter(gene):
    for i in range(0, len(gene),3):  
        codon=gene[i:i+3]
        if codon in ["UAA", "UGA", "UAG"]:
            continue
        else:
            codon_dict[codon]+=1
#fixed & workssss!!! thank goddddd

#פונקציה המחשבת את סכום החומצות האמיניות מכל סוג בעזרת כמות הקודונים ומחזירה מילון של הכמויות של כל אחת מהחומצות האמיניות
def amino_counter():
    for amino in RNA_codon_table:
        sum_amino=0
        for codon in RNA_codon_table[amino]:
            sum_amino += codon_dict[codon]
        amino_dict[amino] = sum_amino

#פונקציה המקבלת מילון של סכום כל אחת מהחומצות האמיניות ומילון של סכום הקודונים מכל סוג, ויוצרת לפיהם פרופיל עם אחוז הקודונים יחסית למופע החומצה האמינית
def create_profile():
    profile_dict=dict()
    for amino in RNA_codon_table:
        for codon in RNA_codon_table[amino]:
            if amino_dict[amino] == 0:
               profile_dict[codon] = 0
            else:
                pre_codon= ((codon_dict[codon]*100)/ amino_dict[amino])
                profile_dict[codon]= pre_codon
    profile_list = list(profile_dict.values())
    return profile_list

# main program

RNA_codon_table = {}
codon_dict = {}
amino_dict = {}
Read_dict()



folder_path = "data/gene"
for file_name in os.listdir(folder_path):
    for codon in codon_dict: codon_dict[codon] = 0
    for amino in amino_dict: amino_dict[amino] = 0
    
    full_path = os.path.join(folder_path, file_name)
    with open (full_path ,"r") as file:
        with open ("results/result_" + file_name, "w") as out_file:
            temp_list=[]
            list_gene = sort_genes(file)
            count_gene=len(list_gene)
            while len(list_gene) > 0.9*count_gene:
                temp_list.append(list_gene.pop(0))
            for seq in list_gene:
                gene= DNA_RNA_Cod(seq)
                codon_counter(gene)
                amino_counter()
            out_file.write("----the profile of this oranism----\n")
            profile = str(create_profile())
            out_file.write(profile + "\n\n")
            out_file.write("----the genes for testing----\n")
            for i in temp_list:
                out_file.write("%s\n" %(i))