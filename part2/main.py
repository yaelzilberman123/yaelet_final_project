import os
import random
from evo_dist import run_full_evo_dist, run_evo_dist
from codon_utils import read_codon_table, dna_to_rna, sort_genes, count_codons, count_amino_acids, create_profile
from classify import run_classify

if __name__ == "__main__":
    # קריאה של טבלת הקודונים והחומצות האמיניות
    codon_table = read_codon_table("data/codon_AA.txt")

    folder_path = "data/gene"

    for file_name in os.listdir(folder_path):
        full_path = os.path.join(folder_path, file_name)

        with open(full_path, 'r') as file:

            # קריאה של הגנים מהקובץ וסידורם ברשימה
            genes = sort_genes(file)

            # ערבוב הגנים בצורה רנדומלית
            random.shuffle(genes)

            # פיצול הנתונים ל90% לאימון ו10% לבדיקת המודל
            split = int(len(genes) * 0.9)

            training_genes = genes[:split]
            test_genes = genes[split:]

            
            total_codon_counts = {codon: 0 for codon in codon_table}

            for gene in training_genes:
                # המרת רצף הדנ"א לרנ"א
                rna = dna_to_rna(gene)

                # ספירת קודונים עבור הגן הנוכחי
                gene_counts = count_codons(rna, codon_table)

                # הוספה לסכום הכולל של כל הקודונים
                for codon in total_codon_counts:
                    total_codon_counts[codon] += gene_counts[codon]
            
            amino_acid_counts = count_amino_acids(total_codon_counts, codon_table)

            profile = create_profile(total_codon_counts, amino_acid_counts, codon_table)

            with open("results/gene_results/result_" + file_name, "w") as result_file:
                
                result_file.write("---- The profile of this organism ----\n")

                # העלאת הפרופיל לקובץ
                result_file.write(str(profile.tolist()) + "\n\n")

                # העלאת הגנים שנבחרו לבדיקה לקובץ
                result_file.write("---- genes for testing ----\n")
                for gene in test_genes:
                    result_file.write(gene + "\n")
    
    run_classify(codon_table)


    evo_distances = {
    "b_subtilis_cds.fna":      3000,
    "e_coli_cds.fna":             0,
    "p_aeruginosa_cds.fna":    1500,
    "s_enterica_cds.fna":       120,
    "v_cholerae_cds.fna":       500,
    "y_enterocolitica_cds.fna": 300,
    }

    run_evo_dist("e_coli_cds.fna", codon_table, evo_distances,
             trendline_exclude=["p_aeruginosa_cds.fna"])