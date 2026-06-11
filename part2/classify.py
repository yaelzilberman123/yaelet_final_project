import os
import random
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from codon_utils import dna_to_rna, count_codons, count_amino_acids, create_profile, find_closest_profile, sort_genes


# פונקציה הטוענת את הפרופילים הקיימים מתיקיית התוצאות
# מחזירה רשימה של פרופילים ורשימה של שמות הקבצים המתאימים
def load_profiles(results_folder):
    existing_profiles = []
    names = []

    for file_name in sorted(os.listdir(results_folder)):
        full_path = os.path.join(results_folder, file_name)

        with open(full_path, 'r') as file:
            lines = file.readlines()

        # שורה 2 היא הפרופיל - מסירים את הסוגריים המרובעים ומפצלים לפי פסיק
        line = lines[1].strip()
        line = line[1:-1]  # הסרת הסוגריים המרובעים

        # המרת המחרוזת לרשימה של מספרים בצורת מחרוזת
        numbers = line.split(', ')

        # ממירים כל מחרוזת למספר עשרוני ויוצרים מערך נומפי
        profile = np.array([float(n) for n in numbers])

        existing_profiles.append(profile)
        names.append(file_name)

    # מחזירים שתי רשימות - פרופילים ושמות - באותו סדר
    return existing_profiles, names




# פונקציה המוצאת את גני הבדיקה בתוך קובץ התוצאות
# גני הבדיקה נמצאים אחרי הכותרת "genes for testing"
def get_test_genes(lines):
    start_index = 0
    for i, line in enumerate(lines):
        if "genes for testing" in line:
            # מתחילים מהשורה שאחרי הכותרת
            start_index = i + 1
            break

    # מסננים שורות ריקות ומחזירים רשימה של גנים תקינים
    return [line.strip() for line in lines[start_index:] if line.strip()]


# פונקציה המסווגת רשימה של גנים ומחזירה את אחוז הדיוק
# correct_profile - הפרופיל הנכון של האורגניזם שממנו לקוחים הגנים
def classify_genes(genes, existing_profiles, correct_profile, codon_table):
    true_predictions = 0

    for gene in genes:
        # בניית פרופיל עבור הגן הנוכחי
        rna = dna_to_rna(gene)
        codon_counts = count_codons(rna, codon_table)
        amino_acid_counts = count_amino_acids(codon_counts, codon_table)
        profile = create_profile(codon_counts, amino_acid_counts, codon_table)

        # מציאת הפרופיל הקרוב ביותר מבין כל הפרופילים הקיימים
        closest = find_closest_profile(existing_profiles, profile)

        # בדיקה אם הפרופיל שנמצא תואם לאורגניזם הנכון
        # np.array_equal משווה שני מערכים איבר איבר
        if np.array_equal(closest, correct_profile):
            true_predictions += 1

    # מחזירים אחוז מתוך סך כל הגנים שנבדקו
    return (true_predictions / len(genes)) * 100


# פונקציה היוצרת מפת חום המשווה את פרופילי הקודונים בין האורגניזמים
# כל שורה היא אורגניזם, כל עמודה היא קודון, הצבע מייצג את האחוז
def generate_heatmap(existing_profiles, names, codon_table, heatmap_path):
    # שמות הקודונים לציר X
    codon_labels = list(codon_table.keys())

    # שמות האורגניזמים לציר Y - מסירים את הקידומת והסיומת משם הקובץ
    organism_labels = [name.replace("result_", "").replace("_cds.fna", "") for name in names]

    # יצירת מטריצה - כל שורה היא פרופיל של אורגניזם אחד
    matrix = np.array(existing_profiles)

    # יצירת תיקיית הפלט אם לא קיימת
    os.makedirs(os.path.dirname(heatmap_path), exist_ok=True)

    # יצירת מפת החום
    # figsize - גודל התמונה במסך (רוחב, גובה באינצ')
    plt.figure(figsize=(20, 4))
    sns.heatmap(
        matrix,
        xticklabels=codon_labels,
        yticklabels=organism_labels,
        cmap="YlOrRd",      # סקאלת צבעים: צהוב = נמוך, אדום = גבוה
        linewidths=0.3,     # קווים דקים בין התאים לקריאות טובה יותר
        annot=False         # False = לא מציג מספרים בתוך התאים (עמוס מדי)
    )

    plt.title("Codon Usage Profile per Organism")
    plt.xlabel("Codon")
    plt.ylabel("Organism")

    # tight_layout מונע חיתוך של תוויות בשוליים
    plt.tight_layout()

    # שמירת התמונה לקובץ ואז הצגה על המסך
    plt.savefig(heatmap_path, dpi=150)
    plt.show()


# פונקציה ראשית המריצה את כל תהליך הסיווג
# ברירות המחדל של הנתיבים מתאימות למבנה התיקיות של הפרויקט
def run_classify(
    codon_table,
    results_folder="results/gene_results",
    final_result_path="results/final_result.txt",
    heatmap_path="results/visualizations/codon_heatmap.png"
):
    # טעינת כל הפרופילים הקיימים מקבצי התוצאות
    existing_profiles, names = load_profiles(results_folder)

    with open(final_result_path, 'w') as out_file:

        # לולאה על כל קובץ תוצאות - i הוא האינדקס של האורגניזם
        for i, file_name in enumerate(sorted(os.listdir(results_folder))):
            full_path = os.path.join(results_folder, file_name)

            with open(full_path, 'r') as file:
                lines = file.readlines()

            # חילוץ גני הבדיקה מהקובץ
            genes = get_test_genes(lines)

            # סיווג הגנים וחישוב אחוז הדיוק
            # existing_profiles[i] הוא הפרופיל הנכון עבור האורגניזם הנוכחי
            accuracy = classify_genes(genes, existing_profiles, existing_profiles[i], codon_table)

            out_file.write(f"{names[i]}: {accuracy:.2f}%\n")
            print(f"{names[i]}: {accuracy:.2f}%")

    # יצירת מפת החום לאחר סיום הסיווג
    generate_heatmap(existing_profiles, names, codon_table, heatmap_path)