import numpy as np

# קובץ לפוקנציון שמועילות לניתוח קודונים

# פונקציה ליצירת מילון קודונים מהקובץ
def read_codon_table(file_path):
    # יצירת המילון
    codon_table = {}

    with open(file_path, 'r') as file:

        for line in file:
            # פיצול השורה למרכיבים
            parts = line.strip().split()

            # לוודא שיש בשורה רק 2 חלקים (קודון וחומצה)
            if len(parts) == 2:
                # הוספת הקודון והחומצה למילון
                codon, amino_acid = parts
                codon_table[codon] = amino_acid
    
    # codon_table - {codon : amino_acid}
    return codon_table

# פונקציה הממירה רצף דיאנאיי לאראנאיי
def dna_to_rna(dna):
    # להפוך אותיות לגדולות ולהחליף בין האותיות
    return dna.upper().replace('T', 'U')

# פונקציה המוודאת שרצף הגן תקין
def sort_genes(file, min_length=300):
    genes = []
    current = []

    for line in file:
        line = line.strip()
        
        # לוודא שהשורה לא ריקה
        if not line:
            continue
        
        if line.startswith(">"):
            gene = "".join(current).upper()

            if len(gene) >= min_length and len(gene) % 3 == 0:
                genes.append(gene)

            current = []
        
        else:
            current.append(line)
    
    # נשאר עוד גן אחרון בקובץ שלא נדבק
    gene = "".join(current).upper()
    if len(gene) >= min_length and len(gene) % 3 == 0:
        genes.append(gene)
    
    # genes - רשימה של גנים תקינים (רצפים ארוכים מ300 שמתחלקים ב3)
    return genes


def count_codons(gene, codon_dict):
    
    # יצירת מילון ספירה עם כל הקודונים שהוגדרו מראש
    counts = {codon: 0 for codon in codon_dict}

    for i in range(0, len(gene), 3):
        codon = gene[i:i+3]

        # לדלג על קודוני עצירה
        if codon in ["UAA", "UGA", "UAG"]:
            continue
        
        if codon in codon_dict:
            counts[codon] += 1

    # counts - {codon : count}
    return counts


# codon counts - {codon : count}
# codon table - {codon : amino_acid}
def count_amino_acids(codon_counts, codon_table):
    amino_counts = {}

    for codon, amino in codon_table.items():

        if amino not in amino_counts:
            amino_counts[amino] = 0
        
        amino_counts[amino] += codon_counts.get(codon, 0)

    # amino_counts - {amino_acid : count}
    return amino_counts

# codon_counts - {codon : count}
# amino_counts - {amino_acid : count}
# codon_table - {codon : amino_acid}
def create_profile(codon_counts, amino_counts, codon_table):
    profile = []

    for codon, amino in codon_table.items():
        if amino_counts[amino] == 0:
            profile.append(0.0)
        else:
            # codon_counts[codon] - כמה פעמים הקודון הזה מופיע
            # amino_counts[amino] - כמה פעמים החומצה האמינית הזו מופיעה
            percentage = (codon_counts[codon] / amino_counts[amino]) * 100
            profile.append(percentage)
    
    return np.array(profile)

# existing_profiles - רשימה של פרופילים קיימים (כל פרופיל הוא מערך של אחוזים)
# new_profile - פרופיל חדש להשוואה (מערך של אחוזים)
def find_closest_profile(existing_profiles, new_profile):
    min_distance = float('inf')
    closest_profile = None

    for profile in existing_profiles:
        # חישוב המרחק בין הפרופיל הקיים לפרופיל החדש
        distance = np.linalg.norm(profile - new_profile)

        if distance < min_distance:
            min_distance = distance
            closest_profile = profile

    # closest_profile - הפרופיל הקיים שהכי קרוב לפרופיל החדש
    return closest_profile


    