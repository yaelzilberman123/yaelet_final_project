import numpy as np

def calc_dist(known_profile, new_profile):
    dist = np.linalg.norm(known_profile - new_profile)
    return dist

# main program

e_coli_profile = np.array([1,2,3,4,5])
gene_profile = np.array([0,1,3,4,2])
stap_aureus_profile = np.array([2,3,4,5,6])

dist_e = calc_dist(e_coli_profile, gene_profile)
dist_stap = calc_dist(stap_aureus_profile, gene_profile)

if dist_e > dist_stap:
    print("e.coli")
elif dist_e < dist_stap:
    print("stap aureus")
else:
    print("both")