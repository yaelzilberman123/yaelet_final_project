import numpy as np


profile_ecoli = np.array([5,10,46,34])
new_profile = np.array([5,10,44,40])

dist = np.linalg.norm(profile_ecoli - new_profile)
print(dist)