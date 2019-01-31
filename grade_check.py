import matplotlib.pyplot as plt
import numpy as np

score = np.array([2,3,2,3,1,2,3,2,3,1,3,2,3,3,3,3,3])
credit = np.array([5,8,5,5,8,5,8,8,5,5,8,5,5,5,5,5,30])

grade = np.sum(score*credit)/np.sum(credit)

print(grade)