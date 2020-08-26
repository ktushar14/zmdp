# Example of results csv file:
#
# seconds,calls_to_solver,lb,ub,regret
# 3,1,10.4318,27.264,16.8322
# 4,9,15.1923,26.8201,11.6278
# 5,17,15.7434,26.2878,10.5444
# 6,26,16.8145,25.9644,9.14995
# 7,35,16.8145,25.8089,8.99444
# 8,40,16.9868,25.7557,8.76888
# 9,47,17.7256,25.7096,7.984
# 10,53,17.7256,25.6493,7.92369

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

file = "../../results/bounds_per_second.csv"
df = pd.read_csv(file)
print(df)