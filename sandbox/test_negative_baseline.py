# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import hplc.quant

# %%
# Generate data
dt = 0.01
time = np.arange(0, 40, dt)
