#%%
import numpy as np
import scipy.stats
import pandas as pd

dt = 0.01
loc = [10, 10.6]
amp = [100, 20]
scale = [1, 1]
time = np.arange(0, 20, dt)
sig1 = amp[0] * scipy.stats.norm(loc[0], scale[0]).pdf(time)
sig2 = amp[1] * scipy.stats.norm(loc[1], scale[1]).pdf(time)

df = pd.DataFrame(np.array([time, sig1, sig2, sig1 + sig2]).T, columns=['time', 'signal_1', 'signal_2', 'signal'])
df.to_csv('./bounding_example.csv', index=False)