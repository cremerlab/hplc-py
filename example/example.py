# %%
from hplc.quant import Chromatogram
from hplc.io import load_chromatogram
example = load_chromatogram('example.csv', cols=['time', 'signal'])
chrom = Chromatogram(example)
chrom.show()
peaks = chrom.fit_peaks()
# %%
chrom.show(time_range=[10, 20])
plt.savefig('./reconstructed_chromatogram.png')
