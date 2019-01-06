import numpy as np
import pandas as pd
from math import log10
from scipy.stats import norm, pearsonr, spearmanr, kendalltau, linregress, wilcoxon
import matplotlib
matplotlib.use("cairo")
import matplotlib.pyplot as plt
plt.style.use('ggplot') # Use the ggplot style
import data_helpers

