import numpy as np
import matplotlib.pyplot as plt

N = 4
men_means = (0.14246744020, 0.14229883960, 0.7202706140, 0.8593171549)
men_std = (0.02018608760, 0.01473439110, 0.02521205370, 0.030318799)

ind = np.arange(N)  # the x locations for the groups
width = 0.35       # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind, men_means, width, color='r', yerr=men_std)

women_means = (0.18048029330, 0.19328810370, 0.82694409230, 0.8852732535)
women_std = (0.0146661147, 0.0105209624, 0.0102292694, 0.0046109057)
rects2 = ax.bar(ind + width, women_means, width, color='y', yerr=women_std)

# add some text for labels, title and axes ticks
ax.set_ylabel('Correlation')
ax.set_title('Prediction power by model')
ax.set_xticks(ind + width / 2)
ax.set_xticklabels(('M1', 'M2', 'M3', 'M4'))

ax.legend((rects1[0], rects2[0]), ('Before', 'After'))
plt.show()
