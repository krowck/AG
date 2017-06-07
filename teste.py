# Mathieu Blondel, February 2012
# License: BSD 3 clause

# Port to Python of examples in chapter 5 of
# "Introductory Statistics with R" by Peter Dalgaard

import numpy as np
from scipy.stats import ttest_1samp, wilcoxon, ttest_ind, mannwhitneyu

group1 = np.array([1.223024,
1.270963,
2.420152,
1.783592,
1.868370,
1.843076,
2.267265,
1.454344,
2.045524,
1.372313
])

group2 = np.array([0.298947,
0.343682,
0.303039,
0.310800,
0.230425,
0.368618,
0.485994,
0.499390,
0.321040,
0.464811,

])

# two-sample t-test
# null hypothesis: the two groups have the same mean
# this test assumes the two groups have the same variance...
# (can be checked with tests for equal variance)
# independent groups: e.g., how boys and girls fare at an exam
# dependent groups: e.g., how the same class fare at 2 different exams
t_statistic, p_value = ttest_ind(group1, group2)

# p_value < 0.05 => alternative hypothesis:
# they don't have the same mean at the 5% significance level
print "two-sample t-test", p_value

# two-sample wilcoxon test
# a.k.a Mann Whitney U
u, p_value = mannwhitneyu(group1, group2)
print "two-sample wilcoxon-test", p_value

t, p_value = wilcoxon(group1, group2)
print p_value

# paired t-test: doing two measurments on the same experimental unit
# e.g., before and after a treatment
t_statistic, p_value = ttest_1samp(group2 - group1, 0)

# p < 0.05 => alternative hypothesis:
# the difference in mean is not equal to 0
print "paired t-test", p_value

# alternative to paired t-test when data has an ordinary scale or when not
# normally distributed
z_statistic, p_value = wilcoxon(group2 - group1)

print "paired wilcoxon-test", p_value