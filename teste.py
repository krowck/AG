# Mathieu Blondel, February 2012
# License: BSD 3 clause

# Port to Python of examples in chapter 5 of
# "Introductory Statistics with R" by Peter Dalgaard

import numpy as np
from scipy.stats import ttest_1samp, wilcoxon, ttest_ind, mannwhitneyu

group1 = np.array([
0.000000,
0.000000,
0.000000,
0.000000,
0.000000,
0.000000,
0.000000,
0.000000,
0.000000,
0.000000,



])

group2 = np.array([
190.062027,
127.591934,
155.237946,
216.693985,
163.897079,
138.015213,
213.000259,
142.128189,
118.371391,
114.123024,



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