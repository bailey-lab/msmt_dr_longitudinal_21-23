'''
a special purpose (and very hacky) script to calculate confidence intervals for
the inverse proportions of mdr1-86Y - takes manually created hardcoded inverted
fractions as input (e.g. 1/20 becomes 19/20) and outputs confidence intervals.
'''
from statsmodels.stats.proportion import proportion_confint

original_values=[[20, 20], [60, 60], [38, 38], [140, 141], [148, 148],
[155, 155], [69, 69], [119, 120], [138, 138], [26, 26], [96, 96], [180, 180],
[64, 64], [191, 191], [177, 179], [79, 79], [214, 215]]


for pair in original_values:
	alt_count, cov_count=pair
	lower_bound, upper_bound=proportion_confint(count=alt_count, nobs=cov_count, alpha=0.1)
	lower_bound, upper_bound=round(lower_bound*100, 1), round(upper_bound*100, 1)
	print(f'{lower_bound}%-{upper_bound}%')