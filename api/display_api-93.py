from numpy import median
ax = sns.barplot(x="day", y="tip", data=tips, estimator=median)
