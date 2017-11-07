from numpy import median
ax = sns.pointplot(x="day", y="tip", data=tips, estimator=median)
