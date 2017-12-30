g = sns.lmplot(x="total_bill", y="tip", row="sex", col="time",
               data=tips, size=3)
