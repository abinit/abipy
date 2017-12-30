g = sns.lmplot(x="total_bill", y="tip", col="day", hue="day",
               data=tips, col_wrap=2, size=3)
