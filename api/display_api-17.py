g = sns.lmplot(x="size", y="total_bill", hue="day", col="day",
               data=tips, aspect=.4, x_jitter=.1)
