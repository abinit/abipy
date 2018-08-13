g = sns.lmplot(x="size", y="total_bill", hue="day", col="day",
               data=tips, height=6, aspect=.4, x_jitter=.1)
