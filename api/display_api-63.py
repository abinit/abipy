ax = sns.stripplot(x="sex", y="total_bill", hue="day",
                   data=tips, jitter=True)
