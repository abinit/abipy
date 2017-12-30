ax = sns.stripplot(x="day", y="total_bill", hue="smoker",
                   data=tips, jitter=True,
                   palette="Set2", split=True)
