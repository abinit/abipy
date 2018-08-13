g = sns.catplot(x="sex", y="total_bill",
                hue="smoker", col="time",
                data=tips, kind="violin", split=True,
                height=4, aspect=.7);
