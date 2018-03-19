g = sns.factorplot(x="sex", y="total_bill",
                   hue="smoker", col="time",
                   data=tips, kind="swarm",
                   size=4, aspect=.7);
