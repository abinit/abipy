g = sns.factorplot(x="time", y="pulse", hue="kind",
                   col="diet", data=exercise,
                   size=5, aspect=.8)
