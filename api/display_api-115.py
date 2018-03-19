g = sns.factorplot(x="class", hue="who", col="survived",
                   data=titanic, kind="count",
                   size=4, aspect=.7);
