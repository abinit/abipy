g = sns.factorplot(x="age", y="embark_town",
                   hue="sex", row="class",
                   data=titanic[titanic.embark_town.notnull()],
                   orient="h", size=2, aspect=3.5, palette="Set3",
                   kind="violin", split=True, cut=0, bw=.2)
