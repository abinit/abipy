import seaborn as sns
sns.set(style="ticks")
exercise = sns.load_dataset("exercise")
g = sns.factorplot(x="time", y="pulse", hue="kind", data=exercise)
