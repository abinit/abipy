import seaborn as sns
sns.set_style("darkgrid")
tips = sns.load_dataset("tips")
ax = sns.pointplot(x="time", y="total_bill", data=tips)
