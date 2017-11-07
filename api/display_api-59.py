ax = sns.violinplot(x="day", y="total_bill", data=tips,
                    inner=None, color=".8")
ax = sns.stripplot(x="day", y="total_bill", data=tips, jitter=True)
