ax =  sns.stripplot("day", "total_bill", "smoker", data=tips,
                   palette="Set2", size=20, marker="D",
                   edgecolor="gray", alpha=.25)
