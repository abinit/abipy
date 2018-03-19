ax = sns.barplot("day", "total_bill", data=tips,
                 linewidth=2.5, facecolor=(1, 1, 1, 0),
                 errcolor=".2", edgecolor=".2")
