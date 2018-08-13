g = sns.lmplot(x="total_bill", y="tip", hue="smoker", data=tips,
               markers=["o", "x"])
