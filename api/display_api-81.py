ax = sns.pointplot(x="time", y="total_bill", hue="smoker",
                   data=tips,
                   markers=["o", "x"],
                   linestyles=["-", "--"])
