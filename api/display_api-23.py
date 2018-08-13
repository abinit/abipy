g = sns.lmplot(x="total_bill", y="tip", hue="smoker", data=tips,
               palette=dict(Yes="g", No="m"))
