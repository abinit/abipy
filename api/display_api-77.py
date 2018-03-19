ax = sns.swarmplot(x="day", y="total_bill", hue="smoker",
                   data=tips, palette="Set2", dodge=True)
