ax = sns.boxplot(x="day", y="total_bill", hue="smoker",
                 data=tips, palette="Set3")
