g = sns.lmplot(x="total_bill", y="tip", row="sex", col="time",
               data=tips, size=3)
g = (g.set_axis_labels("Total bill (US Dollars)", "Tip")
      .set(xlim=(0, 60), ylim=(0, 12),
           xticks=[10, 30, 50], yticks=[2, 6, 10])
      .fig.subplots_adjust(wspace=.02))
