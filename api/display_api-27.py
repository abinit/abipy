g = sns.factorplot(x="who", y="survived", col="class",
                   data=titanic, saturation=.5,
                   kind="bar", ci=None, aspect=.6)
(g.set_axis_labels("", "Survival Rate")
  .set_xticklabels(["Men", "Women", "Children"])
  .set_titles("{col_name} {col_var}")
  .set(ylim=(0, 1))
  .despine(left=True))  #doctest: +ELLIPSIS
# <seaborn.axisgrid.FacetGrid object at 0x...>
