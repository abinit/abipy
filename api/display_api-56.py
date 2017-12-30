g = sns.FacetGrid(tips, col="time", size=4, aspect=.7)
(g.map(sns.violinplot, "sex", "total_bill", "smoker", split=True)
  .despine(left=True)
  .add_legend(title="smoker"))  # doctest: +ELLIPSIS
# <seaborn.axisgrid.FacetGrid object at 0x...>
