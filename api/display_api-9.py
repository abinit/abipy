g = sns.jointplot("petal_length", "sepal_length", data=iris,
                  marginal_kws=dict(bins=15, rug=True),
                  annot_kws=dict(stat="r"),
                  s=40, edgecolor="w", linewidth=1)
