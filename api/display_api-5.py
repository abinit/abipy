g = (sns.jointplot("sepal_length", "sepal_width",
                   data=iris, color="k")
        .plot_joint(sns.kdeplot, zorder=0, n_levels=6))
