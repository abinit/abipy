x, y = np.random.randn(2, 300)
g = (sns.jointplot(x, y, kind="hex", stat_func=None)
        .set_axis_labels("x", "y"))
