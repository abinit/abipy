from scipy.stats import spearmanr
g = sns.jointplot("size", "total_bill", data=tips,
                  stat_func=spearmanr, color="m")
