ax = sns.violinplot(x="orbital_period", y="method",
                    data=planets[planets.orbital_period < 1000],
                    cut=0, scale="width", palette="Set3")
