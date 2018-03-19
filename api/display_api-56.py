planets = sns.load_dataset("planets")
ax = sns.violinplot(x="orbital_period", y="method",
                    data=planets[planets.orbital_period < 1000],
                    scale="width", palette="Set3")
