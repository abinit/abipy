titanic = sns.load_dataset("titanic")
g = sns.factorplot("alive", col="deck", col_wrap=4,
                   data=titanic[titanic.deck.notnull()],
                   kind="count", size=2.5, aspect=.8)
