
import numpy as np
import math

f = open('results.txt', 'w')


def whale_walk(n_max=10000):
    death_rate = 1e-6
    birth_rate = 0.025

    aging = 60

    n_foragers = 100
    age = np.array([i+1 for i in range(n_foragers)], dtype=float)

    for _ in range(n_max):
        # birth
        repab = len(list(x for x in age if 20 <= x <= aging))
        bab = np.random.poisson(lam=repab*birth_rate)
        n_foragers += bab
        age = np.append(age, np.zeros(bab, dtype=float))

        # death
        deads = np.ones(n_foragers)
        for n in range(n_foragers):
            try:
                if np.random.uniform(low=0, high=1) < 1-math.exp(-age[n]*death_rate):
                    deads[n] = 0
            except OverflowError:
                deads[n] = 0
        age = np.delete(age, np.where(deads == 0))
        n_foragers -= np.sum(deads == 0)

        age += 1

    return n_foragers


for n in range(2000, 9000, 1000):
    f.write(str(n) + '\n')
    for _ in range(100):
        k = whale_walk(n)
        f.write(str(k)+',')
        f.flush()
    f.write('\n')
