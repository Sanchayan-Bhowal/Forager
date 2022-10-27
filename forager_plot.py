import forager

realisation = 100
deathrates = [1e-5, 9e-6, 8e-6, 7e-6, 6e-6, 5e-6, 4e-6, 3e-6, 2e-6]
f = open('results.txt', 'w')
for death_rate in deathrates:
    f.write(str(death_rate) + '\n')
    for _ in range(realisation):
        k = forager.whale_walk(death_rate, n_max=10000)
        f.write(str(k)+',')
        f.flush()
    f.write('\n')
f.close()
