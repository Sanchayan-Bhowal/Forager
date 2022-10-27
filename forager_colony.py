import matplotlib.pyplot as plt
import numpy as np
import math

from scipy.spatial import distance

np.random.seed(0)
obedience = 1
n_max = 10000
v = 10
a = 5

h = 0.01

birth_rate = 0.025
death_rate = 1e-5

aging = 60

good = [np.array([-1, 1]), np.array([1.5, 1.5]),
        np.array([1.5, -1.5]), np.array([-1, -1])]
good_x = [blah[0] for blah in good]
good_y = [blah[1] for blah in good]
optimal = np.mean(good, axis=0)

eps = 0.01

n_foragers = 100
age = np.array([i+1 for i in range(n_foragers)], dtype=float)

colony = np.random.uniform(low=-a/2, high=a/2, size=2)

xs = []
ys = []
n_reached = 0
x = np.random.uniform(low=-0.25, high=0.25, size=n_foragers)+colony[0]
y = np.random.uniform(low=-0.25, high=0.25, size=n_foragers)+colony[1]
xs.append(x)
ys.append(y)
i = 0
colonies = [colony]
for j in range(n_max):
    if j % 50 == 0:
        print(j)
    # birth
    repab = len(list(x for x in age if 20 <= x <= aging))
    bab = np.random.poisson(lam=repab*birth_rate)
    x = np.append(x, np.random.uniform(
        low=-0.25, high=0.25, size=bab)+colony[0])
    y = np.append(y, np.random.uniform(
        low=-0.25, high=0.25, size=bab)+colony[1])
    n_foragers += bab
    age = np.append(age, np.zeros(bab, dtype=float))

    for n in range(n_foragers):
        # roaming
        if age[n] > aging:

            theta = np.random.uniform(
                low=-math.pi, high=math.pi, size=1)
            x[n] = x[n] + v * h * math.cos(theta)
            y[n] = y[n] + v * h * math.sin(theta)

            if x[n] > a:
                x[n] = 2.0*a-x[n]

            elif x[n] < -a:
                x[n] = -2.0*a-x[n]

            if y[n] > a:
                y[n] = 2.0*a-y[n]

            elif y[n] < -a:
                y[n] = -2.0*a-y[n]

            try:
                dist = distance.cdist(np.array([(x[n], y[n])]), good)
            except ValueError:
                break
            d = np.argmin(dist[0])
            min_d = dist[0][d]
            if min_d < eps:

                # obedience
                if np.random.uniform(low=0, high=1) < obedience:
                    x[np.argwhere(age <= aging)] -= colony[0]
                    y[np.argwhere(age <= aging)] -= colony[1]
                    colony = (colony*n_reached+good[d])/(n_reached+1)
                    colonies.append(colony)
                    xs.append(x)
                    ys.append(y)
                    print(colony)
                    x[np.argwhere(age <= aging)] += colony[0]
                    y[np.argwhere(age <= aging)] += colony[1]
                    n_reached += 1
                    good.pop(d)

    # death
    deads = np.ones(n_foragers)
    for n in range(n_foragers):
        try:
            if np.random.uniform(low=0, high=1) < 1-math.exp(-age[n]*death_rate):
                deads[n] = 0
        except OverflowError:
            deads[n] = 0
    x = np.delete(x, np.where(deads == 0))
    y = np.delete(y, np.where(deads == 0))
    age = np.delete(age, np.where(deads == 0))
    n_foragers -= np.sum(deads == 0)

    if(np.linalg.norm(colony-optimal) < eps):
        break
    i += 1
    age += 1


plt.xlim(-a, a)
plt.ylim(-a, a)

# plt.scatter(x, y, c='b', s=20, marker="^")
# plt.scatter(init_x, init_y, c='r', s=20, marker="o")
plt.scatter(good_x, good_y, color='g', marker="X", s=200)
print(colonies)
for k in range(len(xs)):
    plt.scatter(xs[k], ys[k], s=20/(k+1))

    try:
        dx = colonies[k+1][0]-colonies[k][0]
        dy = colonies[k+1][1]-colonies[k][1]
        plt.arrow(colonies[k][0], colonies[k][1], dx,
                  dy, width=0.01, head_width=0.05, length_includes_head=True)
    except IndexError:
        pass
plt.show()
