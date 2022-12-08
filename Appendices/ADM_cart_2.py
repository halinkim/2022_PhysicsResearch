import pickle

SIZE = 101
CENTER = SIZE // 2
dx = 0.01
dt = 0.1
M = 0.2

to_ij = [[0, 0], [0, 1], [0, 2], [1, 1], [1, 2], [2, 2]]
_dx = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

def sqrt(x):
return x ** .5

def to_ind(x, y):
return to_ij.index(sorted([x, y]))


def dist(x, y, z):
return (x * x + y * y + z * z) ** .5


def psi(r):
return 1 + M / 2 / r


def lapse(r):
return (2 * r - M) / (2 * r + M)


def dNx(ind, x, y, z):
if ind == 1:
x, y = y, x
elif ind == 2:
x, z = z, x
return 4 * M * x / ((M + 2 * sqrt(x ** 2 + y ** 2 + z ** 2)) ** 2 * sqrt(x ** 2 + y ** 2 + z ** 2))
def dNxx(ind, x, y, z):
if ind == 3:
x, y = y, x
elif ind == 5:
x, z = z, x
return 4 * M * (-x ** 2 * (M + 2 * sqrt(x ** 2 + y ** 2 + z ** 2)) * (x ** 2 + y ** 2 + z ** 2) ** (3 / 2) - 4 * x ** 2 * (
x ** 2 + y ** 2 + z ** 2) ** 2 + (M + 2 * sqrt(x ** 2 + y ** 2 + z ** 2)) * (
x ** 2 + y ** 2 + z ** 2) ** (5 / 2)) / (
(M + 2 * sqrt(x ** 2 + y ** 2 + z ** 2)) ** 3 * (x ** 2 + y ** 2 + z ** 2) ** 3)
def dNxy(ind, x, y, z):
if ind == 2:
y, z = z, y
elif ind == 4:
x, z = z, x
return 4 * M * x * y * (-(M + 2 * sqrt(x ** 2 + y ** 2 + z ** 2)) * (x ** 2 + y ** 2 + z ** 2) - 4 * (
x ** 2 + y ** 2 + z ** 2) ** (3 / 2)) / (
(M + 2 * sqrt(x ** 2 + y ** 2 + z ** 2)) ** 3 * (x ** 2 + y ** 2 + z ** 2) ** (5 / 2))
def ddlapse(ind, r, i, j, k):
_i, _j = to_ij[ind]
x, y, z = (i - CENTER) * dx, (j - CENTER) * dx, (k - CENTER) * dx
if _i == _j:
return dNxx(ind, x, y, z) - sum(Christoffel[_k][ind][i][j][k] * dNx(_k, x, y, z) for _k in range(3))
else:
return dNxy(ind, x, y, z) - sum(Christoffel[_k][ind][i][j][k] * dNx(_k, x, y, z) for _k in range(3))
# if ind == 0:
#     return -16 * M / (2 * r + M) ** 3 - 4 * M / (2 * r + M) ** 2 * Christoffel[0][0][i][j][k]


def transposeMatrix(m):
return list(map(list, zip(*m)))


def getMatrixMinor(m, i, j):
return [row[:j] + row[j + 1:] for row in (m[:i] + m[i + 1:])]


def getMatrixDeternminant(m):
if len(m) == 2:
return m[0][0] * m[1][1] - m[0][1] * m[1][0]
determinant = 0
for c in range(len(m)):
determinant += ((-1) ** c) * m[0][c] * getMatrixDeternminant(getMatrixMinor(m, 0, c))
return determinant


def getMatrixInverse(m):
determinant = getMatrixDeternminant(m)
if determinant == 0:
print('Warning : determinant is zero.')
return [[0] * 3 for _ in range(3)]
if len(m) == 2:
return [[m[1][1] / determinant, -1 * m[0][1] / determinant],
[-1 * m[1][0] / determinant, m[0][0] / determinant]]
cofactors = []
for r in range(len(m)):
cofactorRow = []
for c in range(len(m)):
minor = getMatrixMinor(m, r, c)
cofactorRow.append(((-1) ** (r + c)) * getMatrixDeternminant(minor))
cofactors.append(cofactorRow)
cofactors = transposeMatrix(cofactors)
for r in range(len(cofactors)):
for c in range(len(cofactors)):
cofactors[r][c] = cofactors[r][c] / determinant
return cofactors

initial = 0
tc = 0
errest = 1
evolve = 0
for _tt in range(1):
print('Step :', tc)
# 0rr, 1rtheta, 2rphi, 3thetatheta, 4thetaphi, 5phiphi

if initial:
gamma = [[[[0] * SIZE for _i in range(SIZE)] for _j in range(SIZE)] for _k in range(6)]
K = [[[[0] * SIZE for _i in range(SIZE)] for _j in range(SIZE)] for _k in range(6)]
else:
with open("gamma_%d.txt" % tc, "rb") as f:
gamma = pickle.load(f)
with open("K_%d.txt" % tc, "rb") as f:
K = pickle.load(f)

gamma_inv = [[[[0] * SIZE for _i in range(SIZE)] for _j in range(SIZE)] for _k in range(6)]
pgamma = [[[[[0] * SIZE for _i in range(SIZE)] for _j in range(SIZE)] for _k in range(6)] for _ in range(3)]
meanK = [[[0] * SIZE for _i in range(SIZE)] for _j in range(SIZE)]
R = [[[[0] * SIZE for _i in range(SIZE)] for _j in range(SIZE)] for _k in range(6)]

r0 = M / 2
if initial:
print('Set Inital Gamma')
for i in range(SIZE):
for j in range(SIZE):
for k in range(SIZE):
if i == j == k == CENTER:
continue
else:
r = dist(i - CENTER, j - CENTER, k - CENTER) * dx
gamma[0][i][j][k] = psi(r) ** 4
gamma[3][i][j][k] = psi(r) ** 4
gamma[5][i][j][k] = psi(r) ** 4
with open('gamma_%d.txt' % (0), 'wb') as f:
pickle.dump(gamma, f)
with open('K_%d.txt' % (0), 'wb') as f:
pickle.dump(K, f)

print('Get Inverse Gamma')
for i in range(SIZE):
for j in range(SIZE):
for k in range(SIZE):
if i == j == k == CENTER:
continue

_subgamma = [[gamma[0][i][j][k], gamma[1][i][j][k], gamma[2][i][j][k]],
[gamma[1][i][j][k], gamma[3][i][j][k], gamma[4][i][j][k]],
[gamma[2][i][j][k], gamma[4][i][j][k], gamma[5][i][j][k]]]
_gammainv = getMatrixInverse(_subgamma)
for _ind in range(6):
gamma_inv[_ind][i][j][k] = _gammainv[to_ij[_ind][0]][to_ij[_ind][1]]

print('Get Mean K')
for i in range(SIZE):
for j in range(SIZE):
for k in range(SIZE):
if i == j == k == CENTER:
continue
meanK[i][j][k] = sum(
gamma_inv[to_ind(_i, _j)][i][j][k] * K[to_ind(_i, _j)][i][j][k] for _i in range(3) for _j in range(3))

print('Get Partial Gamma')
for i in range(SIZE):
for j in range(SIZE):
for k in range(SIZE):
if abs(i - CENTER) + abs(j - CENTER) + abs(k - CENTER) <= 1:
continue
if i == 0 or i == SIZE - 1 or j == 0 or j == SIZE - 1 or k == 0 or k == SIZE - 1:
continue

for _sup in range(6):
# _ipgamma = [(gamma[_sup][i + 1][j][k] - gamma[_sup][i - 1][j][k]) / dx / 2,
#             (gamma[_sup][i][j + 1][k] - gamma[_sup][i][j - 1][k]) / dx / 2,
#             (gamma[_sup][i][j][k + 1] - gamma[_sup][i][j][k - 1]) / dx / 2]
for _sub in range(3):
# pgamma[_sub][_sup][i][j][k] = _ipgamma[_sub]
# _case = 0
# _coord = [i, j, k]
# _fixed = _coord[_sub]
# _coordp = _coord[:]
# _coordp[_sub] += 1
# _coordm = _coord[:]
# _coordm[_sub] -= 1
# if _fixed == SIZE - 1 or _coordp == [CENTER, CENTER, CENTER]:
#     _case = 2
# elif _fixed == 0 or _coordm == [CENTER, CENTER, CENTER]:
#     _case = 1
pgamma[_sub][_sup][i][j][k] = (gamma[_sup][i + _dx[_sub][0]][j + _dx[_sub][1]][k + _dx[_sub][2]]
- gamma[_sup][i - _dx[_sub][0]][j - _dx[_sub][1]][k - _dx[_sub][2]]) / dx / 2
# elif _case == 1:
#     pgamma[_sub][_sup][i][j][k] = (2 * gamma[_sup][i + 2 * _dx[_sub][0]][j + 2 * _dx[_sub][1]][k + 2 * _dx[_sub][2]]
#                                        - gamma[_sup][i + 3 * _dx[_sub][0]][j + 3 * _dx[_sub][1]][k + 3 * _dx[_sub][2]]
#                                        - 2 * gamma[_sup][i][j][k]
#                                        + gamma[_sup][i + _dx[_sub][0]][j + _dx[_sub][1]][k + _dx[_sub][2]]) / dx / 2
# else:
#     pgamma[_sub][_sup][i][j][k] = (- 2 * gamma[_sup][i - 2 * _dx[_sub][0]][j - 2 * _dx[_sub][1]][
#                                            k - 2 * _dx[_sub][2]]
#                                        + gamma[_sup][i - 3 * _dx[_sub][0]][j - 3 * _dx[_sub][1]][
#                                            k - 3 * _dx[_sub][2]]
#                                        + 2 * gamma[_sup][i][j][k]
#                                        - gamma[_sup][i - _dx[_sub][0]][j - _dx[_sub][1]][
#                                            k - _dx[_sub][2]]) / dx / 2

if initial:
with open('Pgamma.txt', 'wb') as f:
pickle.dump(pgamma, f)

print('Get Christoffel')
Christoffel = [[[[[0] * SIZE for _i in range(SIZE)] for _j in range(SIZE)] for _k in range(6)] for _l in range(3)]
for i in range(SIZE):
print('Christoffel for : ', i)
for j in range(SIZE):
for k in range(SIZE):
if abs(i - CENTER) + abs(j - CENTER) + abs(k - CENTER) <= 1:
continue
if i == 0 or i == SIZE - 1 or j == 0 or j == SIZE -1 or k == 0 or k == SIZE - 1:
continue

for _sup in range(3):
for _sub in range(6):
# _subgamma = [[gamma[0][i][j][k], gamma[1][i][j][k], gamma[2][i][j][k]],
#              [gamma[1][i][j][k], gamma[3][i][j][k], gamma[4][i][j][k]],
#              [gamma[2][i][j][k], gamma[4][i][j][k], gamma[5][i][j][k]]]
# _gammainv = getMatrixInverse(_subgamma)
_subi, _subj = to_ij[_sub]
Christoffel[_sup][_sub][i][j][k] = sum(
gamma_inv[to_ind(l, _sup)][i][j][k] * (pgamma[_subi][to_ind(l, _subj)][i][j][k] +
pgamma[_subj][to_ind(_subi, l)][i][j][k] -
pgamma[l][to_ind(_subi, _subj)][i][j][k]) for l in
range(3)) / 2
if initial:
with open('Christoffel.txt', 'wb') as f:
pickle.dump(Christoffel, f)

print('Get Ricci Tensor')

def pChristoffel(sup, sub, i, j, k, x):
# _case = 0
# _coord = [i, j, k]
# _fixed = _coord[x]
# _coordp = _coord[:]
# _coordp[x] += 1
# _coordm = _coord[:]
# _coordm[x] -= 1
# if _fixed == SIZE - 1 or _coordp == [CENTER, CENTER, CENTER]:
#     _case = 2
# elif _fixed == 0 or _coordm == [CENTER, CENTER, CENTER]:
#     _case = 1

# if _case == 0:
return (Christoffel[sup][sub][i + _dx[x][0]][j + _dx[x][1]][k + _dx[x][2]] - Christoffel[sup][sub][i - _dx[x][0]][j - _dx[x][1]][k - _dx[x][2]]) / dx / 2
# elif _case == 1:
#     return (2 * Christoffel[sup][sub][i + 2 * _dx[x][0]][j + 2 * _dx[x][1]][k + 2 * _dx[x][2]]
#             - Christoffel[sup][sub][i + 3 * _dx[x][0]][j + 3 * _dx[x][1]][k + 3 * _dx[x][2]]
#             - 2 * Christoffel[sup][sub][i][j][k]
#             + Christoffel[sup][sub][i + _dx[x][0]][j + _dx[x][1]][k + _dx[x][2]])
# else:
#     return (- 2 * Christoffel[sup][sub][i - 2 * _dx[x][0]][j - 2 * _dx[x][1]][k - 2 * _dx[x][2]]
#             + Christoffel[sup][sub][i - 3 * _dx[x][0]][j - 3 * _dx[x][1]][k - 3 * _dx[x][2]]
#             + 2 * Christoffel[sup][sub][i][j][k]
#             - Christoffel[sup][sub][i - _dx[x][0]][j - _dx[x][1]][k - _dx[x][2]])



for i in range(SIZE):
print('Ricci Tensor for :', i)
for j in range(SIZE):
for k in range(SIZE):
if abs(i - CENTER) + abs(j - CENTER) + abs(k - CENTER) <= 2:
continue
# if i == j == k == CENTER or i == 0 or i == SIZE - 1 or j == 0 or j == SIZE - 1 or k == 0 or k == SIZE - 1:
#     continue
if i <= 1 or i >= SIZE - 2 or j <= 1 or j >= SIZE - 2 or k <= 1 or k >= SIZE - 2:
continue

for _ind in range(6):
_i, _j = to_ij[_ind]
R[_ind][i][j][k] = sum(pChristoffel(_k, _ind, i, j, k, _k) for _k in range(3)) \
- sum(pChristoffel(_k, to_ind(_i, _k), i, j, k, _j) for _k in range(3))\
+ sum(
Christoffel[_k][to_ind(_i, _j)][i][j][k] * Christoffel[_l][to_ind(_k, _l)][i][j][k] for _k in
range(3) for _l in range(3)) \
- sum(
Christoffel[_l][to_ind(_i, _k)][i][j][k] * Christoffel[_k][to_ind(_l, _j)][i][j][k] for _k in
range(3) for _l in range(3))

if errest:
print('Error est')
print('mean R')
meanR = [[[0] * SIZE for _i in range(SIZE)] for _j in range(SIZE)]
for i in range(SIZE):
for j in range(SIZE):
for k in range(SIZE):
if i == j == k == CENTER:
continue
meanR[i][j][k] = sum(
gamma_inv[to_ind(_i, _j)][i][j][k] * R[to_ind(_i, _j)][i][j][k] for _i in range(3) for _j in range(3))
print('KK')
KK = [[[0] * SIZE for _i in range(SIZE)] for _j in range(SIZE)]
for i in range(SIZE):
for j in range(SIZE):
for k in range(SIZE):
if i == j == k == CENTER:
continue
KK[i][j][k] = sum(sum(gamma_inv[to_ind(_i, _mu)][i][j][k] * gamma_inv[to_ind(_j, _nu)][i][j][k] *
K[to_ind(_mu, _nu)][i][j][k] for _mu in range(3) for _nu in range(3))
* K[to_ind(_i, _j)][i][j][k] for _i in range(3) for _j in range(3))
err = [[[0] * SIZE for _i in range(SIZE)] for _j in range(SIZE)]
for i in range(SIZE):
for j in range(SIZE):
for k in range(SIZE):
if i == j == k == CENTER:
continue
err[i][j][k] = meanR[i][j][k] + (meanK[i][j][k]) ** 2 - KK[i][j][k]
with open('error_%d.txt' % (tc), 'wb') as f:
pickle.dump(err, f)




# R[_ind][i][j][k] = sum((Christoffel[_k][_ind][i + _dx[_k][0]][j + _dx[_k][1]][k + _dx[_k][2]] -
#                         Christoffel[_k][_ind][i][j][k]) / dx for _k in range(3)) \
#                    - sum((Christoffel[_k][to_ind(_i, _k)][i + _dx[_j][0]][j + _dx[_j][1]][
#                               k + _dx[_j][2]] - Christoffel[_k][to_ind(_i, _k)][i][j][k]) / dx for _k in
#                          range(3)) \
#                    + sum(
#     Christoffel[_k][to_ind(_i, _j)][i][j][k] * Christoffel[_l][to_ind(_k, _l)][i][j][k] for _k in
#     range(3) for _l in range(3)) \
#                    - sum(
#     Christoffel[_l][to_ind(_i, _k)][i][j][k] * Christoffel[_k][to_ind(_l, _j)][i][j][k] for _k in
#     range(3) for _l in range(3))

# with open("gamma.txt", 'wb') as f:
#     pickle.dump(gamma, f)
# with open('gamma_inv.txt', 'wb') as f:
#     pickle.dump(gamma_inv, f)
# with open('pgamma.txt', 'wb') as f:
#     pickle.dump(pgamma, f)
# with open('K.txt', 'wb') as f:
#     pickle.dump(K, f)
# with open('meanK.txt', 'wb') as f:
#     pickle.dump(meanK, f)
# with open('R.txt', 'wb') as f:
#     pickle.dump(R, f)
# with open('Christoffel.txt', 'wb') as f:
#     pickle.dump(Christoffel, f)
#
if evolve:
print('Evolve')
_gamma = [[[[0] * SIZE for _i in range(SIZE)] for _j in range(SIZE)] for _k in range(6)]
_K = [[[[0] * SIZE for _i in range(SIZE)] for _j in range(SIZE)] for _k in range(6)]

for i in range(SIZE):
print('Evolve for :', i)
for j in range(SIZE):
for k in range(SIZE):
# if i == j == k == CENTER:
#     continue
# if i == 0 or i == SIZE - 1 or j == 0 or j == SIZE - 1 or k == 0 or k == SIZE - 1:
#     continue
if i <= 1 or i >= SIZE - 2 or j <= 1 or j >= SIZE - 2 or k <= 1 or k >= SIZE - 2:
for ind in range(6):
_gamma[ind][i][j][k] = gamma[ind][i][j][k]
_K[ind][i][j][k] = K[ind][i][j][k]
continue
r = dist(i - CENTER, j - CENTER, k - CENTER) * dx
if r <= r0:
for ind in range(6):
_gamma[ind][i][j][k] = gamma[ind][i][j][k]
_K[ind][i][j][k] = K[ind][i][j][k]
continue
for ind in range(6):
_i, _j = to_ij[ind]
_gamma[ind][i][j][k] = gamma[ind][i][j][k] - 2 * lapse(r) * K[ind][i][j][k] * dt
_K[ind][i][j][k] = K[ind][i][j][k] - dt * (-ddlapse(ind, r, i, j, k) + lapse(r) * (R[ind][i][j][k]
+ meanK[i][j][k] *
K[ind][i][j][
k] - 2 * sum(
K[to_ind(_i, _k)][i][j][k] * sum(
gamma_inv[to_ind(_k, _u)][i][j][k] * K[to_ind(_u, _j)][i][j][k] for _u in range(3)) for
_k in range(3))))

with open('gamma_%d.txt' % (tc + 1), 'wb') as f:
pickle.dump(_gamma, f)
with open('K_%d.txt' % (tc + 1), 'wb') as f:
pickle.dump(_K, f)
tc += 1


