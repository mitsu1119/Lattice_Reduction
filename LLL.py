import random
import time

# <vs, ws>
def inner_product(vs, ws):
    res = 0
    for i in range(len(vs)):
        res += vs[i] * ws[i]
    return res

# vs - ws
def vecminus(vs, ws):
    res = []
    for i in range(len(vs)):
        res.append(vs[i] - ws[i])
    return res

# a * vs
def scalarMul(vs, a):
    res = []
    for i in range(len(vs)):
        res.append(vs[i] * a)
    return res

# vs -> vs*
def gram_schmidt(vs):
    n = len(vs)
    res = [0] * n
    mus = [([1] * n) for i in range(n)]
    res[0] = vs[0]
    for i in range(1, n):
        v = vs[i]
        for j in range(i):
            mus[i][j] = inner_product(vs[i], res[j]) / inner_product(res[j], res[j])
            v = vecminus(v, scalarMul(res[j], mus[i][j]))
        res[i] = v
    return res, mus

def lll(bs):
    n = len(bs)
    bstars, mus = gram_schmidt(bs)
    norm_bstars = [0] * n
    for i in range(n):
        norm_bstars[i] = inner_product(bstars[i], bstars[i])
    delta = 3 / 4

    k = 1
    while k < n:
        for j in reversed(range(k)):
            mu_kj = mus[k][j]
            if abs(mu_kj) > 0.5:
                q = round(mu_kj)
                bs[k] = vecminus(bs[k], scalarMul(bs[j], q))
                # bstars, mus = gram_schmidt(bs)
                for l in range(j + 1):
                    mus[k][l] = mus[k][l] - q * mus[j][l]

        if norm_bstars[k] >= (delta - mus[k][k - 1] ** 2) * norm_bstars[k - 1]:
            k += 1
        else:
            bs[k], bs[k - 1] = bs[k - 1], bs[k]
            mup = mus[k][k - 1]
            B = norm_bstars[k] + (mup ** 2) * norm_bstars[k - 1]
            mus[k][k - 1] = mup * norm_bstars[k - 1] / B
            norm_bstars[k] = norm_bstars[k] * norm_bstars[k - 1] / B
            norm_bstars[k - 1] = B
            for j in range(k - 1):
                mus[k - 1][j], mus[k][j] = mus[k][j], mus[k - 1][j]
            for j in range(k + 1, n):
                t = mus[j][k]
                mus[j][k] = mus[j][k - 1] - mup * t
                mus[j][k - 1] = t + mus[k][k - 1] * mus[j][k]
            k = max(k - 1, 1)
    return [list(map(int, b)) for b in bs]

# example in wikipedia
vs = [[1, 1, 1], [-1, 0, 2], [3, 5, 6]]

n = random.randrange(10, 200)
vs = [([random.randrange(1, 10000) for i in range(n)]) for j in range(n)]

s = time.time()
l = lll(vs)
t = time.time() - s

print(l)
print("dim =", n)
print("time =", t, "[s]")
