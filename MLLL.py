import random
import time

# vs == [0, 0, ..., 0] 
def is_zero(vs):
    for i in vs:
        if i != 0:
            return False
    return True

# <vs, ws>
def inner_product(vs, ws):
    res = 0
    for i in range(len(vs)):
        res += vs[i] * ws[i]
    return res

# vs + ws
def vecplus(vs, ws):
    res = []
    for i in range(len(vs)):
        res.append(vs[i] + ws[i])
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

# Modified LLL
def lll(bs):
    if len(bs) == 0:
        return []

    assert len(bs) >= len(bs[0])

    n = len(bs[0])
    h = len(bs)
    bstars = [0] * h
    norm_bstars = [0] * h
    mus = [([0] * h) for i in range(h)]
    delta = 0.75

    z = h - 1
    g = 0
    while g <= z:
        if is_zero(bs[g]):
            if g < z:
                bs[g], bs[z] = bs[z], bs[g]
            z -= 1

        # gram schmidt
        bstars[g] = bs[g]
        for j in range(g):
            mus[g][j] = inner_product(bs[g], bstars[j]) / norm_bstars[j]
            bstars[g] = vecminus(bstars[g], scalarMul(bstars[j], mus[g][j]))
        norm_bstars[g] = inner_product(bstars[g], bstars[g])
        mus[g][g] = 1

        if g == 0:
            g = 1 
        else:
            l = g
            k = g
            startagain = False
            while k <= l and not startagain:
                # size reduce (k, k - 1)
                if abs(mus[k][k - 1]) > 0.5:
                    q = round(mus[k][k - 1])
                    bs[k] = vecminus(bs[k], scalarMul(bs[k - 1], q))
                    for L in range(k):
                        mus[k][L] -= q * mus[k - 1][L]
                nu = mus[k][k - 1]
                B = norm_bstars[k] + nu * nu * norm_bstars[k - 1]

                # Lovasz condition
                if B >= delta * norm_bstars[k - 1]:
                    for j in reversed(range(k - 1)):
                        # size reduce (k, j)
                        if abs(mus[k][j]) > 0.5:
                            q = round(mus[k][j])
                            bs[k] = vecminus(bs[k], scalarMul(bs[j], q))
                            for L in range(j + 1):
                                mus[k][L] -= q * mus[j][L]
                    k += 1
                else:
                    if is_zero(bs[k]):
                        if k < z:
                            bs[k], bs[z] = bs[z], bs[k]
                        z -= 1
                        g = k
                        startagain = True
                    else:
                        bs[k - 1], bs[k] = bs[k], bs[k - 1]
                        for j in range(k - 1):
                            mus[k][j], mus[k - 1][j] = mus[k - 1][j], mus[k][j]
                        if B != 0:
                            if norm_bstars[k] == 0:
                                # remove lineary dependent
                                norm_bstars[k - 1] = B
                                bstars[k - 1] = scalarMul(bstars[k - 1], nu)
                                mus[k][k - 1] = 1 / nu
                                for i in range(k + 1, l + 1):
                                    mus[i][k - 1] = mus[i][k - 1] / nu
                            else:
                                t = norm_bstars[k - 1] / B
                                mus[k][k - 1] = nu * t
                                w = bstars[k - 1]
                                bstars[k - 1] = vecplus(bstars[k], scalarMul(w, nu))
                                norm_bstars[k - 1] = B
                                if k <= l:
                                    bstars[k] = vecminus(scalarMul(w, norm_bstars[k] / B), scalarMul(bstars[k], mus[k][k - 1]))
                                    norm_bstars[k] = norm_bstars[k] * t

                                for i in range(k + 1, l + 1):
                                    t = mus[i][k]
                                    mus[i][k] = mus[i][k - 1] - nu * t
                                    mus[i][k - 1] = t + mus[k][k - 1] * mus[i][k]
                        else:
                            norm_bstars[k], norm_bstars[k - 1] = norm_bstars[k - 1], norm_bstars[k]
                            bstars[k], bstars[k - 1] = bstars[k - 1], bstars[k]
                            for i in range(k + 1, l + 1):
                                mus[i][k], mus[i][k - 1] = mus[i][k - 1], mus[i][k]
                        k = max(k - 1, 1)
            if not startagain:
                g += 1

    return [list(map(int, b)) for b in bs]

vs = []
for i in range(6):
    xx = []
    for j in range(4):
        n = random.randrange(-300, 300)
        xx.append(n)
    vs.append(xx[::])
print(vs)
n = len(vs)

s = time.time()
l = lll(vs)
t = time.time() - s

print(l)
print("dim =", n)
print("time =", t, "[s]")
