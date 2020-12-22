import random
import time

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

# output: vs  (|π_k(v)|^2 <= R_{n+1-k}^2)
def coeff_vec(bs, Rs):
    n = len(bs)
    bstars, mus = gram_schmidt(bs)
    norm_bstars = [0] * n
    for i in range(n):
        norm_bstars[i] = inner_product(bstars[i], bstars[i])

    sigmas = [([0] * n) for i in range(n + 1)]

    # r[-1], r[0], r[1], ..., r[n - 1]
    rs = list(range(1, n + 1))
    rs[-1] = 0

    # rhos[k] == |π_k(v)|^2
    rhos = [0] * (n + 1)

    # v = sum (vs[i] * bs[i])
    vs = [0] * n
    vs[0] = 1

    cs = [0] * n
    ws = [0] * n

    # v[cnt + 1] = v[cnt + 2] = ... = 0
    # v[cnt], ..., v[0] != 0
    cnt = 0
    
    k = 0
    while True:
        # |π_k(v)|^2 = sum^{n}_{j = k} ( v_j + sum^{n}_{i = j+1} μ_{i,j} v_i )^2 |b^*|^2
        rhos[k] = rhos[k + 1] + (vs[k] - cs[k]) * (vs[k] - cs[k]) * norm_bstars[k]
        if rhos[k] <= Rs[n - k - 1]:
            if k == 0:
                return vs
            k -= 1
            rs[k - 1] = max(rs[k - 1], rs[k])
            for i in reversed(range(k + 1, rs[k] + 1)):
                # sigmas[i][k] = sum^{n}_{h = i} μ_{h,j} v_h
                sigmas[i][k] = sigmas[i + 1][k] + mus[i][k] * vs[i]
            cs[k] = -sigmas[k + 1][k]
            vs[k] = round(cs[k])
            ws[k] = 1
        else:
            k += 1
            if k == n:
                return []
            rs[k - 1] = k
            if k >= cnt:
                cnt = k
                vs[k] += 1
            else:
                if vs[k] > cs[k]:
                    vs[k] -= ws[k]
                else:
                    vs[k] += ws[k]
                ws[k] += 1

n = 10
vs = [([random.randrange(-200, 200) for i in range(n)]) for j in range(n)]
print("vs:", vs)

s = time.time()
l = coeff_vec(vs, [49500] * n)
t = time.time() - s

print("coeff:", l)
v = scalarMul(vs[0], l[0])
for i in range(1, n):
    v = vecplus(v, scalarMul(vs[i], l[i]))

print("v:", v)
print("dim =", n)
print("time =", t, "[s]")
