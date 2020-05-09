import random

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


def Lagrange(b1, b2):
    B1 = inner_product(b1, b1)
    B2 = inner_product(b2, b2)

    while True:
        if B2 < B1:
            b1, b2 = b2, b1
            B1 = B2
        mu = round(inner_product(b1, b2) / B1)
        b2 = vecminus(b2, scalarMul(b1, mu))
        B2 = inner_product(b2, b2)

        if B2 >= B1:
            break
    return b1, b2

# example
# output ([-3, 0], [-1, -3])
vs = [[1, 1, 1], [-1, 0, 2], [3, 5, 6]]

n = 2
vs = [([random.randrange(1, 10000) for i in range(n)]) for j in range(n)]

print(vs)
l = Lagrange(vs[0], vs[1])

print(l)
