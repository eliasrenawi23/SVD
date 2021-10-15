import math
import sys
from timeit import timeit
import numpy as np


def jacobi_eigenvalue(n, a, it_max):
    #  Parameters:
    #    Input, integer N, the order of the matrix.
    #    Input, real A(N,N), the matrix, which must be square, real,
    #    and symmetric.
    #    Input, integer IT_MAX, the maximum number of iterations.
    #    Output, real V(N,N), the matrix of eigenvectors.
    #    Output, real D(N), the eigenvalues, in descending order.
    #    Output, integer IT_NUM, the total number of iterations.
    #    Output, integer ROT_NUM, the total number of rotations.
    v = np.zeros([n, n])  # eigenvectors
    d = np.zeros(n)  ##eigenvalue

    for j in range(0, n):
        v[j, j] = 1.0
    for i in range(0, n):
        d[i] = a[i, i]

    bw = np.zeros(n)
    zw = np.zeros(n)
    w = np.zeros(n)

    for i in range(0, n):
        bw[i] = d[i]

    it_num = 0
    rot_num = 0

    while (it_num < it_max):
        print("it_num=", it_num)

        it_num = it_num + 1
        #
        #  The convergence threshold is based on the size of the elements in
        #  the strict upper triangle of the matrix.
        #
        thresh = 0.0
        for j in range(0, n):
            for i in range(0, j):
                thresh = thresh + a[i, j] ** 2

        thresh = np.sqrt(thresh) / float(4 * n)

        if (thresh == 0.0):
            break

        for p in range(0, n):
            for q in range(p + 1, n):

                gapq = 10.0 * abs(a[p, q])
                termp = gapq + abs(d[p])
                termq = gapq + abs(d[q])
                #
                #  Annihilate tiny offdiagonal elements.
                #
                if (4 < it_num and termp == abs(d[p]) and termq == abs(d[q])):

                    a[p, q] = 0.0
                #
                #  Otherwise, apply a rotation.
                #
                elif (thresh <= abs(a[p, q])):

                    h = d[q] - d[p]
                    term = abs(h) + gapq

                    if (term == abs(h)):
                        t = a[p, q] / h
                    else:
                        theta = 0.5 * h / a[p, q]
                        t = 1.0 / (abs(theta) + np.sqrt(1.0 + theta * theta))
                        if (theta < 0.0):
                            t = - t

                    c = 1.0 / np.sqrt(1.0 + t * t)
                    s = t * c
                    tau = s / (1.0 + c)
                    h = t * a[p, q]
                    #
                    #  Accumulate corrections to diagonal elements.
                    #
                    zw[p] = zw[p] - h
                    zw[q] = zw[q] + h
                    d[p] = d[p] - h
                    d[q] = d[q] + h

                    a[p, q] = 0.0
                    #
                    #  Rotate, using information from the upper triangle of A only.
                    #
                    for j in range(0, p):
                        g = a[j, p]
                        h = a[j, q]
                        a[j, p] = g - s * (h + g * tau)
                        a[j, q] = h + s * (g - h * tau)

                    for j in range(p + 1, q):
                        g = a[p, j]
                        h = a[j, q]
                        a[p, j] = g - s * (h + g * tau)
                        a[j, q] = h + s * (g - h * tau)

                    for j in range(q + 1, n):
                        g = a[p, j]
                        h = a[q, j]
                        a[p, j] = g - s * (h + g * tau)
                        a[q, j] = h + s * (g - h * tau)
                    #
                    #  Accumulate information in the eigenvector matrix.
                    #
                    for j in range(0, n):
                        g = v[j, p]
                        h = v[j, q]
                        v[j, p] = g - s * (h + g * tau)
                        v[j, q] = h + s * (g - h * tau)

                    rot_num = rot_num + 1

        for i in range(0, n):
            bw[i] = bw[i] + zw[i]
            d[i] = bw[i]
            zw[i] = 0.0
    #
    #  Restore upper triangle of input matrix.
    #
    for j in range(0, n):
        for i in range(0, j):
            a[i, j] = a[j, i]
    #
    #  Ascending sort the eigenvalues and eigenvectors.
    #
    for k in range(0, n - 1):
        m = k
        for l in range(k + 1, n):
            if (d[l] < d[m]):
                m = l

        if (k != m):

            t = d[m]
            d[m] = d[k]
            d[k] = t

            for i in range(0, n):
                w[i] = v[i, m]
                v[i, m] = v[i, k]
                v[i, k] = w[i]

    return v, d


def name():
    A = np.loadtxt("C:\data_8.txt", dtype=np.float, delimiter=',')

    # A = np.array([[1,0,1,0],[0,1,0,1]])
    np.set_printoptions(threshold=sys.maxsize)  # to print all of the matrix

    # print(A)
    r, c = A.shape
    f = open("results.txt", 'w')
    segma = np.zeros([r, c])

    AAT = np.dot(A, A.T)

    print("-------------------------------------------AAT---------------------------------\n")
    # f.write("-------------------------------------------AAT ---------------------------\n")

    # print(AAT)

    vectors1, valu = jacobi_eigenvalue(len(AAT), AAT, 100)
    print("-------------------------------------------valu ---------------------------\n")
    for i in range(r):
        for j in range(c):
            if i == j and valu[i] != 0:
                segma[i, j] = math.sqrt(valu[i])
    # print(valu)
    f.write("------------------------------------------- The U matrix ---------------------------\n")
    print("------------------------------------------- U matrix  ---------------------------\n")
    # print(vectors1)
    f.write(str(vectors1) + "\n")
    # print("-------------------------------------------ATA ---------------------------\n")
    ATA = np.dot(A.T, A)

    # print(ATA)

    vectors2, valu2 = jacobi_eigenvalue(len(ATA), ATA, 100)
    print("-------------------------------------------valu and V matrix ---------------------------\n")
    f.write("-------------------------------------------The V matrix  ---------------------------\n")

    print(valu2)
    print(vectors2)
    f.write(str(vectors1) + "\n")
    # print("---------------------------------")
    # print(vectors2.transpose())
    for i in range(r):  # fill the segma matrix
        for j in range(c):
            if i == j and valu2[i] != 0:
                if segma[i, j] != math.sqrt(valu2[i]):
                    segma[i, j] = math.sqrt(valu2[i])
    f.write("-------------------------------------------The segma matrix  ---------------------------\n")
    f.write(str(segma))
    f.close()

    '''# print("the segma matrix is :\n", segma) # this is for check if we got the right answers 
    npsegma = np.array(segma)
    us = np.dot(vectors1, npsegma)
    usv = np.dot(us, vectors2)
    print(usv) # A =U*SEGMA* V^T'''


print("the total time is : ", timeit(name, number=1))
