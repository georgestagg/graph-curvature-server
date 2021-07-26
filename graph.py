import web  # type: ignore[import]
from numpy.linalg import eigvalsh  # type: ignore[import]
import numpy as np  # type: ignore[import]
import json
from scipy.optimize import linprog  # type: ignore[import]
import scipy  # type: ignore[import]
from curvature import inf, normalised_unweighted_curvature, non_normalised_unweighted_curvature  # type: ignore[import]


def Amat(n, m):
    res = [[0 for i in range(n+m)] for i in range(n*m)]
    for i in range(n):
        for j in range(i*m, (i+1)*m):
            res[j][i] = 1
    for i in range(n, n+m):
        for j in range(n*m):
            if (j % m) == (i-n):
                res[j][i] = 1
    return res


def eta(n, m):
    res = [0 for i in range(n+m)]
    for i in range(1, n):
        res[i] = -1.0/(n-1)
    for i in range(n+1, n+m):
        res[i] = -1.0/(m-1)
    return res


def etap(n, m, p):
    res = [0 for i in range(n+m)]
    res[0] = -p
    res[n] = -p
    for i in range(1, n):
        res[i] = (p-1.0)/(n-1)
    for i in range(n+1, n+m):
        res[i] = (p-1.0)/(m-1)
    return res


def dist(i, j, A):
    if i == j:
        return 0
    if A[i][j] == 1:
        return 1
    for a in range(len(A)):
        if (A[i][a]+A[j][a]) == 2:
            return 2
    return 3


def d(x, y, A):
    res = []
    xnbs = [x]
    ynbs = [y]
    for i in range(len(A)):
        if A[x][i] == 1:
            xnbs.append(i)
        if A[y][i] == 1:
            ynbs.append(i)
    for i in xnbs:
        for j in ynbs:
            res.append(dist(i, j, A))
    return res


def ocurve(x, y, A):
    dx = sum(A[x])
    dy = sum(A[y])
    return 1+scipy.optimize.OptimizeResult.values(
        linprog(c=eta(dx+1, dy+1), A_ub=Amat(dx+1, dy+1), b_ub=d(x, y, A), bounds=(None, None)))[3]


def lazocurve(x, y, A, p):
    dx = sum(A[x])
    dy = sum(A[y])
    return 1+scipy.optimize.OptimizeResult.values(
        linprog(c=etap(dx+1, dy+1, p), A_ub=Amat(dx+1, dy+1), b_ub=d(x, y, A), bounds=(None, None)))[3]


def LLYcurv(x, y, A):
    d = max([sum(A[x]), sum(A[y])])
    return ((1.0*(d+1))/(1.0*d))*lazocurve(x, y, A, 1.0/(d+1))


def etanonnorm(n, m):
    res = [-1.0 for i in range(n+m)]
    res[0] = -m+1
    res[n] = -n+1
    return res


def nonnorm_ocurve(x, y, A):
    dx = sum(A[x])
    dy = sum(A[y])
    return dx+dy+scipy.optimize.OptimizeResult.values(
        linprog(c=etanonnorm(dx+1, dy+1), A_ub=Amat(dx+1, dy+1), b_ub=d(x, y, A), bounds=(None, None)))[3]


urls = (
  '/', 'index'
)
web.config.debug = False


def sign(x):
    return -1 if x < 0 else (1 if x > 0 else 0)


class index:
    def POST(self):
        web.header('Access-Control-Allow-Origin', '*')
        user_data = web.input()
        try:
            AM = json.loads(user_data.am)
            V = json.loads(user_data.v)
            t = json.loads(user_data.t)
        except Exception:
            return '["error0"]'

        if t == 0:
            try:
                ret = range(len(V))
            except Exception:
                return '["error3"]'
        if t == 1:
            try:
                ret = non_normalised_unweighted_curvature(AM, inf)
                ret = [sign(j) for j in ret]
            except Exception:
                return '["error1"]'
        if t == 2:
            try:
                ret = normalised_unweighted_curvature(AM, inf)
                ret = [sign(j) for j in ret]
            except Exception:
                return '["error4"]'
        if t == 3:
            try:
                ret = non_normalised_unweighted_curvature(AM, inf)
            except Exception:
                return '["error5"]'
        if t == 4:
            try:
                ret = normalised_unweighted_curvature(AM, inf)
            except Exception:
                return '["error6"]'
        if t == 5:
            try:
                dimn = json.loads(user_data.d)
                if(dimn == 0):
                    return '["error8"]'
                elif(dimn < 0):
                    return '["error8b"]'
            except Exception:
                return '["error9"]'
            try:
                ret = ret = non_normalised_unweighted_curvature(AM, dimn)
            except Exception:
                return '["error7"]'
        if t == 6:
            try:
                ret = dict()
                ret["AM"] = AM
                ret["ORC"] = [[0 for i in range(len(V))] for j in range(len(V))]

                for i in range(len(V)):
                    for j in range(len(V)):
                        if AM[i][j] == 1:
                            ret["ORC"][i][j] = ocurve(i, j, AM)
            except Exception:
                return '["error10"]'
        if t == 7:
            try:
                dimn = json.loads(user_data.d)
                if(dimn == 0):
                    return '["error8"]'
                elif(dimn < 0):
                    return '["error8b"]'
            except Exception:
                return '["error9"]'
            try:
                ret = normalised_unweighted_curvature(AM, dimn)
            except Exception:
                return '["error11"]'

        if t == 8:
            try:
                ret = dict()
                idlen = json.loads(user_data.idlen)
                ret["AM"] = AM
                ret["ORCI"] = [[0 for i in range(len(V))] for j in range(len(V))]
                if idlen < 0:
                    return '["error13a"]'
                elif idlen == 1:
                    return '["error13b"]'
                elif idlen > 1:
                    return '["error13c"]'
                for i in range(len(V)):
                    for j in range(len(V)):
                        if AM[i][j] == 1:
                            ret["ORCI"][i][j] = lazocurve(i, j, AM, idlen)
            except Exception:
                return '["error12"]'
        
        if t == 9:
            try:
                ret = dict()
                ret["AM"] = AM
                ret["LLYC"] = [[0 for i in range(len(V))] for j in range(len(V))]

                for i in range(len(V)):
                    for j in range(len(V)):
                        if AM[i][j] == 1:
                            ret["LLYC"][i][j] = LLYcurv(i, j, AM)
            except Exception:
                return '["error14"]'

        if t == 10:
            try:
                ret = dict()
                ret["AM"] = AM
                ret["NNLLYC"] = [[0 for i in range(len(V))] for j in range(len(V))]

                for i in range(len(V)):
                    for j in range(len(V)):
                        if AM[i][j] == 1:
                            ret["NNLLYC"][i][j] = nonnorm_ocurve(i, j, AM)
            except Exception:
                return '["error15"]'

        try:
            return json.dumps(ret)
        except Exception:
            return '["error2"]'


if __name__ == "__main__":
    app = web.application(urls, globals())
    app.run()
