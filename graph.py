import web
from numpy.linalg import eigvalsh
import numpy as np
import json
from sympy.matrices import *
from scipy.optimize import linprog
import scipy

def evs(A):
	V=eigvalsh(A)
	for i in range(len(V)):
		V[i] = format(V[i], 'f')
	return V

def intersect(a, b):
    return list(set(a) & set(b))

def one_sphere(A, i):
	n = len(A)
	RES = []
	for j in range(n):
		if A[i][j] == 1:
			RES.append(j)
	return RES

def two_sphere(A, i):
	M = np.array(A)
	M2 = np.dot(M, M)
	n = len(A)
	RES = []
	V = M2[i]
	for j in range(n):
		if V[j] != 0:
			if A[i][j] != 1:
				if i != j:
					RES.append(j)
	return RES

def two_ball(A, i):
	S1 = one_sphere(A, i)
	S2 = two_sphere(A, i)
	O = []
	O.append(i)
	O.extend(S1)
	O.extend(S2)
	n = len(O)
	RES = [[0 for x in range(n)] for x in range(n)]	
	for j in range(n):
		for k in range(n):
			RES[j][k] = A[O[j]][O[k]]
	return RES

def fourGammatwo(A, i):
	B = two_ball(A, i)
	n = len(B)
	RES = [[0 for x in range(n)] for x in range(n)]
	d = sum(B[0])
	S1 = one_sphere(B, 0)
	S2 = two_sphere(B, 0)
	E = [len(intersect(one_sphere(B, j), S1)) for j in range(n)]
	F = [len(intersect(one_sphere(B, j), S2)) for j in range(n)]
	RES[0][0] = d*(d+3)
	for j in S1:
		RES[j][j] = 5-d+3*F[j]+4*E[j]
		RES[0][j] = -3-d-F[j]  
		RES[j][0] = -3-d-F[j]
	for j in S2:
		RES[j][j] = E[j]
		RES[0][j] = E[j] 
		RES[j][0] = E[j]
	for j in S1:
		for k in S1:
			if j != k:
				RES[j][k] = 2 - 4*B[j][k]
	for j in S1:
		for k in S2:
			RES[j][k] = -2*B[j][k]
			RES[k][j] = -2*B[j][k]  
	return RES

def mu_two_ball(A, MU, i):
	S1 = one_sphere(A, i)
	S2 = two_sphere(A, i)
	O = []
	O.append(i)
	O.extend(S1)
	O.extend(S2)
	n = len(O)
	RES = [0 for x in range(n)]
	for j in range(n):
		RES[j] = MU[O[j]]
	return RES

def fourGammatwoMu(A, MU, i):
	B = two_ball(A, i)
	M = mu_two_ball(A, MU, i)
	n = len(B)
	RES = [[0 for x in range(n)] for x in range(n)]
	d = sum(B[0])
	mx = 1.0*M[0]
	S1 = one_sphere(B, 0)
	S2 = two_sphere(B, 0)
	RES[0][0] = (1.0*d*d)/(mx**2)+(3.0/(mx))*sum([1.0/(1.0*M[y]) for y in S1])
	for j in S1:
		RES[j][j] = 3.0/(mx*M[j])+2.0/(mx**2)-(1.0*d)/(mx**2)+(3.0/(mx*M[j]))*sum([B[j][z] for z in S2])+(1.0/mx)*sum([(1.0/(1.0*M[z])+3.0/(1.0*M[j]))*B[j][z] for z in S1])
		RES[0][j] = -3.0/(mx*M[j])-(1.0*d)/(mx**2)-(1.0/(mx*M[j]))*sum([B[j][z] for z in S2])-(1.0/mx)*sum([(-1.0/(1.0*M[z])+1.0/(1.0*M[j]))*B[j][z] for z in S1])
		RES[j][0] = -3.0/(mx*M[j])-(1.0*d)/(mx**2)-(1.0/(mx*M[j]))*sum([B[j][z] for z in S2])-(1.0/mx)*sum([(-1.0/(1.0*M[z])+1.0/(1.0*M[j]))*B[j][z] for z in S1])
	for j in S2:
		RES[j][j] = (1.0/mx)*sum((1.0/M[y])*B[y][j] for y in S1)
		RES[0][j] = (1.0/mx)*sum((1.0/M[y])*B[y][j] for y in S1)
		RES[j][0] = (1.0/mx)*sum((1.0/M[y])*B[y][j] for y in S1)
	for j in S1:
		for k in S1:
			if j != k:
				RES[j][k] = 2.0/(mx**2)-(1.0/mx)*(2.0/(1.0*M[j]) + 2.0/(1.0*M[k]))*B[j][k] 
	for j in S1:
		for k in S2:
			RES[j][k] = (-2.0*B[j][k])/(mx*M[j])
			RES[k][j] = (-2.0*B[j][k])/(mx*M[j])
	return RES

def fourGammatwoNorm(A, i):
	MU = [sum([A[j][k] for k in range(len(A))]) for j in range(len(A))]
	return fourGammatwoMu(A, MU, i)

def curvature_sign(A, i):
	M = fourGammatwo(A,i)
	ev = evs(M)
	if ev[0]<0:
		return -1
	if ev[1] == 0:
		return 0
	else:
		return 1

def curvature_sign_norm(A, i):
	M = fourGammatwoNorm(A,i)
	ev = evs(M)
	if ev[0]<0:
		return -1
	if ev[1] == 0:
		return 0
	else:
		return 1

def fourGamma(A, i):
  B = two_ball(A, i)
  n = len(B)
  RES = [[0 for x in range(n)] for x in range(n)]
  RES[0][0] = 2*sum(B[0])
  for j in range(1,n):
      if B[0][j] == 1:
          RES[j][j] = 2
          RES[0][j] = -2
          RES[j][0] = -2
  return RES

def fourGammaNorm(A, i):
  d = sum(A[i])
  M = (1.0/(1.0*d))*np.array(fourGamma(A,i))
  return M.tolist()

def brange(a, b, n):
    d = (b-a)/float(n)
    res = [a+i*d for i in range(n+1)]
    return res


def  Amat(n,m):
    res = [[0 for i in range(n+m)] for i in range(n*m)]
    for i in range(n):
        for j in range(i*m,(i+1)*m):
            res[j][i] = 1
    for i in range(n,n+m):
        for j in range(n*m):
            if (j % m) == (i-n):
                res[j][i] = 1
    return res

def eta(n,m):
    res = [0 for i in range(n+m)]
    for i in range(1,n):
        res[i] = -1.0/(n-1)
    for i in range(n+1,n+m):
        res[i] = -1.0/(m-1)
    return res

def etap(n,m,p):
    res = [0 for i in range(n+m)]
    res[0] = -p
    res[n] = -p
    for i in range(1,n):
        res[i] = (p-1.0)/(n-1)
    for i in range(n+1,n+m):
        res[i] = (p-1.0)/(m-1)
    return res


def dist(i,j,A):
    if i == j:
        return 0
    if A[i][j] == 1:
        return 1
    for a in range(len(A)):
        if (A[i][a]+A[j][a])==2:
            return 2
    return 3

def d(x,y,A):
    n = sum(A[x])+1
    m = sum(A[y])+1
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
            res.append(dist(i,j,A))
    return res



def curv_calc(A, i):
	M = fourGammatwo(A,i)
	N = fourGamma(A,i)
	ev = evs(M)
	if ev[1] == 0 and ev[0] == 0:
		return 0
	if ev[0] < 0:
		K = 0
		t = 0
		sig = 1
		while sig > 0.00001:
			while t == 0:
				K -= sig
				if evs(np.array(M)-K*np.array(N))[0] == 0:
                        		t += 1
			K+=sig
			t = 0
			sig/=10.0
		return round(K,3)
	if ev[1] > 0:
		K = 0
		t = 0
                sig = 1
                while sig > 0.00001:
                        while t == 0:
                                K += sig
                                if evs(np.array(M)-K*np.array(N))[0] < 0:
                                        t += 1
                        K-=sig
                        t = 0
                        sig/=10.0
                return round(K-(sig*10.0),3)

def curv_calc_norm(A, i):
	M = fourGammatwoNorm(A,i)
	N = fourGammaNorm(A,i)
	ev = evs(M)
	if ev[1] == 0 and ev[0] == 0:
		return 0
	if ev[0] < 0:
		K = 0
		t = 0
		sig = 1
		while sig > 0.00001:
			while t == 0:
				K -= sig
				if evs(np.array(M)-K*np.array(N))[0] == 0:
                        		t += 1
			K+=sig
			t = 0
			sig/=10.0
		return round(K,3)
	if ev[1] > 0:
		K = 0
		t = 0
                sig = 1
                while sig > 0.00001:
                        while t == 0:
                                K += sig
                                if evs(np.array(M)-K*np.array(N))[0] < 0:
                                        t += 1
                        K-=sig
                        t = 0
                        sig/=10.0
                return round(K-(sig*10.0),3)

def dimpart(A, i, n):
	M = two_ball(A,i)
	m = len(M)
	R = [[0 for x in range(m)] for x in range(m)]
	d = sum(A[i])
	R[0][0] = d**2
	for i in range(1, d+1):
		R[0][i] = -d
		R[i][0] = -d
	for i in range(1, d+1):
		for j in range(1,d+1):
			R[i][j] = 1
	RES = (4.0/(1.0*n))*np.array(R)
	return RES

def dimpartNOR(A, i, n):
	M = two_ball(A,i)
	m = len(M)
	R = [[0 for x in range(m)] for x in range(m)]
	d = sum(A[i])
	R[0][0] = d**2
	for i in range(1, d+1):
		R[0][i] = -1.0/(1.0*d)
		R[i][0] = -1.0/(1.0*d)
	for i in range(1, d+1):
		for j in range(1,d+1):
			R[i][j] = 1.0/(1.0*d*d)
	RES = (4.0/(1.0*n))*np.array(R)
	return RES


def dim_curv_calc(A, i, n):
	M = np.array(fourGammatwo(A,i))-dimpart(A, i, n)
	N = fourGamma(A,i)
	ev = evs(M)
	if ev[1] == 0 and ev[0] == 0:
		return 0
	if ev[0] < 0:
		K = 0
		t = 0
		sig = 1
		while sig > 0.00001:
			while t == 0:
				K -= sig
				if evs(np.array(M)-K*np.array(N))[0] >= 0:
                        		t += 1
			K+=sig
			t = 0
			sig/=10.0
		return round(K,3)
	if ev[1] > 0:
		K = 0
		t = 0
                sig = 1
                while sig > 0.00001:
                        while t == 0:
                                K += sig
                                if evs(np.array(M)-K*np.array(N))[0] < 0:
                                        t += 1
                        K-=sig
                        t = 0
                        sig/=10.0
                return round(K-(sig*10.0),3)

def dim_curv_calc_norm(A, i, n):
	M = np.array(fourGammatwoNorm(A,i))-dimpartNOR(A, i, n)
	N = fourGammaNorm(A,i)
	ev = evs(M)
	if ev[1] == 0 and ev[0] == 0:
		return 0
	if ev[0] < 0:
		K = 0
		t = 0
		sig = 1
		while sig > 0.00001:
			while t == 0:
				K -= sig
				if evs(np.array(M)-K*np.array(N))[0] >= 0:
                        		t += 1
			K+=sig
			t = 0
			sig/=10.0
		return round(K,3)
	if ev[1] > 0:
		K = 0
		t = 0
                sig = 1
                while sig > 0.00001:
                        while t == 0:
                                K += sig
                                if evs(np.array(M)-K*np.array(N))[0] < 0:
                                        t += 1
                        K-=sig
                        t = 0
                        sig/=10.0
                return round(K-(sig*10.0),3)

def ocurve(x,y,A):
    dx = sum(A[x])
    dy = sum(A[y])
    return 1+scipy.optimize.OptimizeResult.values(linprog(c = eta(dx+1, dy+1), A_ub = Amat(dx+1,dy+1), b_ub = d(x,y,A), bounds = (None, None)))[3]

def lazocurve(x,y,A,p):
    dx = sum(A[x])
    dy = sum(A[y])
    return 1+scipy.optimize.OptimizeResult.values(linprog(c = etap(dx+1, dy+1, p), A_ub = Amat(dx+1,dy+1), b_ub = d(x,y,A), bounds = (None, None)))[3]



urls = (
  '/', 'index'
)
web.config.debug = False

class index:
    def GET(self):
        user_data = web.input()
  #load data
	try:
		AM = json.loads(user_data.am)
		V  = json.loads(user_data.v)
		t  = json.loads(user_data.t)
	except:
		return user_data.callback+'(["error0"]);'

	#do calcs
	if t == 0:
		try:
			ret = range(len(V));
		except:
			return user_data.callback+'(["error3"]);'
	if t == 1:
		try:
			ret = []
			for j in range(len(V)):
				ret.append(curvature_sign(AM, j))
		except:
			return user_data.callback+'(["error1"]);'
	if t == 2:
		try:
			ret = []
			for j in range(len(V)):
				ret.append(curvature_sign_norm(AM, j))
		except:
			return user_data.callback+'(["error4"]);'
	if t == 3:
		try:
			ret = []
			for j in range(len(V)):
				ret.append(curv_calc(AM, j))
		except:
			return user_data.callback+'(["error5"]);'
	if t == 4:
		try:
			ret = []
			for j in range(len(V)):
				ret.append(curv_calc_norm(AM, j))
		except:
			return user_data.callback+'(["error6"]);'
	if t == 5:
		try:
			dimn  = json.loads(user_data.d)
			if(dimn == 0):
				return user_data.callback+'(["error8"]);'
			elif(dimn < 0):
				return user_data.callback+'(["error8b"]);'
		except:
				return user_data.callback+'(["error9"]);'
		try:
			ret = []
			for j in range(len(V)):
				ret.append(dim_curv_calc(AM, j, dimn))
		except:
			return user_data.callback+'(["error7"]);'
	if t == 6:
		try:
			ret = dict()
			ret["AM"] = AM
			ret["ORC"] = [[0 for i in range(len(V))] for j in range(len(V))]

			for i in range(len(V)):
				for j in range(len(V)):
					if AM[i][j] == 1:
						ret["ORC"][i][j] = ocurve(i,j,AM)
		except:
			return user_data.callback+'(["error10"]);'
	if t == 7:
		try:
			dimn  = json.loads(user_data.d)
			if(dimn == 0):
				return user_data.callback+'(["error8"]);'
			elif(dimn < 0):
				return user_data.callback+'(["error8b"]);'
		except:
				return user_data.callback+'(["error9"]);'
		try:
			ret = []
			for j in range(len(V)):
				ret.append(dim_curv_calc_norm(AM, j, dimn))
		except:
			return user_data.callback+'(["error11"]);'

	#send data back
	try:
		return user_data.callback+"("+json.dumps(ret)+");"
	except:
		return user_data.callback+'(["error2"]);'

if __name__ == "__main__":
    app = web.application(urls, globals())
    app.run()

