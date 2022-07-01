import numpy as np  # type: ignore[import]
import copy
from numpy.linalg import eigvalsh  # type: ignore[import]
from numpy import linalg
from enum import Enum

"""
Ricci Flatness

Erin Law 2022
"""

class Flatness(Enum):
    NotRegular = 1
    LargeDegree = 2
    RFlat = 3
    SFlat = 4
    RSFlat = 5
    RandSFlat = 6
    Flat = 7
    NotFlat = 8
	
def flatnessToString( type ):
    if ( type == Flatness.NotRegular ):
        return ""
    if ( type == Flatness.LargeDegree ):
        return ""
    if ( type == Flatness.RFlat ):
        return "R"
    if ( type == Flatness.SFlat ):
        return "S"
    if ( type == Flatness.RSFlat ):
        return "RS"
    if ( type == Flatness.RandSFlat ):
        return "R+S"
    if ( type == Flatness.Flat ):
        return "F"
    if ( type == Flatness.NotFlat ):
        return "X"
    return "ERR"
    
#finds if a vertex is regular by checking if it's neighbours have the same degree
#Note this does not check if the whole graph is regular
def regular(A,x):
    xnbs=[]
    for j in range(len(A)):
        if A[x][j]==1:
            xnbs.append(j)
    vertdeg=[0 for i in range(len(xnbs))]
    vertdeg.append(sum(A[x]))
    for i in range(len(xnbs)):
        vertdeg[i]=sum(A[i])
    if all(x==vertdeg[0] for x in vertdeg):
        return True
    return False

#returns a matrix of possible values for the Ricci flatness matrix
#given an adjacency matrix and vertex
def ChoiceTable(A,x):
    xnbs=[]
    for j in range(len(A)):
        if A[x][j]==1:
            xnbs.append(j)
    D = len(xnbs)
    RFmat= []
    Res = [ [ [] for i in range(D) ] for j in range(D) ]
    for a in range(D):
        for b in range(D):
            vert=[]
            i = xnbs[a]
            n = xnbs[b]
            for j in range(len(A)):
                if A[i][j]==1 and A[n][j]==1:
                    vert.append(j)
            if  vert==[]:
                RFmat.append(0)
            else:
                RFmat.append(vert)
            Res[a][b] = vert
    return Res
 
#Detects whether more than one of the same single value exists in a row
def singleValueChecker(CT):
    N = len(CT)
    for i in range(N):
        singleValues = []
        row = CT[i]
        for j in range(N):
            if (len(row[j]) == 1 ):
                u = row[j][0]
                if u in singleValues:
                    return False
                else:
                    singleValues.append(u)
    return True
 
#If a single value exists in a row then remove it from all other entries in the row
def singleValueRemover(CT):
    N = len(CT)
    for i in range(N):
        singleValues = []
        row = CT[i]
        for j in range(N):
            if (len(row[j]) == 1 ):
                singleValues.append( row[j][0] )
        for s in singleValues:
            for j in range(N):
                if ( len(row[j]) > 1 ):
                    if s in row[j]:
                        row[j].remove(s)
    return 0

#repeats removing single values until all repeats of single values are gone
def singleValueRecursive(CT):
    cont = True
    while cont:
        CTcopy = copy.deepcopy(CT)
        singleValueRemover(CT)
        if not singleValueChecker(CT):
            return False
        CT = CT.transpose()
        singleValueRemover(CT)
        if not singleValueChecker(CT):
            return False
        CT = CT.transpose()
        if CT.all() == CTcopy.all():
            cont = False
    return True

#checks if the grid is entirely single values
def isSolved(CT):
    N = len(CT)
    for i in range(N):
        for j in range(N):
            if ( len(CT[i][j] ) > 1  ):
                return False
    return True
    
def isSSolved(CT):
    N = len(CT)
    for i in range(N):
        for j in range(N):
            if (len(CT[i][j])>1) or CT[i][j]!=CT[j][i]:
                return False
    return True
 
def isRSolved(CT):
	N = len(CT)
	for i in range(N):
		for j in range(N):
			if (len(CT[i][j])>1):
				return False
	for i in range(N):
		if (CT[i][i][0] != 0):
			return False
	return True

def solve(CT):
    #find an entry of length>1
    rIndex, cIndex = 0,0
    N = len(CT)
    minSize = N
    for i in range(N):
        for j in range(N):
            length = len( CT[i][j] )
            if (length>1) and (length < minSize):
            	minSize = length
    ex = False
    for i in range(N):
        for j in range(N):
            if (len(CT[i][j]) == minSize):
                rIndex, cIndex = i,j
                ex = True
                break
        if ex:
            break
    valuesToTry=copy.deepcopy(CT[rIndex][cIndex])
    for value in valuesToTry:
        CTcopy=copy.deepcopy(CT)
        CTcopy[i][j]=[value]
        if not singleValueChecker(CTcopy):
            continue
        if not singleValueRecursive(CTcopy):
            continue
        if isSolved(CTcopy):
            return True
        if solve(CTcopy):
            return True
    return False


def Ssolve(CT):
    #find an entry of length>1
    rIndex, cIndex = 0,0
    N = len(CT)
    minSize = N
    for i in range(N):
        for j in range(N):
            length = len( CT[i][j] )
            if (length>1) and (length < minSize):
            	minSize = length
    ex = False
    for i in range(N):
        for j in range(N):
            if (len(CT[i][j]) == minSize):
                rIndex, cIndex = i,j
                ex = True
                break
        if ex:
            break
    valuesToTry=copy.deepcopy(CT[rIndex][cIndex])
    for value in valuesToTry:
        CTcopy=copy.deepcopy(CT)
        CTcopy[rIndex][cIndex]=[value]
        if value not in CTcopy[cIndex][rIndex]:
            return False
        CTcopy[cIndex][rIndex]=[value]
        if not singleValueChecker(CTcopy):
            continue
        if not singleValueRecursive(CTcopy):
            continue
        if isSolved(CTcopy):
            return True
        if Ssolve(CTcopy):
            return True
    return False

#Check for R-Ricci flatness which means we place zeroes on the diagonals
def RFlat(CT2):
    CT = copy.deepcopy(CT2)
    N = len(CT)
    for i in range(N):
        CT[i][i] = [0]
    if not singleValueRecursive(CT):
        return False
    #At this point we need to start making choices in our graph
    if isSolved(CT):
        return True
    return solve(CT)

def SFlat(CT2):
    CT = copy.deepcopy(CT2)
    if not singleValueRecursive(CT):
        return False
    if isSSolved(CT):
        return True
    return Ssolve(CT)

def RSFlat(CT2):
    CT = copy.deepcopy(CT2)
    N=len(CT)
    for i in range(N):
        for j in range(N):
            CT[i][i]= [0]
    if not singleValueRecursive(CT):
        return False
    if isSSolved(CT):
        return True
    return Ssolve(CT)
 
#checks if a vertex is Ricci flat and if it has R/S flatness
def RicciFlat(A,x):
    if sum(A[x])>10:
        return Flatness.LargeDegree
    if not regular(A,x):
        return Flatness.NotRegular
    choiceTable = np.array(ChoiceTable(A,x),dtype=object)
    if not singleValueRecursive(choiceTable):
        return Flatness.NotFlat
    if isSolved(choiceTable):
        if isSSolved(choiceTable):
        	if isRSolved(choiceTable):
        		return Flatness.RSFlat
        	return Flatness.SFlat
        if isRsolved(choiceTable):
        	return Flatness.RFlat
        return Flatness.Flat
    if RSFlat(choiceTable):
        return Flatness.RSFlat
    if RFlat(choiceTable):
        if SFlat(choiceTable):
            return Flatness.RandSFlat
        return Flatness.RFlat
    if SFlat(choiceTable):
        return Flatness.SFlat
    if solve(choiceTable):
        return Flatness.Flat
    return Flatness.NotFlat
 
def RicciFlatGraph(A):
    vec=[0 for i in range(len(A))]
    for i in range(len(vec)):
        vec[i]=flatnessToString(RicciFlat(A,i))
    return vec

"""
Steinerberger Curvature

Erin Law 2022
"""

def distanceMatrix(A):
    A=np.array(A)
    n = len(A)
    D = copy.deepcopy(A)
    An= copy.deepcopy(A)
    for x in range(n):
        An=A@An
        for i in range(n):
            for j in range(i+1):
                if An[i,j]>0 and D[i,j]==0 and i!=j:
                    D[i,j]=D[j,i]=x+2
    return D

def steinerbergerCurvature(A):
    n = len(A)
    vec = np.array([n for i in range(n)])
    D = distanceMatrix(A)
    isConnected = True
    for i in range(1, n):
        if D[0][i]==0:
            isConnected = False
    Di = linalg.pinv(D)
    curvature = Di@vec
    if not isConnected:
        componentSizes=[1 for i in range(n)]
        for i in range(n):
            for j in range(n):
                if D[i][j] > 0:
                    componentSizes[i] += 1
        for i in range(n):
            curvature[i] = (1.0*curvature[i]* componentSizes[i])/(1.0*n )
    return curvature

"""
Graph curvature calculator functions, written by Ben Snodgrass.

This code is based off formulae given in 'Bakry-Émery curvature on graphs as
an eigenvalue problem' by Cushing et al. (2021).
"""

inf = float('inf')


def non_normalised_unweighted_curvature(A, n):
    """
    Curvature calculator for an arbitrary simple, unweighted graph as a function
    of the adjacency matrix and curvature dimension.

    A is the adjacency matrix of the graph G each vertex has a arbitrarily
    assigned vertex number, determined by position in A.
    """

    q = len(A)
    curvlst = []

    # Switch to Numpy to increase calculation speed
    # In this case, A[i, j] = p_ij as mu[x] = 1 (mu is the measure) for all x, and there is no weighting
    A = np.array(A, dtype=float)

    # list of one-spheres of the vertices
    onesps = [[] for i in range(q)]
    for i in range(q):
        for j in range(q):
            if A[i, j] == 1:
                onesps[i].append(j)
    # number of nearest neighbours of each vertex, in this case also the degree of each vertex
    lenonesps = [len(onesps[i]) for i in range(q)]
    # twosp is a matrix whose (i, j)th element = 0 iff i is in the two ball of j and vice versa
    # if i and j have a common neighbour (and so could be in each other's two balls), A[i].A[j] must be positive
    A_2 = np.matmul(A, A)
    twosp = copy.copy(A_2)
    for i in range(q):
        for j in range(q):
            # a point is not in its own two-ball
            if i == j:
                twosp[i, j] = 0
            # also not in two ball if the vertices are adjacent
            if A[i, j] == 1:
                twosp[i, j] = 0
            # more than one common neighbour is irrelevant to being in two-sphere or not
            if twosp[i, j] != 0:
                twosp[i, j] = 1
            # create the following matrices to perform summations. Use np.matmul rather than list
            # comprehension sums as they are far quicker this performs first summation in eq. A.11
    sum2 = np.matmul(A, twosp)
    # create matrix of reciprocals of p^(2)_ij with terms corresponding to vertices not in two_sphere(x)
    # replaced by zero. This corresponds to only summing over z in S2(x)
    recipA_2 = twosp*A_2
    for i in range(q):
        for j in range(q):
            if recipA_2[i, j] != 0:
                recipA_2[i, j] = 1/recipA_2[i, j]
    # this performs the last summation in eq. A.11
    sum3 = np.matmul(recipA_2, A)
    # (i, j)th element of first matrix is the list [p[i, z]*p[j, z] for z in range(q)],
    # which, dotted with recipA_2[x], gives the summation in A.12
    sum4 = np.matmul([[[A[i, z]*A[j, z] for z in range(q)] for j in range(q)] for i in range(q)], recipA_2)
    for x in range(q):
        m = lenonesps[x]
        # isolated points have default curvature 0
        if m == 0:
            curvlst.append(0)
        else:
            onesp = onesps[x]
            # in the non-normalised, non-weighted case, (v_0)(v_0)^T is simplly a matrix of 1's
            A_n = [[-2/n for i in range(m)] for j in range(m)]
            for i in range(m):
                # formula A.11 can be dramatically simplified by noting that p_xy = 1 for all y~x
                # both p-terms in eq. A.13 are 1 and can be ignored
                # d_x, the degree, is just m in this case
                A_n[i][i] += 5/2-1/2*m+2*A_2[x, onesp[i]]+3/2*sum2[onesp[i], x]-2*sum3[x, onesp[i]]
            for i in range(m):
                for j in range(m):
                    if i != j:
                        # again, p_xy = 1 for y~x simplifies eq. A.2
                        A_n[i][j] += 1-2*A[onesp[i], onesp[j]]-2*sum4[onesp[i], onesp[j], x]
            # eigvalsh returns list of eigenvalues of A_n, smallest first, so eigvalsh(A_n)[0] is the
            # smallest eigenvalue of A_n, which is the curvature
            curvlst.append(round(eigvalsh(A_n)[0], 3))
    # returns list of curvature values at all vertices, ordered by vertex index
    return curvlst


def normalised_unweighted_curvature(A, n):
    """
    Normalised Curvature calculator for an arbitrary simple, unweighted graph as a function
    of the adjacency matrix and curvature dimension.

    A is the adjacency matrix of the graph G. Each vertex has a arbitrarily
    assigned vertex number, determined by position in A.
    """
    #Switch to Numpy to increase calculation speed
    A = np.array(A, dtype = float)
    q = len(A)
    #list of one-spheres of the vertices
    onesps = [[] for i in range(q)]
    for i in range(q):
        for j in range(q):
            if A[i, j] == 1:
                onesps[i].append(j)
    #number of nearest neighbours of each vertex
    lenonesps = [len(onesps[i]) for i in range(q)]
    curvlst = []
    #p_xy values for all combinations of vertices x, y in G
    r = lenonesps[0]
    #if graph is r-regular, calculation can be simplified
    regular = True
    z = 1
    while z < q and regular == True:
        if lenonesps[z] != r:
            regular = False
        z += 1
    #isolated points have 0 curvature
    if regular == True:
        if r == 0:
            return [0 for i in range(q)]
        #In non-normalised case, the curvature is 1/r*(curvature in non-normalised case), so the
        #calculation can be performed exactly the same as non-normalised case, with 1/r factor added at end
        #Code is copied directly from non-normalised case, with annotations mostly unchanged. It treats it
        #as non-normalised curvature, so p_xy = A_xy in this section
        else:
            #p_xy = 1/r*A_xy, so A^2 = A_2 = r**2*p_2
            A_2 = np.matmul(A, A)
            twosp = copy.copy(A_2)
            for i in range(q):
                for j in range(q):
                    #a point is not in its own two-ball
                    if i == j:
                        twosp[i, j] = 0
                    #also not in two ball if the vertices are adjacent
                    if A[i, j] == 1:
                        twosp[i, j] = 0
                    #more than one common neighbour is irrelevant to being in two-sphere or not
                    if twosp[i, j] != 0:
                        twosp[i, j] = 1
            #create the following matrices to perform summations. Use np.matmul rather than list
            #comprehension sums as it is far quicker
            #this performs second summation in eq. A.11
            sum1 = 3/2*np.matmul(A, twosp)
            #create matrix of reciprocals of p^(2)_ij with terms corresponding to vertices not
            #in two_sphere(x) removed. The zero entries do not figure in the calculation anyway
            #so the values are irrelevant
            recipA_2 = twosp*A_2
            for i in range(q):
                for j in range(q):
                    if recipA_2[i, j] != 0:
                        recipA_2[i, j] = 1/recipA_2[i, j]
            #this performs the last summation in eq. A.11
            sum3 = 2*np.matmul(recipA_2, A)
            #(i, j)th element of first matrix is the list [A[i, z]*A[j, z] for z in range(q)],
            #which, dotted with recipA_2[x], gives the summation in A.12
            sum4 = 2*np.matmul([[[A[i, z]*A[j, z] for z in range(q)] for j in range(q)] for i in range(q)], recipA_2)
            for x in range(q):
                onesp = onesps[x]
                #in the non-normalised, non-weighted case, (v_0)(v_0)^T is simplly a matrix of 1's
                A_n = [[-2/n for i in range(r)] for j in range(r)]
                for i in range(r):
                    #formula A.11 can be dramatically simplified by noting that p_xy = 1 for all y~x
                    #both p-terms in eq. A.13 are 1 and can be ignored
                    #d_x, the degree, is r
                    A_n[i][i] += 5/2-1/2*r+2*A_2[x, onesp[i]]+sum1[onesp[i], x]-sum3[x, onesp[i]]
                for i in range(r):
                    for j in range(r):
                        if i != j:
                            #again, p_xy = 1 for y~x simplifies eq. A.12
                            A_n[i][j] += 1-2*A[onesp[i], onesp[j]]-sum4[onesp[i], onesp[j], x]
                #eigvalsh returns list of eigenvalues of A_n, smallest first, so eigvalsh(A_n)[0]
                #is the smallest eigenvalue of A_n, which is the curvature
                #1/r adds required normalisation
                curvlst.append(1/r*eigvalsh(A_n)[0])
        #returns list of curvature values at all vertices, ordered by vertex index
        return list(np.around(curvlst, 3))

    else:
        #create matrix of reciprocals of p^(2)_ij. If p^(2)_ij = 0, the will not figure in the
        #calculation anyway so the value is irrelevant
        P = A*[[1/lenonesps[i]] for i in range(q)]
        #matrix p_2 = p^2 has the (x, z)th element identical to p^(2)_xz
        P_2 = np.matmul(P, P)
        #twosp is a matrix whose (i, j)th element = 0 iff i is in the two ball of j and vice versa
        #if i and j have a common neighbour (and so could be in each other's two balls), p[i].p[j]
        #must be positive, hence use of p_2
        twosp = copy.copy(P_2)
        for i in range(q):
            for j in range(q):
                #a point is not in its own two-ball
                if i == j:
                    twosp[i, j] = 0
                #also not in two ball if the vertices are adjacent
                if A[i, j] == 1:
                    twosp[i, j] = 0
                #more than one common neighbour is irrelevant to being in two-sphere or not
                if twosp[i, j] != 0:
                    twosp[i, j] = 1
        recipP_2 = twosp*P_2
        for i in range(q):
            for j in range(q):
                if recipP_2[i, j] != 0:
                    recipP_2[i, j] = 1/recipP_2[i, j]
        #create the following matrices to perform summations. Use np.matmul rather than list
        #comprehension sums as they are far quicker
        #this performs first summation in eq. A.11
        sum1 = 3/2*np.matmul(P, twosp)
        #used in first part of second summation in A.11. Multiplication by A 'filters out'
        #terms not corresponding to vertices in the one-sphere of a point
        sum2 = 3/2*np.matmul(P, A)
        #the performs last summation in eq. A.11
        sum3 = 2*np.matmul(recipP_2, np.transpose(np.square(P)))
        #(i, j)th element of this is the list [p[i, z]*p[j, z] for z in range(q)], which is
        #used in conjunction with recipp_2 in line ** to perform summation in eq. A.12
        sum4 = 2*np.matmul([[[P[i, z]*P[j, z] for z in range(q)] for j in range(q)] for i in range(q)], np.transpose(recipP_2))
        for x in range(q):
            m = lenonesps[x]
            #isolated points have default curvature 0
            if m == 0:
                curvlst.append(0)
            else:
                onesp = onesps[x]
                #create -2/n*(v_0)(v_0)^T component of A_n
                A_n = [[-2/n*(P[x, onesp[i]]*P[x, onesp[j]])**0.5 for i in range(m)] for j in range(m)]
                #add A_infinity terms, using the above 'sum' matrices
                for i in range(m):
                    #since d_x = mu[x] (d_x is vertex degree, mu is the vertex weighting
                    # in this case, d_x/mu[x] in eq. A.11 is simply 1
                    #on-diagonal terms use A.11 and A.13
                    #divide through by p(x, onesp[i]) as in formula A.13
                    try:
                        A_n[i][i] += P[x, onesp[i]]-1/2+3/2*P[onesp[i], x]+sum1[onesp[i], x]+sum2[onesp[i], x]+1/2*P_2[x, onesp[i]]/P[x, onesp[i]]-P[x, onesp[i]]*sum3[x, onesp[i]]
                    except:
                        return onesp[i]
                for i in range(m):
                    for j in range(m):
                        if i != j:
                            #off-diagonal terms use equations A.12 and A.13
                            A_n[i][j] += (P[x, onesp[i]]*P[x, onesp[j]]-P[x, onesp[i]]*P[onesp[i], onesp[j]]-P[x, onesp[j]]*P[onesp[j], onesp[i]]-P[x, onesp[i]]*P[x, onesp[j]]*sum4[onesp[i], onesp[j], x])/(P[x, onesp[i]]*P[x, onesp[j]])**0.5
                #eigvalsh returns list of eigenvalues of A_n, smallest first, so
                #eigvalsh(A_n)[0] is the smallest eigenvalue of A_n, which is the curvature
                curvlst.append(eigvalsh(A_n)[0])
        #returns list of curvature values at all vertices, ordered by vertex index
        return list(np.around(curvlst, 3))
