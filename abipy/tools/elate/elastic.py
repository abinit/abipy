# -*- coding: utf-8 -*-

import json
import math

import numpy as np
from scipy import optimize

__author__ = "Romain Gaillac and François-Xavier Coudert"
__version__ = "2025.08.26"
__license__ = "MIT"


def dirVec(theta, phi):
    """Return a unit vector associated with angles theta and phi"""
    return [math.sin(theta)*math.cos(phi), math.sin(theta)*math.sin(phi), math.cos(theta)]


def dirVec1(theta, phi, chi):
    """Return the first unit vector associated with angles (theta, phi, chi)"""
    return [math.sin(theta)*math.cos(phi), math.sin(theta)*math.sin(phi), math.cos(theta)]


def dirVec2(theta, phi, chi):
    """Return the second unit vector associated with angles (theta, phi, chi)"""
    return [math.cos(theta)*math.cos(phi)*math.cos(chi) - math.sin(phi)*math.sin(chi),
            math.cos(theta)*math.sin(phi)*math.cos(chi) + math.cos(phi)*math.sin(chi),
            - math.sin(theta)*math.cos(chi)]


# Functions to minimize/maximize
def minimize(func, dim):
    """Find the global minimum of a function with 2 (theta, phi) or 3 (theta, phi, chi) parameters"""

    if dim == 2:
        # Mechanical properties are centrosymmetric so we reduce the range for phi
        r = ((0, np.pi), (0, np.pi))
        n = 25
    elif dim == 3:
        # Mechanical properties are centrosymmetric so we reduce the range for phi and chi
        r = ((0, np.pi), (0, np.pi), (0, np.pi))
        n = 10

    # With good heuristics for our number of points, brute-force is not too bad
    # because we are in low dimension.
    return optimize.brute(func, r, Ns=n, full_output=True, finish=optimize.fmin)[0:2]


def maximize(func, dim):
    """Find the global maximum of a function with 2 (theta, phi) or 3 (theta, phi, chi) parameters"""

    res = minimize(lambda x: -func(x), dim)
    return (res[0], -res[1])


class Elastic:
    """
    An elastic tensor, along with methods to access it.

    It will raise a ValueError when the input is invalid, and a TypeError
    if the input is valid but corresponds to a 2D material.
    """

    def __init__(self, s):
        """Initialize the elastic tensor from a string"""

        # Argument can be a 6-line string, a list of list, or a string representation of the list of list

        # If the argument is JSON, decode it
        if isinstance(s, str):
            try:
                s = json.loads(s)
            except:
                pass

        if isinstance(s, np.ndarray) or isinstance(s, list):
            # If we already have a list or numpy array, fine
            mat = s
        elif isinstance(s, str):
            # Remove braces and pipes
            s = s.replace("|", " ").replace("(", " ").replace(")", " ")

            # Remove empty lines
            lines = [line for line in s.split('\n') if line.strip()]
            if len(lines) == 3:
                raise TypeError("is a 2D material")
            if len(lines) != 6:
                raise ValueError("should have three or six rows")

            # Convert to list of list of floats
            try:
                mat = [list(map(float, line.split())) for line in lines]
            except:
                raise ValueError("not all entries are numbers")
        else:
            raise ValueError("invalid argument as matrix")

        # Make it into a square matrix
        try:
            mat = np.array(mat)
        except:
            # Is it upper triangular?
            if list(map(len, mat)) == [6,5,4,3,2,1]:
                mat = [[0]*i + mat[i] for i in range(6)]
                mat = np.array(mat)
            if list(map(len, mat)) == [3,2,1]:
                mat = [[0]*i + mat[i] for i in range(3)]
                mat = np.array(mat)

            # Is it lower triangular?
            if list(map(len, mat)) == [1,2,3,4,5,6]:
                mat = [mat[i] + [0]*(5-i) for i in range(6)]
                mat = np.array(mat)
            if list(map(len, mat)) == [1,2,3]:
                mat = [mat[i] + [0]*(2-i) for i in range(3)]
                mat = np.array(mat)

        if not isinstance(mat, np.ndarray):
            raise ValueError("should be a square or triangular matrix")
        if mat.shape == (3,3):
            raise TypeError("is a 2D material")
        if mat.shape != (6,6):
            raise ValueError("should be a square or triangular matrix")

        # Check that is is symmetric, or make it symmetric
        if np.linalg.norm(np.tril(mat, -1)) == 0:
            mat = mat + np.triu(mat, 1).transpose()
        if np.linalg.norm(np.triu(mat, 1)) == 0:
            mat = mat + np.tril(mat, -1).transpose()
        if np.linalg.norm(mat - mat.transpose()) > 1e-3:
            raise ValueError("should be symmetric, or triangular")
        elif np.linalg.norm(mat - mat.transpose()) > 0:
            # It was almost symmetric: symmetrize it completely
            mat = 0.5 * (mat + mat.transpose())

        # Store it
        self.CVoigt = mat

        # Put it in a more useful representation
        try:
            self.SVoigt = np.linalg.inv(self.CVoigt)
        except:
            raise ValueError("matrix is singular")

        VoigtMat = [[0, 5, 4], [5, 1, 3], [4, 3, 2]]
        def SVoigtCoeff(p,q):
            return 1. / ((1+p//3)*(1+q//3))

        self.Smat = [[[[ SVoigtCoeff(VoigtMat[i][j], VoigtMat[k][l]) * self.SVoigt[VoigtMat[i][j]][VoigtMat[k][l]]
                         for i in range(3) ] for j in range(3) ] for k in range(3) ] for l in range(3) ]
        return

    def is2D(self):
        return False

    def isOrthorhombic(self):
        def iszero(x): return (abs(x) < 1.e-3)
        return (iszero(self.CVoigt[0][3]) and iszero(self.CVoigt[0][4]) and iszero(self.CVoigt[0][5])
                and iszero(self.CVoigt[1][3]) and iszero(self.CVoigt[1][4]) and iszero(self.CVoigt[1][5])
                and iszero(self.CVoigt[2][3]) and iszero(self.CVoigt[2][4]) and iszero(self.CVoigt[2][5])
                and iszero(self.CVoigt[3][4]) and iszero(self.CVoigt[3][5]) and iszero(self.CVoigt[4][5]))

    def isCubic(self):
        def iszero(x): return (abs(x) < 1.e-3)
        return (iszero(self.CVoigt[0][3]) and iszero(self.CVoigt[0][4]) and iszero(self.CVoigt[0][5])
                and iszero(self.CVoigt[1][3]) and iszero(self.CVoigt[1][4]) and iszero(self.CVoigt[1][5])
                and iszero(self.CVoigt[2][3]) and iszero(self.CVoigt[2][4]) and iszero(self.CVoigt[2][5])
                and iszero(self.CVoigt[3][4]) and iszero(self.CVoigt[3][5]) and iszero(self.CVoigt[4][5])
                and iszero(self.CVoigt[0][0] - self.CVoigt[1][1]) and iszero(self.CVoigt[0][0] - self.CVoigt[2][2])
                and iszero(self.CVoigt[0][0] - self.CVoigt[1][1]) and iszero(self.CVoigt[0][0] - self.CVoigt[2][2])
                and iszero(self.CVoigt[3][3] - self.CVoigt[4][4]) and iszero(self.CVoigt[3][3] - self.CVoigt[5][5])
                and iszero(self.CVoigt[0][1] - self.CVoigt[0][2]) and iszero(self.CVoigt[0][1] - self.CVoigt[1][2]))

    def Young(self, x):
        a = dirVec(x[0], x[1])
        r = sum([ a[i]*a[j]*a[k]*a[l] * self.Smat[i][j][k][l]
                  for i in range(3) for j in range(3) for k in range(3) for l in range(3) ])
        return 1/r

    def Young_2(self, x, y):
        a = dirVec(x, y)
        r = sum([ a[i]*a[j]*a[k]*a[l] * self.Smat[i][j][k][l]
                  for i in range(3) for j in range(3) for k in range(3) for l in range(3) ])
        return 1/r

    def LC(self, x):
        a = dirVec(x[0], x[1])
        r = sum([ a[i]*a[j] * self.Smat[i][j][k][k]
                  for i in range(3) for j in range(3) for k in range(3) ])
        return 1000 * r

    def LC_2(self, x, y):
        a = dirVec(x, y)
        r = sum([ a[i]*a[j] * self.Smat[i][j][k][k]
                  for i in range(3) for j in range(3) for k in range(3) ])
        return 1000 * r

    def shear(self, x):
        a = dirVec(x[0], x[1])
        b = dirVec2(x[0], x[1], x[2])
        r = sum([ a[i]*b[j]*a[k]*b[l] * self.Smat[i][j][k][l]
                  for i in range(3) for j in range(3) for k in range(3) for l in range(3) ])
        return 1/(4*r)

    def Poisson(self, x):
        a = dirVec(x[0], x[1])
        b = dirVec2(x[0], x[1], x[2])
        r1 = sum([ a[i]*a[j]*b[k]*b[l] * self.Smat[i][j][k][l]
                   for i in range(3) for j in range(3) for k in range(3) for l in range(3) ])
        r2 = sum([ a[i]*a[j]*a[k]*a[l] * self.Smat[i][j][k][l]
                   for i in range(3) for j in range(3) for k in range(3) for l in range(3) ])
        return -r1/r2

    def averages(self):
        A = (self.CVoigt[0][0] + self.CVoigt[1][1] + self.CVoigt[2][2]) / 3
        B = (self.CVoigt[1][2] + self.CVoigt[0][2] + self.CVoigt[0][1]) / 3
        C = (self.CVoigt[3][3] + self.CVoigt[4][4] + self.CVoigt[5][5]) / 3
        a = (self.SVoigt[0][0] + self.SVoigt[1][1] + self.SVoigt[2][2]) / 3
        b = (self.SVoigt[1][2] + self.SVoigt[0][2] + self.SVoigt[0][1]) / 3
        c = (self.SVoigt[3][3] + self.SVoigt[4][4] + self.SVoigt[5][5]) / 3

        KV = (A + 2*B) / 3
        GV = (A - B + 3*C) / 5

        KR = 1 / (3*a + 6*b)
        GR = 5 / (4*a - 4*b + 3*c)

        KH = (KV + KR) / 2
        GH = (GV + GR) / 2

        return [ [KV, 1/(1/(3*GV) + 1/(9*KV)), GV, (1 - 3*GV/(3*KV+GV))/2],
                 [KR, 1/(1/(3*GR) + 1/(9*KR)), GR, (1 - 3*GR/(3*KR+GR))/2],
                 [KH, 1/(1/(3*GH) + 1/(9*KH)), GH, (1 - 3*GH/(3*KH+GH))/2] ]

    def eigenvalues(self):
        return np.sort(np.linalg.eig(self.CVoigt)[0])

    def shear2D(self, x):
        ftol = 0.001
        xtol = 0.01
        def func1(z): return self.shear([x[0], x[1], z[0]])
        r1 = optimize.minimize(func1, np.pi/2.0, args=(), method = 'Powell', options={"xtol":xtol, "ftol":ftol})
        def func2(z): return -self.shear([x[0], x[1], z[0]])
        r2 = optimize.minimize(func2, np.pi/2.0, args=(), method = 'Powell', options={"xtol":xtol, "ftol":ftol})
        return (float(r1.fun), -float(r2.fun))

    def shear3D(self, x, y, guess1 = np.pi/2.0, guess2 = np.pi/2.0):
        tol = 0.0005
        def func1(z): return self.shear([x, y, z[0]])
        r1 = optimize.minimize(func1, guess1, args=(), method = 'L-BFGS-B')
        def func2(z): return -self.shear([x, y, z[0]])
        r2 = optimize.minimize(func2, guess2, args=(), method = 'L-BFGS-B')
        return (float(r1.fun), -float(r2.fun), float(r1.x[0]), float(r2.x[0]))

    def Poisson2D(self, x):
        ftol = 0.001
        xtol = 0.01
        def func1(z): return self.Poisson([x[0], x[1], z[0]])
        r1 = optimize.minimize(func1, np.pi/2.0, args=(), method = 'Powell', options={"xtol":xtol, "ftol":ftol})
        def func2(z): return -self.Poisson([x[0], x[1], z[0]])
        r2 = optimize.minimize(func2, np.pi/2.0, args=(), method = 'Powell', options={"xtol":xtol, "ftol":ftol})
        return (min(0,float(r1.fun)), max(0,float(r1.fun)), -float(r2.fun))

    def Poisson3D(self, x, y, guess1 = np.pi/2.0, guess2 = np.pi/2.0):
        tol = 0.005
        def func1(z): return self.Poisson([x, y, z[0]])
        r1 = optimize.minimize(func1, guess1, args=(), method = 'L-BFGS-B')
        def func2(z): return -self.Poisson([x, y, z[0]])
        r2 = optimize.minimize(func2, guess2, args=(), method = 'L-BFGS-B')
        return (min(0,float(r1.fun)), max(0,float(r1.fun)), -float(r2.fun), float(r1.x[0]), float(r2.x[0]))


class ElasticOrtho(Elastic):
    """An elastic tensor, for the specific case of an orthorhombic system"""

    def __init__(self, arg):
        """Initialize from a matrix, or from an Elastic object"""
        if isinstance(arg, str):
            Elastic.__init__(self, arg)
        elif isinstance(arg, Elastic):
            self.CVoigt = arg.CVoigt
            self.SVoigt = arg.SVoigt
            self.Smat = arg.Smat
        else:
            raise TypeError("ElasticOrtho constructor argument should be string or Elastic object")

    def Young(self, x):
        ct2 = math.cos(x[0])**2
        st2 = 1 - ct2
        cf2 = math.cos(x[1])**2
        sf2 = 1 - cf2
        s11 = self.Smat[0][0][0][0]
        s22 = self.Smat[1][1][1][1]
        s33 = self.Smat[2][2][2][2]
        s44 = 4 * self.Smat[1][2][1][2]
        s55 = 4 * self.Smat[0][2][0][2]
        s66 = 4 * self.Smat[0][1][0][1]
        s12 = self.Smat[0][0][1][1]
        s13 = self.Smat[0][0][2][2]
        s23 = self.Smat[1][1][2][2]
        return 1/(ct2**2*s33 + 2*cf2*ct2*s13*st2 + cf2*ct2*s55*st2 + 2*ct2*s23*sf2*st2 + ct2*s44*sf2*st2 + cf2**2*s11*st2**2 + 2*cf2*s12*sf2*st2**2 + cf2*s66*sf2*st2**2 + s22*sf2**2*st2**2)

    def LC(self, x):
        ct2 = math.cos(x[0])**2
        cf2 = math.cos(x[1])**2
        s11 = self.Smat[0][0][0][0]
        s22 = self.Smat[1][1][1][1]
        s33 = self.Smat[2][2][2][2]
        s12 = self.Smat[0][0][1][1]
        s13 = self.Smat[0][0][2][2]
        s23 = self.Smat[1][1][2][2]
        return 1000 * (ct2 * (s13 + s23 + s33) + (cf2 * (s11 + s12 + s13) + (s12 + s22 + s23) * (1 - cf2)) * (1 - ct2))

    def shear(self, x):
        ct = math.cos(x[0])
        ct2 = ct*ct
        st2 = 1 - ct2
        cf = math.cos(x[1])
        sf = math.sin(x[1])
        sf2 = sf*sf
        cx = math.cos(x[2])
        cx2 = cx*cx
        sx = math.sin(x[2])
        sx2 = 1 - cx2
        s11 = self.Smat[0][0][0][0]
        s22 = self.Smat[1][1][1][1]
        s33 = self.Smat[2][2][2][2]
        s44 = 4 * self.Smat[1][2][1][2]
        s55 = 4 * self.Smat[0][2][0][2]
        s66 = 4 * self.Smat[0][1][0][1]
        s12 = self.Smat[0][0][1][1]
        s13 = self.Smat[0][0][2][2]
        s23 = self.Smat[1][1][2][2]
        r = (
            ct2*ct2*cx2*s44*sf2 + cx2*s44*sf2*st2*st2 + 4*cf**3*ct*cx*(-2*s11 + 2*s12 + s66)*sf*st2*sx
            + 2*cf*ct*cx*sf*(ct2*(s44 - s55) + (4*s13 - 4*s23 - s44 + s55 - 4*s12*sf2 + 4*s22*sf2 - 2*s66*sf2)*st2)*sx
            + s66*sf2*sf2*st2*sx2 + cf**4*st2*(4*ct2*cx2*s11 + s66*sx2)
            + ct2*(2*cx2*(2*s33 + sf2*(-4*s23 - s44 + 2*s22*sf2))*st2 + s55*sf2*sx2)
            + cf**2*(ct2*ct2*cx2*s55 + ct2*(-2*cx2*(4*s13 + s55 - 2*(2*s12 + s66)*sf2)*st2 + s44*sx2)
                     + st2*(cx2*s55*st2 + 2*(2*s11 - 4*s12 + 2*s22 - s66)*sf2*sx2))
            )
        return 1/r

    def Poisson(self, x):
        ct = math.cos(x[0])
        ct2 = ct*ct
        st2 = 1 - ct2
        cf = math.cos(x[1])
        sf = math.sin(x[1])
        cx = math.cos(x[2])
        sx = math.sin(x[2])
        s11 = self.Smat[0][0][0][0]
        s22 = self.Smat[1][1][1][1]
        s33 = self.Smat[2][2][2][2]
        s44 = 4 * self.Smat[1][2][1][2]
        s55 = 4 * self.Smat[0][2][0][2]
        s66 = 4 * self.Smat[0][1][0][1]
        s12 = self.Smat[0][0][1][1]
        s13 = self.Smat[0][0][2][2]
        s23 = self.Smat[1][1][2][2]

        return (
    (-(ct**2*cx**2*s33*st2) - cf**2*cx**2*s13*st2*st2 - cx**2*s23*sf**2*st2*st2 + ct*cx*s44*sf*st2*(ct*cx*sf + cf*sx) -
            ct**2*s23*(ct*cx*sf + cf*sx)**2 - cf**2*s12*st2*(ct*cx*sf + cf*sx)**2 - s22*sf**2*st2*(ct*cx*sf + cf*sx)**2 +
            cf*ct*cx*s55*st2*(cf*ct*cx - sf*sx) - cf*s66*sf*st2*(ct*cx*sf + cf*sx)*(cf*ct*cx - sf*sx) -
            ct**2*s13*(cf*ct*cx - sf*sx)**2 - cf**2*s11*st2*(cf*ct*cx - sf*sx)**2 - s12*sf**2*st2*(cf*ct*cx - sf*sx)**2)/
            (ct**4*s33 + 2*cf**2*ct**2*s13*st2 + cf**2*ct**2*s55*st2 + 2*ct**2*s23*sf**2*st2 + ct**2*s44*sf**2*st2 +
            cf**4*s11*st2*st2 + 2*cf**2*s12*sf**2*st2*st2 + cf**2*s66*sf**2*st2*st2 + s22*sf**4*st2*st2)
        )


class Elastic2D:
    """Elastic tensor for 2D material, along with methods to access it"""

    def __init__(self, s):
        """Initialize the elastic tensor from a string"""

        # Argument can be a 3-line string, a list of list, or a string representation of the list of list
        try:
            if isinstance(json.loads(s), list):
                s = json.loads(s)
        except:
            pass

        if isinstance(s, str):
            # Remove braces and pipes
            s = s.replace("|", " ").replace("(", " ").replace(")", " ")

            # Remove empty lines
            lines = [line for line in s.split('\n') if line.strip()]
            if len(lines) != 3:
                raise ValueError("should have three rows")

            # Convert to float
            try:
                mat = [list(map(float, line.split())) for line in lines]
            except:
                raise ValueError("not all entries are numbers")
        elif isinstance(s, list):
            # If we already have a list, simply use it
            mat = s
        else:
            raise ValueError("invalid argument as matrix")

        # Make it into a square matrix
        try:
            mat = np.array(mat)
        except:
            # Is it upper triangular?
            if list(map(len, mat)) == [3,2,1]:
                mat = [ [0]*i + mat[i] for i in range(3) ]
                mat = np.array(mat)

            # Is it lower triangular?
            if list(map(len, mat)) == [1,2,3]:
                mat = [ mat[i] + [0]*(2-i) for i in range(3) ]
                mat = np.array(mat)

        if not isinstance(mat, np.ndarray):
            raise ValueError("should be a square or triangular matrix")
        if mat.shape != (3,3):
            raise ValueError("should be a square or triangular matrix")

        # Check that is is symmetric, or make it symmetric
        if np.linalg.norm(np.tril(mat, -1)) == 0:
            mat = mat + np.triu(mat, 1).transpose()
        if np.linalg.norm(np.triu(mat, 1)) == 0:
            mat = mat + np.tril(mat, -1).transpose()
        if np.linalg.norm(mat - mat.transpose()) > 1e-3:
            raise ValueError("should be symmetric, or triangular")
        elif np.linalg.norm(mat - mat.transpose()) > 0:
            mat = 0.5 * (mat + mat.transpose())

        # Store it
        self.CVoigt = mat

        # Put it in a more useful representation
        try:
            self.SVoigt = np.linalg.inv(self.CVoigt)
        except:
            raise ValueError("matrix is singular")

        # Store the important components
        self.s11 = self.SVoigt[0][0]
        self.s12 = self.SVoigt[0][1]
        self.s16 = self.SVoigt[0][2]
        self.s22 = self.SVoigt[1][1]
        self.s26 = self.SVoigt[1][2]
        self.s66 = self.SVoigt[2][2]
        return

    def is2D(self):
        return True

    def Young(self, theta):
        ct = math.cos(theta)
        st = math.sin(theta)

        return 1/(self.s11*ct**4 + self.s22*st**4 + 2*self.s16*ct**3*st + 2*self.s26*ct*st**3 + (2*self.s12+self.s66)*ct**2*st**2)

    def shear(self, theta):
        ct = math.cos(theta)
        st = math.sin(theta)

        calc = ((self.s11 + self.s22 - 2*self.s12)*ct**2*st**2
                + self.s66/4*(ct**4 + st**4 - 2*st**2*ct**2)
                + self.s16*(st**3*ct - ct**3*st)
                + self.s26*(ct**3*st - st**3*ct))
        return 1 / (4 * calc)

    def Poisson(self, theta):
        ct = math.cos(theta)
        st = math.sin(theta)

        num = ((self.s11 + self.s22 - self.s66)*ct**2*st**2
               + self.s12*(ct**4 + st**4)
               + self.s16*(st**3*ct - ct**3*st)
               + self.s26*(ct**3*st - st**3*ct))
        denom = ((2*self.s12 + self.s66)*ct**2*st**2
                 + self.s11*ct**4 + self.s22*st**4
                 + 2*self.s16*st**3*ct
                 + 2*self.s26*ct**3*st)
        return -num/denom

    def eigenvalues(self):
        return np.sort(np.linalg.eig(self.CVoigt)[0])
