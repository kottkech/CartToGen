import numpy
from scipy.spatial import distance

class Atom:
    def __init__(self, index, atom, res, x, y, z):
        self.index = index
        self.atom = atom
        self.res = res
        self.pos = (float(x),float(y),float(z))
        self.left = []
        self.right = []
        self.bl = None
        self.theta = None
        self.phi = None
        self.psi = None
        self.Omega = None

class Backbone:
    def __init__(self):
        self.chain = []
        self.length = 0

    def add(self, atom):
        self.chain += [atom]

        if self.length >= 1:
            if self.length >= 2:
                if self.length >= 3:
                    self.chain[self.length].left += [self.chain[self.length - 3]]
                    self.chain[self.length - 3].right += [self.chain[self.length]]
                self.chain[self.length].left += [self.chain[self.length - 2]]
                self.chain[self.length - 2].right += [self.chain[self.length]]
            self.chain[self.length].left += [self.chain[self.length-1]]
            self.chain[self.length - 1].right += [self.chain[self.length]]

        self.length += 1

    def testPrint(self):
        for atom in self.chain:
            if len(atom.left) > 0:
                print(atom.left[len(atom.left)-1].atom + atom.left[len(atom.left)-1].index + "->", end='')
            print(atom.atom + atom.index, end='')
            if len(atom.right) > 0:
                print("<-" + atom.right[0].atom + atom.right[0].index, end='')
            print(" = " + str(atom.bl))

    def getBl(self):
        for atom in self.chain:
            if len(atom.left) > 0:
                atom.bl = distance.euclidean(atom.pos,atom.left[len(atom.left)-1].pos)


f = open("1qcq.pdb", "r")

bb = Backbone()

for line in f:
    arrIn = line.split()
    if arrIn[0] == "ATOM":
        if arrIn[2] == "N" or arrIn[2] == "CA" or arrIn[2] == "C":
            newAtom = Atom(arrIn[1], arrIn[2], arrIn[3], arrIn[6], arrIn[7], arrIn[8])
            bb.add(newAtom)

f.close()

bb.getBl()
bb.testPrint()



