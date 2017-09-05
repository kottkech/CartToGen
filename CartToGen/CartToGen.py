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
        self.dihedral = None

class Backbone:
    def __init__(self):
        self.chain = []
        self.length = 0

    def add(self, atom):
        self.chain += [atom]
	
        if self.length >= 1:
            self.chain[self.length].left += [self.chain[self.length-1]]
            if self.length >= 2:
                self.chain[self.length].left += [self.chain[self.length - 2]]
                if self.length >= 3:
                    self.chain[self.length].left += [self.chain[self.length - 3]]
                    self.chain[self.length - 3].right += [self.chain[self.length]]
                self.chain[self.length - 2].right += [self.chain[self.length]]
            self.chain[self.length - 1].right += [self.chain[self.length]]

        self.length += 1

    def getBl(self):
        for atom in self.chain:
            if len(atom.left) > 0:
                atom.bl = distance.euclidean(atom.pos,atom.left[0].pos)

    def getAngles(self):
        for atom in self.chain:
            if len(atom.left) > 0:
		if len(atom.right) > 0:
                    b1 = numpy.subtract(atom.pos, atom.left[0].pos)
                    b2 = numpy.subtract(atom.right[0].pos, atom.pos)
                    atom.theta = numpy.degrees(numpy.arccos(numpy.dot(b1,b2) / (numpy.linalg.norm(b1) * numpy.linalg.norm(b2))))

    def getDihedrals(self):
        for atom in self.chain:
            if len(atom.left) > 1:
                if len(atom.right) > 0:
                    b1 = numpy.subtract(atom.left[0].pos, atom.left[1].pos)
                    b2 = numpy.subtract(atom.pos, atom.left[0].pos)
                    b3 = numpy.subtract(atom.right[0].pos, atom.pos)

                    atom.dihedral = numpy.sign(numpy.dot(numpy.cross(b1,b2),b3)) * numpy.degrees(
                        numpy.arccos(
                            numpy.dot(numpy.cross(b1,b2),numpy.cross(b2,b3)) / (numpy.linalg.norm(numpy.cross(b1,b2))*numpy.linalg.norm(numpy.cross(b2,b3)))
                        )
                    )

    def printToFile(self):
        blt = open("blAndTheta.txt", "r+")
        blt.truncate()
        blt.write("%-8s %-5s %-15s %-10s\n" % ("Index", "Atom", "Bond_Length", "Theta"))
        dh = open("dihedrals.txt", "r+")
        dh.truncate()
        dh.write("%-8s %-8s %-8s %-8s %-8s\n" % ("Index", "Residue", "Phi", "Psi", "Omega"))
        count=0
        first=True
        for atom in self.chain:
            if atom.theta is not None:
                blt.write("%-8s %-5s %-15.3f %-10.3f\n" % (atom.index, atom.atom, atom.bl, atom.theta))
            else:
                blt.write("%-8s %-5s %-15s %-10s\n" % (atom.index, atom.atom, "N/A", "N/A"))
                
            if not first:
                if atom.atom == "CA":
                    count += 1
                    dh.write("\n%-8s %-8s " % (count, atom.res)) 
                if atom.dihedral is not None:
                    dh.write("%-8.3f " % (atom.dihedral))
                else:
                    dh.write("%-8s " % ("N/A"))
            first=False

        blt.close()
        dh.close()

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
bb.getAngles()
bb.getDihedrals()
bb.printToFile()

