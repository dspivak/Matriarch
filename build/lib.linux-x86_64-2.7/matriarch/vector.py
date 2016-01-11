#VECTOR OPERATIONS

import math

#class to let us apply a vector function to a vector without having to change references
class triple(list):
    def alter(self, transformation):
        self.extend(transformation(self))
        self.pop(0)
        self.pop(0)
        self.pop(0)
        return self
    def evaluate(self, transformation):
        return triple(transformation(self))
    def change(self, newvalue):
        self.extend(newvalue)
        self.pop(0)
        self.pop(0)
        self.pop(0)

#basically just an exp(i*theta)
def angleToVector(angle):
    return [math.cos(angle), math.sin(angle), 0]

#inverse of angleToVector
def vectorToAngle(vector):
    output = math.acos(vector[0])
    if vector[1] < 0:
        return 2*math.pi - output
    else:
        return output

def translate(vector):
    return lambda x: [x[0]+vector[0],x[1]+vector[1],x[2]+vector[2]]

def minus(vec):
    output = []
    for i in range(0,len(vec)):
        output.append(-vec[i])
    return output

#returns a rotation about an axis
def rotateAxis(axis, angle, fixedpoint=[0,0,0]):
    def rotate(point, axis, angle):
        [a,b,c] = unit(axis)
        [u,v,w] = normalvector([a,b,c])
        M = [[a*a+(-c*v+b*w)*(math.sin(angle)*u+math.cos(angle)*(-c*v+b*w))+u*(math.cos(angle)*u-math.sin(angle)*(-c*v+b*w)),a*b+(c*u-a*w)*(math.sin(angle)*u+math.cos(angle)*(-c*v+b*w))+v*(math.cos(angle)*u-math.sin(angle)*(-c*v+b*w)),a*c+(-b*u+a*v)*(math.sin(angle)*u+math.cos(angle)*(-c*v+b*w))+w*(math.cos(angle)*u-math.sin(angle)*(-c*v+b*w))],[a*b+(-c*v+b*w)*(math.sin(angle)*v+math.cos(angle)*(c*u-a*w))+u*(math.cos(angle)*v-math.sin(angle)*(c*u-a*w)),b*b+(c*u-a*w)*(math.sin(angle)*v+math.cos(angle)*(c*u-a*w))+v*(math.cos(angle)*v-math.sin(angle)*(c*u-a*w)),b*c+(-b*u+a*v)*(math.sin(angle)*v+math.cos(angle)*(c*u-a*w))+w*(math.cos(angle)*v-math.sin(angle)*(c*u-a*w))],[a*c+u*(-math.sin(angle)*(-b*u+a*v)+math.cos(angle)*w)+(-c*v+b*w)*(math.cos(angle)*(-b*u+a*v)+math.sin(angle)*w),b*c+v*(-math.sin(angle)*(-b*u+a*v)+math.cos(angle)*w)+(c*u-a*w)*(math.cos(angle)*(-b*u+a*v)+math.sin(angle)*w),c*c+w*(-math.sin(angle)*(-b*u+a*v)+math.cos(angle)*w)+(-b*u+a*v)*(math.cos(angle)*(-b*u+a*v)+math.sin(angle)*w)]]
        output = [0,0,0]
        for i in range(0,3):
            for j in range(0,3):
                output[i]+= M[i][j]*point[j]
        return output
    return lambda x: translate(fixedpoint)(rotate(translate(minus(fixedpoint))(x), axis, angle))

#returns composite of rotations
def rotatecomp(axis, normal, thetaNormal, thetaAnormal, torsion):
    #1 radian rotation thetaNormal, thetaAnormal, torsion topologically generate the group SO(3)
    #thetaNormal is the rotation about normal
    #thetaAnormal is rotation about translation x normal
    #torsion is rotation about translation
    #this returns the image of point after the prescribed rotations then translating (if translate = True)
    #composition order is torsional angle then thetaAnormal then thetaNormal then translate (if translate = True)
    #Torsion
    def aux(vector):
        output = rotateaxis(axis, torsion)(vector)
        #Anormal
        output = rotateaxis(crossproduct(axis, normal), thetaAnormal)(output)
        #Normal
        output = rotateaxis(normal, thetaNormal)(output)
        return output
    return aux

#moves the fixed point of a Euclidean automorphism
def changefixedpoint(transformation, fixedpoint = [0,0,0]):
    return lambda x: translate(fixedpoint)(transformation(translate(minus(fixedpoint))(x)))
    
def innerproduct(vec1, vec2):
    output = 0
    for i in range(0,3):
        output += vec1[i] * vec2[i]
    return output
    
def orthogonalComponent(vec, thingToWhichOutputIsOrthogonal):
    innerProductValue = innerproduct(vec, thingToWhichOutputIsOrthogonal)
    lengthSquared = norm(thingToWhichOutputIsOrthogonal) * norm(thingToWhichOutputIsOrthogonal)
    return translate(scale(thingToWhichOutputIsOrthogonal,-innerProductValue/lengthSquared))(vec)

def orthogonal(source1, source2, target1, target2):
    #given two pairs of orthogonal vectors, constructs the unique element of SO_3 sending sourcei to a positive multiple of targeti for i=1,2
    source = [unit(source1),unit(source2)]
    target = [unit(target1),unit(target2)]
    source.append(crossproduct(source[0],source[1]))
    target.append(crossproduct(target[0],target[1]))
    M = []
    for i in range(0,3):
        temp = []
        for j in range(0,3):
            value = 0
            for k in range(0,3):
                value+= source[k][j]*target[k][i]
            temp.append(value)
        M.append(temp)
    def transformation(x):
        output = []
        for i in range(0,3):
            value = 0
            for j in range(0,3):
                value += M[i][j]*x[j]
            output.append(value)
        return output
    return transformation

#figures out where the next nitrogen is supposed to be in an amino acid
def Nlocation(newC, newCA, newO):
    CtoCA = [1.246, 0.427, 0.763]
    CtoO = [0.011, -0.068, -1.23]
    CtoN = [-1.076, -0.297, 0.728]
    tminusC = translate(minus(newC))
    newCtoCA, newCtoO = tminusC(newCA), tminusC(newO)
    transformation = orthogonal(CtoCA, crossproduct(CtoCA,CtoO), newCtoCA, crossproduct(newCtoCA,newCtoO))
    return translate(newC)(transformation(CtoN))

def normalvector(vector):
    #returns a normal to non-zero vector
    if vector[0] != 0:
        u = (-vector[1]-vector[2])/vector[0]
        output = []
        output.append(u/math.sqrt(u*u+2))
        for i in range(1,3):
            output.append(1/math.sqrt(u*u+2))
    else:
        temp = []
        for i in range(0,3):
            temp.append(vector[(i+1)%3])
        temp2 = normalvector(temp)
        output = []
        for i in range(0,3):
            output.append(temp2[(i-1)%3])
    return output

def crossproduct(vec1, vec2):
    output = []
    for i in range(0,3):
        output.append(vec1[(i+1)%3]*vec2[(i+2)%3]-vec1[(i+2)%3]*vec2[(i+1)%3])
    return output

def norm(vec):
    normsquared = 0
    for i in range(0,3):
        normsquared += vec[i]*vec[i]
    return math.sqrt(normsquared)
    
def scale(vec, scalar):
    output = []
    for i in range(0,3):
        output.append(vec[i] * scalar)
    return output

def unit(vec):
    #returns the unit vector that is a positive multiple of vec
    output = []
    length = norm(vec)
    for i in range(0,3):
        output.append(vec[i]/length)
    return output

#projects vector onto a unit vector
def project(vector, onto):
    #returns the projection of "vector" onto "onto"
    temp = unit(onto)
    dotproduct = 0
    for i in range(0,3):
        dotproduct += vector[i]*temp[i]
    output = []
    for i in range(0,3):
        output.append(dotproduct*temp[i])
    return output

#computes the unique rigid motion sending point to point 2
#whose derivative sends v to v2 and w to w2
def T(point,v,w,point2,v2,w2):
    vectortransform = orthogonal(v,w,v2,w2)
    return lambda x: translate(point2)(vectortransform(translate(minus(point))(x)))

#returns a pair of a rigid motion and its derivative as in the T(...)
def Tlist(list1,list2):
    return [T(list1[0],list1[1],list1[2],list2[0],list2[1],list2[2]), orthogonal(list1[1],list1[2],list2[1],list2[2])]
