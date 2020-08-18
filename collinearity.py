import numpy as np
import math as mt
import json as js
import sys
from triangulation import modelTriangulation

def getRotationMatrix(angle):
	omega = angle[0]
	phi = angle[1]
	kappa = angle[2]
	m = np.zeros((3,3))
	m[0,0] = mt.cos(phi) * mt.cos(kappa)
	m[0,1] = mt.sin(omega) * mt.sin(phi) * mt.cos(kappa) + mt.cos(omega)*mt.sin(kappa)
	m[0,2] = -1 * mt.cos(omega) * mt.sin(phi) * mt.cos(kappa) + mt.sin(omega)*mt.sin(kappa)
	m[1,0] = -1 * mt.cos(phi) * mt.sin(kappa)
	m[1,1] = -1 * mt.sin(omega) * mt.sin(phi) * mt.sin(kappa) + mt.cos(omega)*mt.cos(kappa)
	m[1,2] = mt.cos(omega) * mt.sin(phi) * mt.sin(kappa) + mt.sin(omega)*mt.cos(kappa)
	m[2,0] = mt.sin(phi)
	m[2,1] = -1 * mt.sin(omega)*mt.cos(phi)
	m[2,2] = mt.cos(omega)*mt.cos(phi)
	return m


def getRsq(m, rel):
	r = m[0] @ rel
	s = m[1] @ rel
	q = m[2] @ rel
	return(r,s,q)

def blist1(staticVars, dynamicVars, mode, index):
	#point = pointVar["pos"]
	pointOnGround = dynamicVars[(5 + index*3): (5 + index*3 + 3)]
	if mode == 1:
		camera = dynamicVars[3:5]
		xl = np.array([staticVars["baseline"]])
		camera = np.concatenate((xl, camera))
	if mode == 0:
		camera = staticVars["cameraPositionLeft"]
	#camera = dynamicVars["pos"]
	rel = pointOnGround - camera

	if mode == 1:
		m = getRotationMatrix(dynamicVars[0:3])
	if mode == 0:
		m = staticVars["rotationMatrix"]

	f = staticVars["focalLength"]
	r, s, q = getRsq(m, rel)

	if mode == 0:
		angles = staticVars["anglesLeft"]
	if mode == 1:
		angles = dynamicVars[0:3]

	omega = angles[0]
	phi = angles[1]
	kappa = angles[2]
	fq2 = (f/mt.pow(q,2))
	b11 = fq2*(r*(-1*m[2,2]*rel[1] + m[2,1]*rel[2]) - q*(-1*m[0,2]*rel[1] + m[0,1]*rel[2]))

	b12 = fq2*(r*(mt.cos(phi)*rel[0] + mt.sin(omega)*mt.sin(phi)*rel[1] - mt.cos(omega)*mt.sin(phi)*rel[2]))
	b12 = b12  + fq2*(-1*q*(-1*mt.sin(phi)*mt.cos(kappa)*rel[0] + mt.sin(omega)*mt.cos(phi)*mt.cos(kappa)*rel[1]))
	b12 = b12  + fq2*(q*mt.cos(omega)*mt.cos(phi)*mt.cos(kappa)*rel[2])

	b13 = -1*(f/q)*(m[1] @ rel)
	b14 = fq2*(r*m[2,0] - q*m[0,0])

	b15 = fq2*(r*m[2,1] - q*m[0,1])
	b16 = fq2*(r*m[2,2] - q*m[0,2])

	pointVars = staticVars["pointVars"]
	if mode == 0:
		photoPos = pointVars[index]["left"]
	if mode == 1:
		photoPos = pointVars[index]["right"]

	pp = staticVars["principalPoint"]

	relPhoto = photoPos - pp
	J = relPhoto[0] + f*(r/q)
	l = [b11, b12, b13, -1*b14, -1*b15, -1*b16, J]
	return l

def blist2(staticVars, dynamicVars, mode, index):
	#point = pointVar["pos"]
	pointOnGround = dynamicVars[(5 + index*3): (5 + index*3 + 3)]
	if mode == 1:
		camera = dynamicVars[3:5]
		xl = np.array([staticVars["baseline"]])
		camera = np.concatenate((xl, camera))
	if mode == 0:
		camera = staticVars["cameraPositionLeft"]
	#camera = dynamicVars["pos"]
	rel = pointOnGround - camera

	if mode == 1:
		m = getRotationMatrix(dynamicVars[0:3])
	if mode == 0:
		m = staticVars["rotationMatrix"]

	f = staticVars["focalLength"]
	r, s, q = getRsq(m, rel)

	if mode == 0:
		angles = staticVars["anglesLeft"]
	if mode == 1:
		angles = dynamicVars[0:3]
	#point = pointVar["pos"]
	#camera = cameraVars["pos"]
	#rel = point - camera
	#m = cameraVars["rotationMatrix"]
	#f = cameraVars["focalLength"]
	#r, s, q = getRsq(m, rel)
	#angles = cameraVars["angles"]
	omega = angles[0]
	phi = angles[1]
	kappa = angles[2]
	fq2 = (f/mt.pow(q,2))

	b21 = fq2*(s*(-1*m[2,2]*rel[1] + m[2,1]*rel[2]) - q*(-1*m[1,2]*rel[1] + m[1,1]*rel[2]))

	b22 = fq2*(s*(mt.cos(phi)*rel[0] + mt.sin(omega)*mt.sin(phi)*rel[1] - mt.cos(omega)*mt.sin(phi)*rel[2]))
	b22 = b22 + fq2*(-1*q*(mt.sin(phi)*mt.sin(kappa)*rel[0] - mt.sin(omega)*mt.cos(phi)*mt.sin(kappa)*rel[1]))
	b22 = b22 + fq2*(-1*q*(mt.cos(omega)*mt.cos(phi)*mt.sin(kappa)*rel[2]))

	b23 = (f/q)*(m[0] @ rel)
	b24 = fq2 * (s*m[2,0] - q*m[1,0])
	b25 = fq2 * (s*m[2,1] - q*m[1,1])
	b26 = fq2 * (s*m[2,2] - q*m[1,2])

	pointVars = staticVars["pointVars"]
	if mode == 0:
		photoPos = pointVars[index]["left"]
	if mode == 1:
		photoPos = pointVars[index]["right"]

	pp = staticVars["principalPoint"]

	relPhoto = photoPos - pp
	K = relPhoto[1] + f*(s/q)
	
	l = [b21, b22, b23, -1*b24, -1*b25, -1*b26, K]
	return l

def collinearity(staticVars, dynamicVars):
	print("Static Vars")
	print(staticVars)
	print("Initial dynamicVars")
	print(dynamicVars)

	threshold = 0.001
	count = 0

	while(True):
		B = getMatrixBee(dynamicVars, staticVars) 
		count += 1
		E = B[:,(B.shape[1] - 1)]
		B = B[:,0:(B.shape[1] - 1)]
		delta = np.linalg.lstsq(B, E, rcond=None)[0]
		makeChanges(dynamicVars, delta)
		#print("Iteration: " + str(count))
		#print("dynamicVars")
		#print(dynamicVars)
		if np.linalg.norm(delta) < threshold:
			break

	print("FINAL")
	print(dynamicVars)
	output = {}
	output["anglesRight"] = dynamicVars[0:3].tolist()
	yzofCameraRight = staticVars["baseline"]
	cameraPositionRight = np.concatenate((np.array([yzofCameraRight]), dynamicVars[3:5]))
	output["cameraPositionRight"] = cameraPositionRight.tolist()
	output["anglesLeft"] = staticVars["anglesLeft"]
	output["cameraPositionLeft"] = staticVars["cameraPositionLeft"].tolist()

	output["baseline"] = (cameraPositionRight - staticVars["cameraPositionLeft"]).tolist()
	#output["baseline"] = (output["cameraPositionRight"] - output["cameraPositionLeft"]).tolist()
	#output["cameraPositionRight"] = staticVars["cameraPositionLeft"] + baselineVec.tolist()
	#output["cameraPositionRight"] = output["cameraPositionRight"].tolist()
	return output

def triangulateModel():
	fil = open('outputOfRelativeOrientation.json', 'r')
	outputRel = js.loads(fil.read())
	fil.close()

	mr = getRotationMatrix(outputRel["anglesRight"])
	mr = np.transpose(mr)
	ml = getRotationMatrix(outputRel["anglesLeft"])
	ml = np.transpose(ml)

	fil = open('pointsForModelSpaceTriangulationIm1Im2.json', 'r')
	points = js.loads(fil.read())
	fil.close()

	fil = open('cameraParameters.json', 'r')
	cp = js.loads(fil.read())
	fil.close()

	f = cp['focalLength']
	pl = cp['pixLength']
	pp = pl*np.array(cp['principalPoint'])

	#print(pp)
	#print(pp)
	#print(pl)

	for point in points:

		pointLeft = point["left"]
		pointLeft = np.array([pl*pointLeft[0] - pp[0], pl*pointLeft[1] - pp[1], -1*f])
		pointLeft = ml @ pointLeft
		point["left"] = pointLeft
	
		pointRight = point["right"]
		pointRight = np.array([pl*pointRight[0] - pp[0], pl*pointRight[1] - pp[1], -1*f])
		pointRight = ml @ pointRight
		point["right"] = pointRight

	positionOfCameraOne = getPositionOfCameraLeft()
	scaledBaseline = getBaseline()
	scaledBaseline = scaledBaseline/np.linalg.norm(scaledBaseline)
	modelSpace = modelTriangulation(points, positionOfCameraOne, scaledBaseline)
	#print(modelSpace)
	sep = ','
	print('photo' + sep*7 + 'ModelSpace')
	print('left' + sep*3 + 'right')
	for i in range(len(points)):
		st = str(points[i]['left'][0]) +sep + str(points[i]['left'][1]) + sep*2 + str(points[i]['right'][0]) +sep+ str(points[i]['right'][1])
		st = st +  sep*3 + str(modelSpace[i][0]) + sep + str(modelSpace[i][1]) + sep + str(modelSpace[i][2])
		print(st)
	print('omega, phi, kappa')
	angles = outputRel['anglesRight']
	print(str(angles[0]) + sep + str(angles[1]) + sep + str(angles[2]))

def makeChanges(dynamicVars, delta):
	dynamicVars += delta

def getInitialEstimate():
	#Static Vars will contain baseline
	#cameraPositionLeft
	#rotationMatrix
	#focalLength
	#anglesLeft  (left)
	#principalPoint
	f = open("cameraParameters.json", 'r')
	staticVars = js.loads(f.read())
	f.close()

	f = open('Points.json', 'r')
	pointVars = js.loads(f.read())
	f.close()

	pixLength = staticVars["pixLength"]
	n = len(pointVars)
	for i in range(n):
		pointVars[i]["left"] = pixLength*np.array(pointVars[i]["left"])
		pointVars[i]["right"] = pixLength*np.array(pointVars[i]["right"])
	
	print(pointVars)
	staticVars["pointVars"] = pointVars

	staticVars["anglesLeft"] = [0,0,0]
	angles = staticVars["anglesLeft"]
	staticVars["rotationMatrix"] = getRotationMatrix(angles)
	staticVars["cameraPositionLeft"] = np.array(staticVars["position"])
	staticVars["principalPoint"] = pixLength*np.array(staticVars["principalPoint"])
#bEE
	staticVars["baseline"] = 100.0
	noOfPoints = len(staticVars["pointVars"])
	dynamicVars = np.zeros(5 + noOfPoints*3)
	dynamicVars[0:3] = np.array(staticVars["anglesLeft"])
	position = staticVars["cameraPositionLeft"]
	dynamicVars[3:5] = position[1:3]
	points = staticVars["pointVars"]
	for i in range(noOfPoints):
		#write here the initial estimates
		dynamicVars[5 + i*3] = points[i]["left"][0]
		dynamicVars[5 + i*3 + 1] = points[i]["left"][1]
		dynamicVars[5 + i*3 + 2] = 0.0
	
	return (dynamicVars, staticVars)



def getMatrixBee(dynamicVars, staticVars):
	
	P = getMatrixP(dynamicVars, staticVars)
	Q = getMatrixQ(dynamicVars, staticVars)
	R = getMatrixR(dynamicVars, staticVars)
	S = getMatrixS(dynamicVars, staticVars)

	upper = np.concatenate((P, Q), axis=1)
	lower = np.concatenate((R, S), axis=1)

	return np.concatenate((upper, lower))

def getMatrixP(dynamicVars, staticVars):
	noOfPoints = len(staticVars["pointVars"])
	P = np.zeros((noOfPoints*2, 5))
	return P

def getMatrixR(dynamicVars, staticVars):
	noOfPoints = len(staticVars["pointVars"])
	R = np.zeros((noOfPoints*2, 5))
	mode = 1
	for i in range(noOfPoints):
		b1 = blist1(staticVars, dynamicVars, mode, i)
		del b1[3]
		del b1[5]
		b2 = blist2(staticVars, dynamicVars, mode, i)
		del b2[3]
		del b2[5]
		R[i*2] = np.array(b1)
		R[i*2+1] = np.array(b2)
	return R

def getMatrixQ(dynamicVars, staticVars):
	noOfPoints = len(staticVars["pointVars"])
	Q = np.zeros((noOfPoints*2, noOfPoints*3 + 1))
	mode = 0
	for i in range(noOfPoints):
		b1 = blist1(staticVars, dynamicVars, mode, i)
		b1 = -1 * np.array(b1[3:7])
		b2 = blist2(staticVars, dynamicVars, mode, i)
		b2 = -1 * np.array(b2[3:7])
		Q[(i*2), (i*3):(i*3+3)] = b1[0:3]
		Q[(i*2), noOfPoints*3] = -1 * b1[3]
		Q[(i*2 + 1), (i*3):(i*3+3)] = b2[0:3]
		Q[(i*2 + 1), noOfPoints*3] = -1 * b2[3]

	return Q

def getMatrixS(dynamicVars, staticVars):
	noOfPoints = len(staticVars["pointVars"])
	S = np.zeros((noOfPoints*2, noOfPoints*3 + 1))
	mode = 1
	for i in range(noOfPoints):
		b1 = blist1(staticVars, dynamicVars, mode, i)
		b1 = -1 * np.array(b1[3:7])
		b2 = blist2(staticVars, dynamicVars, mode, i)
		b2 = -1 * np.array(b2[3:7])
		S[(i*2), (i*3):(i*3+3)] = b1[0:3]
		S[(i*2), noOfPoints*3] = -1 * b1[3]
		S[(i*2 + 1), (i*3):(i*3+3)] = b2[0:3]
		S[(i*2 + 1), noOfPoints*3] = -1 * b2[3]

	return S

def getPositionOfCameraLeft():
	f = open('cameraParameters.json', 'r')
	cp = js.loads(f.read())
	f.close()
	return np.array(cp['position'])

def getBaseline():
	f = open('outputOfRelativeOrientation.json', 'r')
	ro = js.loads(f.read())
	f.close()
	return np.array(ro['baseline'])

if __name__ == '__main__':
	Main()
