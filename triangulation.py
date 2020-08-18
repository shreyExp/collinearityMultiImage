import numpy as np
import math as mt
import json as js
import sys
#from coplanarity2 import coplanarity

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


def solve(point, positionOfCameraOne, baseline):
	a1 = point["left"]
	a2 = point["right"]

	baselineAdd = positionOfCameraOne
	A = np.zeros((3,3))
	A[0] = a1
	A[1] = -1 * a2
	A[2] = baseline

	A = np.transpose(A)
	B = A[:,2]
	A = A[:,0:2]
	scales = np.linalg.lstsq(A, B, rcond=None)[0]

	return (scales[0] * a1 + baselineAdd, baseline + scales[1] * a2 + baselineAdd)

	

def modelTriangulation(points, positionOfCameraOne, baseline):
	objSpaceList = []
	for point in points:
		objSpace = solve(point, positionOfCameraOne, baseline)
		objSpaceList.append(objSpace[0])
	return objSpaceList



def baselineTriangulation(positionOfCameraOne, baseline, modelpoints, vectors):

	#fl = cameraVars["focalLength"]
	#ppl = np.array(cameraVars["principalPointLeft"])
	#ppr = np.array(cameraVars["principalPointRight"])


	#ml = cameraVars["rotationMatrixLeft"]
	#ml = np.transpose(ml)
	#mr = cameraVars["rotationMatrixRight"]
	#mr = np.transpose(mr)

	#b = cameraVars["baseline"]
	#baselineAdd = cameraVars["baselineAdd"]

	n = len(modelpoints)
	s = (n*3, n+2)
	m = np.zeros(s)
	#
	#rotatedVectors = []

	#f = open("scales.json", 'r')
	#previousScale = js.loads(f.read())
	#f.close()
	i = 0
	for point in modelpoints:
		#pointLeft = point["left"]
		#pointRight = point["right"]
		#r1 = np.array([pointLeft[0] - ppl[0], pointLeft[1] - ppl[1], -1*fl])
		#r2 = np.array([pointRight[0] - ppr[0], pointRight[1] - ppr[1], -1*fl])

		#a1 = previousScale[i]*ml @ r1
		#a2 = mr @ r2
		#rotatedVectors.append(a2)
		##a1 += cameraVars["baselineAdd"] 
		##a2 += cameraVars["baselineAdd"] 

		m[i*3:(i+1)*3, i+1] = vectors[i]
		m[i*3:(i+1)*3, 0] = baseline
		m[i*3:(i+1)*3, n+1] = point - positionOfCameraOne
		i += 1

	A = m[:,:n+1]
	B = m[:,n+1]
	scalesOut = np.linalg.lstsq(A,B, rcond=None)
	residual = scalesOut[1]
	scales = scalesOut[0]
	print("scales")
	print(scales)
	print("baseline")
	print(scales[0]*baseline)
	sep = ','
	outputPoints = []
	for i in range(n):
		position = scales[i+1] * vectors[i]
		position += scales[0]*baseline + positionOfCameraOne
		outputPoints.append(list(position))
		print(position)
		#stri = str(pointCoplan[i]["left"][0]) + sep + str(pointCoplan[i]["left"][1])+ sep*3
		#stri += str(pointCoplan[i]["right"][0])+ sep + str(pointCoplan[i]["right"][1])+ sep*5
		#stri += str(position[0])+ sep + str(position[1])+ sep + str(position[2])
		#
		#print(stri)
	print("Scaled baseline:")
	print(scales[0]*baseline)
	print("scale for baseline: ")
	print(scales[0])
	print("scales for other points in order: ")
	print(scales[1:])
	print("residual")
	print(residual)
	result = {}
	result["positionRightCamera"] = list(positionOfCameraOne + scales[0]*baseline)
	result["points"] = outputPoints
	return result

	


if __name__ == '__main__':
	Main()
