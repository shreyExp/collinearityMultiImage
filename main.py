import numpy as np
import json as js
import sys
from collinearity import collinearity, getRotationMatrix
from triangulation import modelTriangulation, baselineTriangulation
import math as mt

def Main():
	modelNo = int(sys.argv[1])
	doCollinearity(modelNo)

	if modelNo > 1:
		doBaseTriangulation(modelNo)

	mode = 2
	doModelTriangulation(mode, modelNo)


def doModelTriangulation(mode, modelNo):

	cameraVars = getCameraVarsForCollinearity(modelNo)
	ppl = cameraVars["principalPointLeft"]
	ppr = cameraVars["principalPointRight"]
	f = cameraVars["focalLength"]

	suffix = "Im" + str(modelNo) + "Im" + str(modelNo + 1) + ".json"
	fil = open("outputOfCollinearity" + suffix, 'r')
	outputCoplanarity = js.loads(fil.read())
	fil.close()

	mr = getRotationMatrix(outputCoplanarity["anglesRight"])
	mr = np.transpose(mr)
	ml = getRotationMatrix(outputCoplanarity["anglesLeft"])
	ml = np.transpose(ml)

	if mode == 1:
		fil = open("pointsForModelSpaceTriangulationForFurtherBLT" + suffix, "r")
	elif mode == 2:
		fil = open("pointsForModelSpaceTriangulation" + suffix, "r")
	pointsNew = js.loads(fil.read())
	fil.close()

	for point in pointsNew:
		pointLeft = point["left"]
		pointLeft = np.array([pointLeft[0] - ppl[0], pointLeft[1] - ppl[1], -1 * f])
		pointLeft = ml @ pointLeft
		point["left"] = pointLeft
		pointRight = point["right"]
		pointRight = np.array([pointRight[0] - ppr[0], pointRight[1] - ppr[1], -1 * f])
		pointRight = mr @ pointRight
		point["right"] = pointRight
	positionOfCameraOne =  getPositionOfCameraOne(modelNo)
	scaledBaseline = getPositionRightCamera(modelNo)
	scaledBaseline = scaledBaseline - positionOfCameraOne
	#print(pointsNew)
	modelSpace = modelTriangulation(pointsNew, positionOfCameraOne, scaledBaseline)
	#print(modelSpace)
	modeldic = {}
	modeldic["modelpoints"] = modelSpace
	modeldic["rightCameraPos"] = positionOfCameraOne + scaledBaseline
	suffix = "Im" + str(modelNo) + "Im" + str(modelNo + 1) + ".json"
	if mode == 1:
		suffix = "ForFurtherBLT" + suffix
	saveModelOutput(suffix, modeldic)

def doCollinearity(modelNo):
	imageLeftNo = modelNo
	imageRightNo = imageLeftNo + 1
	prefix = "pointsForCollinearity"
	suffix = "Im" + str(imageLeftNo) + "Im" + str(imageRightNo) + ".json"

	f = open(prefix + suffix, 'r')
	points = js.loads(f.read())
	f.close()
	cameraVars = getCameraVarsForCollinearity(modelNo) 

	pixLength = cameraVars["pixLength"]
	scalePoints(points, pixLength)

	staticVars = {}
	staticVars["pointVars"] = points
	staticVars["principalPoint"] = cameraVars["principalPointLeft"]
	staticVars["focalLength"] = cameraVars["focalLength"]

	#if modelNo == 1:
	staticVars["anglesLeft"] = [0,0,0]
	staticVars["cameraPositionLeft"] = np.array([0,0,staticVars["focalLength"]])

	#elif modelNo > 1:
	#	fil = open('outputOfCollinearityIm' + str(modelNo - 1) + "Im" + str(modelNo) + '.json', 'r')
	#	ocoll = js.loads(fil.read())
	#	fil.close()
	#	staticVars["anglesLeft"] = ocoll['anglesRight']
	#	if modelNo == 2:
	#		staticVars["cameraPositionLeft"] = np.array(ocoll['cameraPositionRight'])
	#	elif modelNo > 2:
	#		fil = open('outputOfBaselineTriangulationIm' + str(modelNo -1) + 'Im' + str(modelNo) + '.json', 'r')
	#		obt = js.loads(fil.read())
	#		fil.close()
	#		staticVars["cameraPositionLeft"] = obt['cameraPositionRight']


	angles = staticVars["anglesLeft"]
	staticVars["rotationMatrix"] = getRotationMatrix(angles)
	staticVars["principalPoint"] = pixLength*np.array(staticVars["principalPoint"])
#bEE
	staticVars["baseline"] = 1.0
	noOfPoints = len(staticVars["pointVars"])
	dynamicVars = np.zeros(5 + noOfPoints*3)
	dynamicVars[0:3] = np.array(staticVars["anglesLeft"])
	position = staticVars["cameraPositionLeft"]
	#Initializing position of the right camera with position of Left camera
	dynamicVars[3:5] = position[1:3]
	points = staticVars["pointVars"]
	for i in range(noOfPoints):
		#write here the initial estimates
		dynamicVars[5 + i*3] = points[i]["left"][0]
		dynamicVars[5 + i*3 + 1] = points[i]["left"][1]
		dynamicVars[5 + i*3 + 2] = 0.0
	#end of paste

	outputCollinearity = collinearity(staticVars, dynamicVars)
	outputCollinearity["baseline"] = outputCollinearity["baseline"]/np.linalg.norm(outputCollinearity["baseline"]) 
	baseline = outputCollinearity["baseline"]
	outputCollinearity["anglesRight"] = getAddedAnglesRight(outputCollinearity["anglesRight"], modelNo)
	#print("Printing whole of the output of collinearity dictionary")
	#print(outputCollinearity)
	suffix = "Im" + str(modelNo) + "Im" + str(modelNo + 1) + ".json"
	saveCollinearity(suffix, outputCollinearity)

def getAddedAnglesRight(anglesRight, modelNo):
	fname = "outputOfCollinearityIm" + str(modelNo - 1) + "Im" + str(modelNo) + ".json"
	f = open(fname, 'r')
	oc = js.loads(f.read())
	f.close()
	anglesRight += np.array(oc['anglesRight'])
	return anglesRight


def getPositionRightCamera(modelNo):
	if modelNo == 1:
		suffix = "Im" + str(modelNo) +  "Im" + str(modelNo + 1) + ".json"
		fil = open("outputOfCollinearity" + suffix, 'r')
		oc = js.loads(fil.read())
		fil.close()
		return np.array(oc["baseline"])
	if modelNo > 1:
		suffix = "Im" + str(modelNo) +  "Im" + str(modelNo + 1) + ".json"
		fil = open("outputOfBaselineTriangulation" + suffix, 'r')
		ob = js.loads(fil.read())
		fil.close()
		return np.array(ob["positionRightCamera"])

def scalePoints(points, pixLength):
	for p in points:
		p["left"] = pixLength * np.array(p["left"])
		p["right"] = pixLength * np.array(p["right"])

def doBaseTriangulation(modelNo):

	prefix = "modelspaceForFurtherBLT"
	suffix = "Im" + str(modelNo - 1) + "Im" + str(modelNo) + ".json"
	fil = open(prefix + suffix, "r")
	model = js.loads(fil.read())
	fil.close()
	modelpoints = model["modelpoints"]
	positionOfCameraOne = model["rightCameraPos"]
	model = []
	for mod in modelpoints:
		model.append(np.array(mod))

	fil = open("cameraParameters" + "Im" + str(modelNo + 1) + ".json","r")
	cameraVars = js.loads(fil.read())
	fil.close()
	pp = cameraVars["principalPoint"]
	fl = cameraVars["focalLength"]
	prefix = "outputOfCollinearity"
	suffix = "Im" + str(modelNo) + "Im" + str(modelNo + 1) + ".json"
	fil = open(prefix + suffix, "r")
	cop = js.loads(fil.read())
	fil.close()
	anglesRight = np.array(cop["anglesRight"])
	baseline = np.array(cop["baseline"])
	rmat = getRotationMatrix(anglesRight)
	rmat = np.transpose(rmat)
	prefix = "pointsForBaselineTriangulation"
	fil = open(prefix + suffix, "r")
	points23 = js.loads(fil.read()) 
	fil.close()
	vectors2d = []
	for points in points23:
		vectors2d.append(np.array(points["right"]))
	vectors3d = []
	for vec in vectors2d:
		vec = vec - pp
		vec = np.array([vec[0], vec[1], -1 * fl])
		vec = rmat @ vec
		vectors3d.append(vec)

	#print(model)
	result = baselineTriangulation(positionOfCameraOne, baseline, model, vectors3d)
	#print("RESSSULT")
	#print(result)
	saveBaselineTriangulation(modelNo, result)

def getPositionOfCameraOne(modelNo):
	if modelNo == 1:
		return np.array([0,0,0])
	elif modelNo == 2:
		suffix = "Im" + str(modelNo - 1) + "Im" + str(modelNo) + ".json"
		fil = open("outputOfCollinearity" + suffix, 'r')
		oc = js.loads(fil.read())
		fil.close()
		return np.array(oc["baseline"])
	elif modelNo > 2:
		suffix = "Im" + str(modelNo -1) + "Im" + str(modelNo) + ".json"
		fil = open("outputOfBaselineTriangulation" + suffix, 'r')
		ob = js.loads(fil.read())
		fil.close()
		return np.array(ob["positionRightCamera"])
	

def saveBaselineTriangulation(modelNo, result):
	suffix = "Im" + str(modelNo) + "Im" + str(modelNo + 1) + ".json"
	fil = open("outputOfBaselineTriangulation" + suffix, 'w')
	fil.write(js.dumps(result))
	fil.close()

def saveCollinearity(suffix, outputCollinearity):
	storeVar = {}
	storeVar["baseline"] = list(outputCollinearity["baseline"])
	storeVar["anglesRight"] = list(outputCollinearity["anglesRight"])
	storeVar["anglesLeft"] = list(outputCollinearity["anglesLeft"])
	storeVar["cameraPositionRight"] = list(outputCollinearity["cameraPositionRight"])
	storeVar["cameraPositionLeft"] = list(outputCollinearity["cameraPositionLeft"])
	fil = open("outputOfCollinearity" + suffix, "w")
	fil.write(js.dumps(storeVar))
	fil.close()

def saveModelOutput(suffix, modeldic):
	modeldicOut = {}
	modeldicOut["modelpoints"] = []
	for point in modeldic["modelpoints"]:
		point = list(point)
		modeldicOut["modelpoints"].append(point)
	modeldicOut["rightCameraPos"] = list(modeldic["rightCameraPos"])
	fil = open("modelspace" + suffix, "w")	
	fil.write(js.dumps(modeldicOut))
	fil.close()
	
def getCameraVarsForCollinearity(modelNo):
	cameraVars = {}
	filename = "cameraParametersIm" + str(modelNo) + ".json"	
	f = open(filename, 'r')
	cl = js.loads(f.read())
	f.close()
	cameraVars["principalPointLeft"] = cl["principalPoint"]
	cameraVars["focalLength"] = cl["focalLength"]
	cameraVars["pixLength"] = cl["pixLength"]

	filename = "cameraParametersIm" + str(modelNo + 1) + ".json"	
	f = open(filename, 'r')
	cr = js.loads(f.read())
	f.close()
	#print(cr)
	cameraVars["principalPointRight"] = cr["principalPoint"]
	if modelNo > 1:
		suffix = "Im" + str(modelNo-1) + "Im" + str(modelNo) + ".json"
		filename = "outputOfCollinearity" + suffix
		f = open(filename, 'r')
		oc = js.loads(f.read())
		f.close()
		cameraVars["anglesLeft"] =  oc["anglesRight"]
	else:
		cameraVars["anglesLeft"] = [0,0,0]
	return cameraVars

if __name__ == '__main__':
	Main()
