#!/usr/bin/python3

import numpy as np
import json as js

f = open('modelspaceIm1Im2.json', 'r')
d = js.loads(f.read())
f.close()
model1 = d["modelpoints"]

f = open('modelspaceIm2Im3.json', 'r')
d = js.loads(f.read())
f.close()
model2 = d["modelpoints"]
modelnp1 = []
for m in model1:
	modelnp1.append(np.array(m))


modelnp2 = []
for m in model2:
	modelnp2.append(np.array(m))

n = len(modelnp1)
distance = []
for i in range(n):
	distance.append(np.linalg.norm(modelnp1[i] - modelnp2[i]))

s = ","
print("ModelOne" + s*4 + "ModelTwo" + s*4 + "Norm(Distance)")
for i in range(n):
	str1 = str(modelnp1[i][0]) + s + str(modelnp1[i][1]) + s + str(modelnp1[i][2])
	str2 = str(modelnp2[i][0]) + s + str(modelnp2[i][1]) + s + str(modelnp2[i][2])
	str3 = str(distance[i])
	print(str1+ s*2 + str2 + s*2 + str3)
