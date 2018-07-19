# usage : python3 combine.py > allObjs

import os.path

objNum = 5
lines = ''
boneList = ['clavicle','scapula','humerus','radius','ulna']

for bone in boneList:
	numVertices = 0
	verticies = ''
	numFaces = 0
	faces = ''
	with open(bone+'.obj', 'r') as f:
		for l in f:
			if l.split(' ')[0] == 'v':
				verticies += l
				numVertices += 1
			if l.split(' ')[0] == 'f':
				faces += l
				numFaces += 1
	lines += str(numVertices) + '\n' + verticies + str(numFaces) + '\n' + faces

lines = str(objNum) + '\n' + lines
print(lines,end='')