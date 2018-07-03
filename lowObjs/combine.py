# usage : python3 combine.py > allObjs

import os.path

objNum = 0
lines = ''
while os.path.exists(str(objNum)+'.obj'):
	objLines = ''
	numLines = 0
	with open(str(objNum)+'.obj', 'r') as f:
		for l in f:
			if l.split(' ')[0] == 'f' or l.split(' ')[0] == 'v':
				objLines += l
				numLines += 1
	lines += str(numLines) + '\n' + objLines
	objNum += 1
lines = str(objNum) + '\n' + lines
print(lines,end='')