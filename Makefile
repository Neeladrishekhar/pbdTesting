
# ./tetgen -pRzf ../upperArmModels/arm/tricep.ply
# tempParticles is tempfile.txt that we get from blender
# Y forward Z up for exporting from blender for objects in lowObjs

# ../apitrace trace ../../../pbdTesting/stressTest < ../../../pbdTesting/interpolatedParticles
# ../apitrace dump-images --call-nos=no stressTest.trace

# triangle no-collision strain PBD :-
# Rendered 855 frames in 1206.16 secs, average of 0.708859 fps

CC := g++
CFLAGS := -Wall -I/usr/include/eigen3
LFLAGS :=
GLFLAGS := -lglut -lGLU -lGL -lm -fPIC

OUT := neighbours objInter collisionReduce choice stressTest
POLYOBJS := $(wildcard lowObjs/*.obj)

.PHONY : all run runChoice runAlpha runCAlpha clean

all: $(OUT) interpolatedParticles

neighbours:Particle.o neighbours.o
	$(CC) $(LFLAGS) $^ -o $@

objInter:Particle.o objInter.o
	$(CC) $(LFLAGS) $^ -o $@

collisionReduce:Particle.o collisionReduce.o
	$(CC) $(LFLAGS) $^ -o $@

choice:Particle.o choice.o
	$(CC) $(LFLAGS) $^ -o $@

stressTest:Particle.o stressTest.o
	$(CC) $(LFLAGS) $^ $(GLFLAGS) -o $@

%.o: %.cpp
	$(CC) -c $< $(CFLAGS) -o $@

# collisionReduce:Particle.o collisionReduce.o
# 	$(CC) $(LFLAGS) $^ -fopenmp -o $@

# choice:Particle.o choice.o
# 	$(CC) $(LFLAGS) $^ -o $@

# stressTest:Particle.o stressTest.o
# 	$(CC) $(LFLAGS) $^ $(GLFLAGS) -fopenmp -o $@

# %.o: %.cpp
# 	$(CC) -c -fopenmp $< $(CFLAGS) -o $@

meta/interpolatedParticlesAlpha: interpolatedParticles meta/collisionReduceOut
	cat $^ > $@

meta/collisionReduceOut: collisionReduce interpolatedParticles
	echo "reduction Start ..." && ./collisionReduce < interpolatedParticles > $@ && echo "... reduction End"

meta/choiceInterParticles: choice interpolatedParticles
	./choice < interpolatedParticles > $@

meta/choiceInterParticlesAlpha: meta/choiceInterParticles meta/collChoiceReduceOut
	cat $^ > $@

meta/collChoiceReduceOut: collisionReduce meta/choiceInterParticles
	echo "reduction Start ..." && ./collisionReduce < meta/choiceInterParticles > $@ && echo "... reduction End"

interpolatedParticles: meta/particles meta/particleInterpolInfo
	cat $^ > $@

meta/particleInterpolInfo: objInter meta/particleIntermediate
	./objInter < meta/particleIntermediate > $@

meta/particleIntermediate: meta/particles lowObjs/allObjs
	cat $^ > $@

meta/particles: meta/tempParticles meta/tempNeighbours
	cat $^ > $@

meta/tempNeighbours: neighbours meta/tempParticles
	./neighbours < meta/tempParticles > $@

lowObjs/allObjs: $(POLYOBJS)
	python3 lowObjs/combine.py > $@

run: stressTest interpolatedParticles
	./stressTest < interpolatedParticles

runChoice: stressTest meta/choiceInterParticles
	./stressTest l < meta/choiceInterParticles

runAlpha: stressTest meta/interpolatedParticlesAlpha
	./stressTest c < meta/interpolatedParticlesAlpha

runCAlpha: stressTest meta/choiceInterParticlesAlpha
	./stressTest c l < meta/choiceInterParticlesAlpha

clean:
	rm -f *.o $(OUT)