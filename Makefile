CXX=g++
PGCC=pgc++
INCLUDE=/usr/local/include/trng
LIB=trng4
TARGET=modlemc
AOPTS= pgc++ -fast -ta=tesla:cc75 -Minfo=accel
CC_ACC_M= -fast -ta=tesla:cc75,managed -Minfo=accel -lcurand -Mcuda




all: $(TARGET) GPUmodle fillRand.o

$(TARGET): modlemc.cpp
	$(CXX) -o $(TARGET) modlemc.cpp -I$(INCLUDE) -l$(LIB)

fillRand.o: fillRand.cu
	nvcc -c fillRand.cu -o fillRand.o

GPUmodle: TommyGPU.cpp fillRand.cu
	$(PGCC) ${CC_ACC_M} -o GPUmodle TommyGPU.cpp

######### clean
clean:
	rm -f $(TARGET) GPUmodle fillRand.o