CXX=g++
PGCC=pgc++
INCLUDE=/usr/local/include/trng
LIB=trng4
TARGET=modlemc
AOPTS= -fast -ta=tesla:cc75 -Minfo=accel



all: $(TARGET) GPUmodle

$(TARGET): modlemc.cpp
	$(CXX) -o $(TARGET) modlemc.cpp -I$(INCLUDE) -l$(LIB)

GPUmodle: TommyGPU.cpp
	$(PGCC) ${AOPTS} -o GPUmodle TommyGPU.cpp -I$(INCLUDE) -l$(LIB)

######### clean
clean:
	rm -f $(TARGET) GPUmodle