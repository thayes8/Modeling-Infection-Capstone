CXX=g++
INCLUDE=/usr/local/include/trng
LIB=trng4
TARGET=modlemc

all: $(TARGET)

$(TARGET): modlemc.cpp
	$(CXX) -o $(TARGET) modlemc.cpp -I$(INCLUDE) -l$(LIB)

######### clean
clean:
	rm -f $(TARGET)