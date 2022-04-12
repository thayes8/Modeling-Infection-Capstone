
all: modlemc1d

modldemc1d: modldemc1d.c infect.c
	gcc -o modldemc1d modldemc1d.c infect.c

clean:
	rm -f modlemc1d