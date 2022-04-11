TARGETS = modldemc1d.c

all: $(TARGETS)

modldemc1d: modldemc1d.c
	gcc -o modldemc1d modldemc1d.c

clean:
	rm -f $(TARGETS)