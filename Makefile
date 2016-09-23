PROGRAM_NAME = fastq_filterer

default: build

filter.o: src/filter.c
	gcc -c src/filter.c

build: filter.o
	gcc filter.o -o $(PROGRAM_NAME) -lz

clean:
	rm $(PROGRAM_NAME) filter.o

check:
	bash test/run_tests.sh
