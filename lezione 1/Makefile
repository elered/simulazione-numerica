CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe

es1 : esercizio1.o random.o
	$(CC) random.o esercizio1.o -o es1

es1giusto : giusto.o random.o
	$(CC) random.o giusto.o -o es1giusto

esercizio1.2 : esercizio1.2.o random.o
	$(CC) random.o esercizio1.2.o -o esercizio1.2

esercizio1.3 : esercizio1.3.o random.o
	$(CC) random.o esercizio1.3.o -o esercizio1.3

main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

giusto.o : giusto.cpp random.h
	$(CC) -c giusto.cpp -o giusto.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out
