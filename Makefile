CC = gcc
CFLAGS = -Wall -Wextra -O2
TARGET = diagonalization

all:
	$(CC) $(CFLAGS) src/main.c -o $(TARGET) -lm

clean:
	rm -f $(TARGET)
