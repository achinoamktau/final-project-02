# Makefile for building spkmeans executable

# Compiler options
CC = gcc
CFLAGS = -Wall -g -Wextra -Werror -pedantic-errors
LDFLAGS = -lm 
# Source files
SOURCES = spkmeans_functions.c spkmeans.c
HEADERS = spkmeans_functions.h

# Object files
OBJECTS = $(SOURCES:.c=.o)

# Build target
TARGET = spkmeans

# Build rules
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $(TARGET) $(LDFLAGS)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@ 

clean:
	rm -f $(TARGET) $(OBJECTS)