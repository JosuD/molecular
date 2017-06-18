# .:: Makefile ::.

# Definiciones comunes.
CC:=gcc

# En principio, se pueden compilar dos versiones: una de `debug'
# (sin optimizaciones y con símbolos para `gdb') y otra optimizada.
#
# Nota: puede usarse una macro (dentro del código) para indicar código
# que deba compilarse sólo cuando se busca generar la versión de `debug'.
# Estoy tomando este nombre como `_DEBUG'.

TARGET:=molecular
TARGET_DBG:=$(addsuffix _dbg,$(TARGET))

OBJS:=$(patsubst %.c,%.o,$(wildcard *.c))
OBJS_DBG:=$(patsubst %.c,%_dbg.o,$(wildcard *.c))

CFLAGS:=-std=c99 -I=/usr/include/gsl -O3 -Wall -finput-charset=UTF-8
CFLAGS_DBG:=-std=c99 -I=/usr/include/gsl -D_DEBUG -O0 -Wall -ggdb -finput-charset=UTF-8

# Sería genial que nos dejaran usar la GSL...
LIBS:=-lm -lgsl -lgslcblas
HEADERS:=$(wildcard *.h)

.PHONY: clean

all: $(TARGET) $(TARGET_DBG)

# Optimizada.
$(TARGET): $(OBJS)
	$(CC) $^ $(LIBS) -o $@

$(OBJS): %.o : %.c $(HEADERS)
	$(CC) -c $(CFLAGS) $< -o $@

# Debug.
$(TARGET_DBG): $(OBJS_DBG)
	$(CC) $^ $(LIBS) -o $@

$(OBJS_DBG): %_dbg.o : %.c $(HEADERS)
	$(CC) -c $(CFLAGS_DBG) $< -o $@

# Clean
clean:
	rm -f $(TARGET) $(TARGET_DBG) $(OBJS) $(OBJS_DBG) *~ *.*~
