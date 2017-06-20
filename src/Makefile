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

TARGET_PRUEBAS:=pruebas
TARGET_PRUEBAS_DBG:=$(addsuffix _dbg,$(TARGET_PRUEBAS))

OBJS:=$(patsubst %.c,%.o,$(wildcard *.c))
OBJS_DBG:=$(patsubst %.c,%_dbg.o,$(wildcard *.c))

CFLAGS:=-std=c99 -I=/usr/include/gsl -O3 -Wall -finput-charset=UTF-8
CFLAGS_DBG:=-std=c99 -I=/usr/include/gsl -D_DEBUG -O0 -Wall -ggdb -finput-charset=UTF-8

LIBS:=-lm -lgsl -lgslcblas
HEADERS:=$(wildcard *.h)

# Agregar otra receta si se añade otro archivo de pruebas con un main() definido.

.PHONY: clean

all: $(TARGET) $(TARGET_DBG) $(TARGET_PRUEBAS) $(TARGET_PRUEBAS_DBG)

# Targets.
$(TARGET): $(filter-out pruebas.o, $(OBJS))
	$(CC) $^ $(LIBS) -o $@

$(TARGET_DBG): $(filter-out pruebas_dbg.o, $(OBJS_DBG))
	$(CC) $^ $(LIBS) -o $@

$(TARGET_PRUEBAS): $(filter-out main.o, $(OBJS))
	$(CC) $^ $(LIBS) -o $@

$(TARGET_PRUEBAS_DBG): $(filter-out main_dbg.o, $(OBJS_DBG))
	$(CC) $^ $(LIBS) -o $@

# Optimizada.
$(OBJS): %.o : %.c $(HEADERS)
	$(CC) -c $(CFLAGS) $< -o $@

# Debug.
$(OBJS_DBG): %_dbg.o : %.c $(HEADERS)
	$(CC) -c $(CFLAGS_DBG) $< -o $@

# Clean
clean:
	rm -f $(TARGET_PRUEBAS) $(TARGET_PRUEBAS_DBG) $(TARGET) $(TARGET_DBG) $(OBJS) $(OBJS_DBG) *~ *.*~
