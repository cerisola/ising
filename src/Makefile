.PHONY: default all clean distclean temperature_sweep outputdir

CC = gcc

LD = $(CC)

CFLAGS = -std=c99 -Wall -pedantic -O3
LDFLAGS = -lm

SRCS_TSWEEP = ising_temperature_sweep.c
OBJS_TSWEEP = $(subst .c,.o,$(SRCS_TSWEEP))

SRCS_COMMON = $(filter-out $(SRCS_TSWEEP), $(wildcard *.c))
OBJS_COMMON = $(subst .c,.o,$(SRCS_COMMON))

OUTDIR = ../bin
TARGET_TSWEEP = $(OUTDIR)/ising_temperature_sweep

default: all

all: temperature_sweep

temperature_sweep: $(OBJS_COMMON) $(OBJS_TSWEEP) outputdir
	$(CC) -o $(TARGET_TSWEEP) $(OBJS_COMMON) $(OBJS_TSWEEP) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $^ -o $@

outputdir:
	@mkdir -p $(OUTDIR)

clean:
	$(RM) $(OBJS_COMMON) $(OBJS_TSWEEP)

distclean: clean
	$(RM) $(TARGET_TSWEEP)
