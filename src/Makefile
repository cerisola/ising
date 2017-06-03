.PHONY: default all clean distclean temperature_sweep outputdir

CC = gcc

LD = $(CC)

CFLAGS = -std=c99 -Wall -pedantic -O3
LDFLAGS = -lm

SRCS_TSWEEP = ising_temperature_sweep.c
OBJS_TSWEEP = $(subst .c,.o,$(SRCS_TSWEEP))

SRCS_AUTOCORR = ising_autocorrelation.c
OBJS_AUTOCORR = $(subst .c,.o,$(SRCS_AUTOCORR))

SRCS_CONTDATA = ising_continuous_data.c
OBJS_CONTDATA = $(subst .c,.o,$(SRCS_CONTDATA))

SRCS_COMMON = $(filter-out $(SRCS_TSWEEP) $(SRCS_AUTOCORR) $(SRCS_CONTDATA), $(wildcard *.c))
OBJS_COMMON = $(subst .c,.o,$(SRCS_COMMON))

OUTDIR = ../bin
TARGET_TSWEEP = $(OUTDIR)/ising_temperature_sweep
TARGET_AUTOCORR = $(OUTDIR)/ising_autocorrelations
TARGET_CONTDATA = $(OUTDIR)/ising_continuous_data

default: all

all: temperature_sweep autocorrelations continuous_data

temperature_sweep: $(OBJS_COMMON) $(OBJS_TSWEEP) outputdir
	$(CC) -o $(TARGET_TSWEEP) $(OBJS_COMMON) $(OBJS_TSWEEP) $(LDFLAGS)

autocorrelations: $(OBJS_COMMON) $(OBJS_AUTOCORR) outputdir
	$(CC) -o $(TARGET_AUTOCORR) $(OBJS_COMMON) $(OBJS_AUTOCORR) $(LDFLAGS)

continuous_data: $(OBJS_COMMON) $(OBJS_CONTDATA) outputdir
	$(CC) -o $(TARGET_CONTDATA) $(OBJS_COMMON) $(OBJS_CONTDATA) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $^ -o $@

outputdir:
	@mkdir -p $(OUTDIR)

clean:
	$(RM) $(OBJS_COMMON) $(OBJS_TSWEEP) $(OBJS_AUTOCORR) $(OBJS_CONTDATA)

distclean: clean
	$(RM) $(TARGET_TSWEEP) $(TARGET_AUTOCORR) $(TARGET_CONTDATA)
