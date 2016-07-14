# run `make DEF=...` to add extra defines
OUT = zerndeco
LDFLAGS = -lm -lgsl -lgslcblas -lpng -lcfitsio -fopenmp
SRCS = $(wildcard *.c)
CC = gcc
DEFINES = $(DEF) -D_GNU_SOURCE -DEBUG
CFLAGS = -std=gnu99 -fopenmp -Wall -Werror -Wextra $(DEFINES)
OBJS = $(SRCS:.c=.o)

all : $(OUT)

$(OUT) : $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $(OUT)

$(SRCS) : %.c : %.h $(INDEPENDENT_HEADERS)
	touch $@

%.h: ;

clean:
	/bin/rm -f *.o *~

depend:
	$(CC) -MM $(SRCS)
