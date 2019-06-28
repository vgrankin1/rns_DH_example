NAME=rnsdh
OBJ=main.o

LIBS=-lgmp
#-lpthread





CFLAGS=-g
CXXFLAGS=-g




all: $(NAME)


$(NAME): $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS) $(LDLIBS)






clean:
	rm -f *.o



