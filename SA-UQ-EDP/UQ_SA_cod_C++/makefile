obj = $(wildcard *.o)
output:	$(obj)
	g++ $(obj).o -o output	
%.o:	%.cpp
	g++ -c $< -o $@
run: output
	./output
