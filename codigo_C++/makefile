obj = $(wildcard %.o)
%.o:	%.cpp 
	g++ -c $< 
output:	main.o HCV_Model.o
	g++ *.o -o output	
run: output
	./output
clean:
	rm *.o
	rm output