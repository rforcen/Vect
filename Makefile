all:
	@g++ -std=c++2a -O3 main.cpp Vect.cpp -o Vect -pthread -lm 
	@./Vect