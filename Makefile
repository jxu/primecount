release: primecount.cpp
	g++ -O3 -Wall -Wextra -o primecount primecount.cpp

debug: primecount.cpp
	g++ -g -Wall -Wextra -D_GLIBCXX_DEBUG -DDEBUG -fsanitize=signed-integer-overflow -o primecount primecount.cpp
