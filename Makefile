release: primecount.cpp
	g++ -O3 -Wall -Wextra -std=c++20 -o primecount primecount.cpp

debug: primecount.cpp
	g++ -g -Wall -Wextra -std=c++20 -D_GLIBCXX_DEBUG -DDEBUG -fwrapv -o primecount primecount.cpp
