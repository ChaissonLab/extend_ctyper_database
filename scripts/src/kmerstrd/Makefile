CXX := g++
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra

TARGET := KmerStrd
SRCS := src/main.cpp

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) -o $@ $^ -Wl,-rpath=$(CONDA_PREFIX)/lib 

clean:
	rm -f $(TARGET)

