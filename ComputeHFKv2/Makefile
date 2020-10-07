ifeq ($(shell uname -s),Darwin)
CXXFLAGS := -std=c++11 -O3 -g -mmacosx-version-min=10.9
else
CXXFLAGS := -std=c++11 -O3 -g
endif

BUILD := build
LIB := ../lib
OBJS := $(addprefix $(BUILD)/, HFKLib.o Utility.o Min.o Max.o Crossing.o Simplify.o \
 HomologyRank.o Report.o Diagrams.o Alternating.o KnotFloer.o )
HDRS := Alg.h Diagrams.h

$(BUILD)/%.o : %.cpp
	g++ -c $(CXXFLAGS) $< -o $@


lib: $(OBJS) | $(LIB)
	ar -ruv $(LIB)/libhfk.a $(OBJS)

$(OBJS): $(HDRS) | $(BUILD)

$(BUILD):
	mkdir $(BUILD)

$(LIB):
	mkdir $(LIB)

.PHONY: clean
clean:
	rm -rf $(BUILD)


