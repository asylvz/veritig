VERITIG_VERSION := "0.2"
VERITIG_UPDATE := "March 31, 2026"
VERITIG_DEBUG := 0
BUILD_DATE := "$(shell date)"

CXX=g++
CXXFLAGS = -O2 -Wall -std=c++17 -DVERITIG_VERSION=\"$(VERITIG_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DVERITIG_UPDATE=\"$(VERITIG_UPDATE)\" -DVERITIG_DEBUG=$(VERITIG_DEBUG)
TARGET_EXEC := veritig
BUILD_DIR := ./build
SRC_DIRS := ./src

SRCS := $(shell find $(SRC_DIRS) -name '*.cpp')
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@

$(BUILD_DIR)/%.cpp.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)

.PHONY: clean
