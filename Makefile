VERITIG_VERSION := "1.1.0"
VERITIG_UPDATE := "June 18, 2026"
VERITIG_DEBUG := 0
BUILD_DATE := "$(shell date)"

CXX=g++
CXXFLAGS = -O2 -Wall -std=c++17 -DVERITIG_VERSION=\"$(VERITIG_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DVERITIG_UPDATE=\"$(VERITIG_UPDATE)\" -DVERITIG_DEBUG=$(VERITIG_DEBUG)
TARGET_EXEC := veritig
BUILD_DIR := ./build
SRC_DIRS := ./src

SRCS := $(shell find $(SRC_DIRS) -name '*.cpp')
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

# -MMD -MP make the compiler emit a .d file listing the headers each object depends on,
# so editing a header rebuilds every object that includes it. Supported by both GCC and
# clang. Without this, changing a struct in a header leaves stale objects with the old
# layout, which links cleanly and misbehaves at run time.
$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@

$(BUILD_DIR)/%.cpp.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)

-include $(DEPS)

.PHONY: clean
