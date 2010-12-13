CILK = cilk++
CXX = g++
CFLAGS := -Wall -g
SRCS   := \
	icache.cpp \
	bmp_io.cpp \
	light_source.cpp \
	raytracer.cpp \
	scene_object.cpp \
	util.cpp \
	random.cpp \
	photonmap.cpp \
	perlin.cpp \
	main.cpp
OBJS := $(SRCS:.cpp=.o)

TEST_SRCS := \
	util_test.cpp
# Add more test source files here.

# Objects to link into test binary.  Can't link raytracer.o into final test
# binary, because it defines main.
TEST_OBJS := $(TEST_SRCS:.cpp=.o) $(filter-out main.o,$(OBJS))
# Rename test objs to .test.o instead of .o so we can build non-cilk object
# files for testing.
TEST_OBJS := $(TEST_OBJS:.o=.test.o)

# Points to the root of Google Test, relative to where this file is.
GTEST_DIR = gtest

# We start our mode with "mode" so that we don't have leading whitespace from
# the +=.
NEWMODE := mode

# This is the usual DEBUG mode trick.  Add DEBUG=1 to your make command line to
# build without optimizations and with assertions.
ifeq ($(DEBUG),1)
  CFLAGS := -DDEBUG -O0 $(CFLAGS)
  NEWMODE += debug
else
  CFLAGS := -O3 -ffast-math -DNDEBUG $(CFLAGS)
  NEWMODE += nodebug
endif

# Perf relies on being able to take quick stack traces by walking the stack
# with the frame pointer.  If you want to profile your code with perf, add
# PROFILE=1 to your make command line.
ifeq ($(PROFILE),1)
  NOCILK := 1
  CFLAGS := -fno-omit-frame-pointer $(CFLAGS)
  NEWMODE += profile
endif

# Add NOCILK=1 to your make command if you want to build without the cilk
# runtime.
ifeq ($(NOCILK),1)
  CILKFLAGS := -fcilk-stub
  LDFLAGS := -fcilk-stub $(LDFLAGS)
  NEWMODE += nocilk
endif

# Add GPROF=1 to your make command if you want to add -pg to CFLAGS
# and LDFLAGS.
ifeq ($(GPROF),1)
  CFLAGS := -pg $(CFLAGS)
  LDFLAGS := -pg $(LDFLAGS)
  NEWMODE += gprof
endif

CILKFLAGS += $(CFLAGS)

# If the new mode doesn't match the old mode, write the new mode to .buildmode.
# This forces a rebuild of all the object files, because they depend on
# .buildmode.
OLDMODE := $(shell cat .buildmode 2> /dev/null)
ifneq ($(OLDMODE),$(NEWMODE))
  $(shell echo "$(NEWMODE)" > .buildmode)
endif

all: raytracer all_tests

# We only really need all_tests, but it would be weird if we had a raytracer
# binary from an old build configuration, so to be safe, we depend on all.
test: all
	./all_tests

# Rule for computing header dependencies.
%.d : %.cpp
	$(CILK) -M $< > $@

# Rule for compiling cpp files.
$(OBJS) : %.o : %.cpp .buildmode Makefile
	$(CILK) $(CILKFLAGS) -c $< -o $@

# Rule for linking object files into the final binary.
raytracer: $(OBJS)
	$(CILK) $(LDFLAGS) $(OBJS) -o $@

# Rule for compiling test cpp files.
$(TEST_OBJS) : %.test.o : %.cpp .buildmode Makefile
	$(CILK) -fcilk-stub -I$(GTEST_DIR)/include $(CFLAGS) -c $< -o $@

# Rule for linking object files into the final binary.
all_tests: $(TEST_OBJS) gtest_main.a
	$(CXX) -lpthread $(TEST_OBJS) gtest_main.a -o $@

# Includes compiler-generated .d files which track dependencies from object
# files to headers.
-include $(SRCS:.cpp=.d)

clean:
	rm -f *.o *.d raytracer .buildmode

## Google Test's targets.

# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
gtest-all.o : $(GTEST_SRCS_)
	$(CXX) -I$(GTEST_DIR) -I$(GTEST_DIR)/include $(CFLAGS) -c \
		$(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) -I$(GTEST_DIR)/include $(CFLAGS) -c $(GTEST_DIR)/src/gtest_main.cc

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^
