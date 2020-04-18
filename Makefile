# defaults
CXXFLAGS += -g -O3
CXXFLAGS += -Wall -Werror -fopenmp
FITS_SUPPORT = 1
CGAL_SUPPORT = 1

$(BUILD_PATH)/%.o: %.cpp

# modify defaults here if necessary
-include features.mk

# FITS support
ifeq ($(FITS_SUPPORT), 1)
    LDFLAGS += -lCCfits -lcfitsio
endif

# have CGAL for Voronoi diagrams?
ifeq ($(CGAL_SUPPORT), 1)
    LDFLAGS += -lCGAL -LCGAL_Core -lgmp
endif

CXXFLAGS += -std=c++11 -I.

MAKEFILES = \
    Makefile \
    features.mk \

TEST_SOURCES = \
    test_photo.cpp \
    test_tools.cpp \
    test_imt.cpp \

TEST_OBJECTS = $(TEST_SOURCES:%.cpp=%.o)

ALL_SOURCES = \
    $(TEST_SOURCES) \
    runtests.cpp \
    sersic.cpp \
    ppanalysis.cpp \
    lodepng.cpp \
    

BINARIES = \
	banana \
    runtests \
    sersic \
    ppanalysis
default: banana
	
all: .ts.mk.hpp .depend $(BINARIES)

# hack to make any code recompile if Makefile changes
.ts.mk.hpp: $(MAKEFILES)
	@touch $@

# automagic dependencies
.depend: $(ALL_SOURCES) *.hpp .ts.mk.hpp
	$(CXX) $(CXXFLAGS) -MM $(ALL_SOURCES) >.depend
-include .depend

# create features.mk if absent
features.mk:
	@touch $@

imganalysis: imganalysis.o
	$(CXX) $(CXXFLAGS) -o $@ imganalysis.o $(LDFLAGS)

banana: banana.o lodepng.o
	$(CXX) $(CXXFLAGS) -o $@ banana.o lodepng.o $(LDFLAGS)

convertRegions: convertRegions.o
	$(CXX) $(CXXFLAGS) -0 $@ convertRegions.o $(LDFLAGS)

ppanalysis: ppanalysis.o
	$(CXX) $(CXXFLAGS) -o $@ ppanalysis.o $(LDFLAGS)

runtests: $(TEST_OBJECTS) runtests.o
	$(CXX) $(CXXFLAGS) -o $@ $(TEST_OBJECTS) runtests.o $(LDFLAGS)

sersic: sersic.o
	$(CXX) $(CXXFLAGS) -o $@ sersic.o $(LDFLAGS)

test: runtests
	./runtests

clean:
	rm -f $(BINARIES) *.o .ts.mk.hpp

.PHONY: default all test clean
