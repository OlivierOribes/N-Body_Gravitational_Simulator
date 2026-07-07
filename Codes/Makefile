# ---------------------------------------------------------------------------
# Build rules for the N-body solar system simulation.
#
# Layout:
#   Projet_Version_V8.cpp       -> final program (23 bodies, matrix method, Yoshida 4)
#   ProjetcodeV0.cpp            -> first stable version (matrix AND loop force methods)
#   Projet_Version_V6.cpp       -> intermediate version (NASA trajectory comparison)
#   Projet_Version_V7_norm.cpp  -> adds norm/angle/orbital-radius diagnostics
#   tools/*.cpp                 -> standalone data-extraction / post-processing helpers
#
# Eigen (header-only) is looked up in this order:
#   1. ./Eigen                      (vendored copy dropped next to the sources)
#   2. `pkg-config --cflags eigen3` (installed via package manager)
#   3. /usr/include/eigen3          (common default location on Linux)
# ---------------------------------------------------------------------------

CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra

EIGEN_INC := $(shell \
	if [ -d ./Eigen ]; then echo "-I."; \
	elif pkg-config --exists eigen3 2>/dev/null; then pkg-config --cflags eigen3; \
	elif [ -d /usr/include/eigen3 ]; then echo "-I/usr/include/eigen3"; \
	fi)

TARGET := simulation

.PHONY: all v0 v6 v7 tools run clean check-eigen

all: check-eigen $(TARGET)

check-eigen:
ifeq ($(strip $(EIGEN_INC)),)
	$(error Eigen headers not found. Install libeigen3-dev (or equivalent) \
	or drop an Eigen/ folder next to the Makefile. See README for details)
endif

# Final version -> executable "simulation"
$(TARGET): Projet_Version_V8.cpp
	$(CXX) $(CXXFLAGS) $(EIGEN_INC) $< -o $@

# Earlier report versions, built on demand for reference/reproducibility
v0: ProjetcodeV0.cpp
	$(CXX) $(CXXFLAGS) $(EIGEN_INC) $< -o $@

v6: Projet_Version_V6.cpp
	$(CXX) $(CXXFLAGS) $(EIGEN_INC) $< -o $@

v7: Projet_Version_V7_norm.cpp
	$(CXX) $(CXXFLAGS) $(EIGEN_INC) $< -o $@

# Post-processing helpers (no Eigen dependency)
tools: tools/dataextraction tools/differenceprog tools/extraction_diffnorm

tools/dataextraction: tools/dataextraction.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

tools/differenceprog: tools/differenceprog.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

tools/extraction_diffnorm: tools/Extraction_diffnorm.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) v0 v6 v7 tools/dataextraction tools/differenceprog tools/extraction_diffnorm
