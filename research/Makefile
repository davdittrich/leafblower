# Epic-J WU-1 stub
EIGEN_INCLUDE := $(shell Rscript -e "cat(system.file('include', package='RcppEigen'))")

PKG_CXXFLAGS = -std=c++17 -O3 -mavx2 -fopenmp-simd -DSTRICT_R_HEADERS -fPIC -I$(EIGEN_INCLUDE)
PKG_LIBS = -llapack -lblas

.PHONY: all clean asan

all: leafblower_research.so

leafblower_research.so: cp_calib.cpp ipm_calib.cpp research_bridge.cpp
	PKG_CXXFLAGS="$(PKG_CXXFLAGS)" PKG_LIBS="$(PKG_LIBS)" R CMD SHLIB cp_calib.cpp ipm_calib.cpp research_bridge.cpp -o leafblower_research.so

clean:
	rm -f *.o *.so

asan: PKG_CXXFLAGS += -fsanitize=address,undefined -g -O1
asan: leafblower_research.so
