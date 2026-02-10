all: scripts/KmerStrd scripts/kmer_searcher


export
CPLUS_INCLUDE_PATH+=$(CONDA_PREFIX)/include
#export CPLUS_INCLUDE_PATH
CPATH+=$(CONDA_PREFIX)/include
#export CPATH
export LIBRARY_PATH := $(CONDA_PREFIX)/lib:$(LIBRARY_PATH)
export LD_LIBRARY_PATH := $(CONDA_PREFIX)/lib:$(LD_LIBRARY_PATH)


scripts/kmer_searcher:
	cd scripts/src/kmersearcher && make && cp kmer_searcher ../../


scripts/KmerStrd:
	cd scripts/src/kmerstrd && make && cp KmerStrd ../../

