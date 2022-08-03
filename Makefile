include make.inc
DIR=`pwd`
subdir=baselib sprel

all: baselib sprel

baselib: 
	( cd src/$@ ; $(MAKE) )

sprel: baselib 
	( cd src/$@ ; $(MAKE) )
