export CC=/opt/local/bin/clang
export CXX=/opt/local/bin/clang++

all: build
	python setup.py install --user

build: python/*.pyx src/*.cc src/*.hh
	python setup.py build_ext -v

clean:
	python setup.py clean --all

.PHONY: build clean
