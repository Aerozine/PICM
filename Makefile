build:
	cmake -B build-dbg -G Ninja -DCMAKE_BUILD_TYPE=Debug; cmake --build build-dbg

build-fast:
	cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release; cmake --build build 

run-green:
	./build/bin/PIC -c test/taylorgreen.json

run-test:
	./build-dbg/bin/PIC -c test/test.json

run-square:
	./build-dbg/bin/PIC -c test/test-square.json

uniform:
	./build/bin/PIC -c test/test-uniform.json

cylinder:
	./build/bin/PIC -c test/test-large-cylinder.json

source:
	./build/bin/PIC -c test/test-source.json


run-fast:
	./build/bin/PIC -c test/test.json

format:
	find . -name "*.cpp" -o -name "*.hpp" | xargs clang-format -i --style=LLVM

clean:
	rm -rf build* results*
