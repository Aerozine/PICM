build:
	cmake -B build-dbg -G Ninja -DCMAKE_BUILD_TYPE=Debug; cmake --build build-dbg

build-fast:
	cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release; cmake --build build 

run-green:
	./build/bin/PIC -c taylorgreen.json

run-test:
	./build-dbg/bin/PIC -c test.json

run-uni:
	./build-dbg/bin/PIC -c test-uniform.json

run-fast:
	./build/bin/PIC -c test.json

format:
	find . -name "*.cpp" -o -name "*.hpp" | xargs clang-format -i --style=LLVM

clean:
	rm -rf build* results*
