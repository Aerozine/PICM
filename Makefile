build:
	cmake -B build-dbg -G Ninja -DCMAKE_BUILD_TYPE=Debug; cmake --build build-dbg

build-fast:
	cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release; cmake --build build 

run-green:
	./build/bin/PIC -c taylorgreen.json

run-test:
	./build-dbg/bin/PIC -c test.json

run-fast:
	./build/bin/PIC -c test.json

run-first:
	./build-dbg/bin/PIC -c first.json

view:
	paraview --script=parashow.py 

format:
	find . -name "*.cpp" -o -name "*.hpp" | xargs clang-format -i --style=LLVM

clean:
	rm -rf build build-dbg results results_taylor_green
