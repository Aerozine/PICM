build:
	cmake -B build-dbg -G Ninja -DCMAKE_BUILD_TYPE=Debug; cmake --build build-dbg
build_fast:
	cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release; cmake --build build 

run2:
	./build/bin/PIC -c test.json

run:
	./build-dbg/bin/PIC -c test.json

view:
	paraview --script=parashow.py 
format:
	find . -name "*.cpp" -o -name "*.hpp" | xargs clang-format -i --style=LLVM
