build:
	cmake -B build-dbg -G Ninja -DCMAKE_BUILD_TYPE=Debug; cmake --build build-dbg

build-fast:
	cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release; cmake --build build 

run-green:
	./build/bin/PIC -c test/SemiLagrangian/taylorgreen.json

run-test:
	./build-dbg/bin/PIC -c test/SemiLagrangian/test.json

run-square:
	./build-dbg/bin/PIC -c test/SemiLagrangian/test-square.json

uniform:
	./build-dbg/bin/PIC -c test/SemiLagrangian/test-uniform.json

particles:
	./build-dbg/bin/PIC -c test/PIC/particles.json

cylinder:
	./build/bin/PIC -c test/SemiLagrangian/test-large-cylinder.json

source:
	./build-dbg/bin/PIC -c test/SemiLagrangian/test-source.json

run-fast:
	./build/bin/PIC -c test/SemiLagrangian/test.json

format:
	find . -name "*.cpp" -o -name "*.hpp" | xargs clang-format -i --style=LLVM

clean:
	rm -rf build* results*
memcheck:
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose ./build-dbg/bin/PIC -c test/test-large-cylinder.json
