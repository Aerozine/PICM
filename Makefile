build:
	cmake -B build-dbg -G Ninja -DCMAKE_BUILD_TYPE=Debug; cmake --build build-dbg

build-fast:
	cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release; cmake --build build 

sl-square:
	./build-dbg/bin/PIC -c test/SemiLagrangian/test-square.json

sl-uniform:
	./build-dbg/bin/PIC -c test/SemiLagrangian/test-uniform.json

sl-small-cylinder:
	./build-dbg/bin/PIC -c test/SemiLagrangian/small-von-karman.json

sl-large-cylinder:
	./build-dbg/bin/PIC -c test/SemiLagrangian/large-von-karman.json

pic-default:
	./build-dbg/bin/PIC -c test/PIC/particles.json

pic-uniform:
	./build-dbg/bin/PIC -c test/PIC/uniform.json

format:
	find . -name "*.cpp" -o -name "*.hpp" | xargs clang-format -i --style=LLVM

clean:
	rm -rf build* results* cmake-build*
memcheck:
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose ./build-dbg/bin/PIC -c test/test-large-cylinder.json
