dbg:
	cmake -B build-dbg -G Ninja -DCMAKE_BUILD_TYPE=Debug; cmake --build build-dbg

release:
	cmake -B build-release -G Ninja -DCMAKE_BUILD_TYPE=Release; cmake --build build-release 

sl-square:
	./build/bin/PIC -c test/SemiLagrangian/test-square.json

sl-uniform:
	./build/bin/PIC -c test/SemiLagrangian/test-uniform.json

sl-small-cylinder:
	./build/bin/PIC -c test/SemiLagrangian/small-von-karman.json

sl:
	./build/bin/PIC -c test/bin/debug.json


sl-large-cylinder:
	./build/bin/PIC -c test/SemiLagrangian/large-von-karman.json

pic-default:
	./build/bin/PIC -c test/PIC/particles.json

pic-uniform:
	./build/bin/PIC -c test/PIC/uniform.json

heatsink:
	find test -name '*.json' | sort | while read -r f; do \
		echo "Running: $$f"; \
		./build-release/bin/PIC "$$f"; \
	done
heatsink-PIC:
	find test/PIC -name '*.json' | sort | while read -r f; do \
		echo "Running: $$f"; \
		./build-release/bin/PIC "$$f"; \
	done
 heatsink-APIC:
	find test/APIC -name '*.json' | sort | while read -r f; do \
		echo "Running: $$f"; \
		./build-release/bin/PIC "$$f"; \
	done
 heatsink-SL:
	find test/SL -name '*.json' | sort | while read -r f; do \
		echo "Running: $$f"; \
		./build-release/bin/PIC "$$f"; \
	done
 heatsink-FLIP:
	find test/FLIP -name '*.json' | sort | while read -r f; do \
		echo "Running: $$f"; \
		./build-release/bin/PIC "$$f"; \
	done

format:
	find . -name "*.cpp" -o -name "*.hpp" | xargs clang-format -i --style=LLVM

clean:
	find . -maxdepth 1 -type d \( -name 'build*' -o -name 'cmake-build*' \) -prune -exec rm -rf {} +
	if [ -d results ]; then \
		find results -type d \( -name raw -o -name runs -o -name configs -o -name plots \) -prune -exec rm -rf {} +; \
		find results -type f \( -name '*.vti' -o -name '*.vtp' -o -name '*.pvd' -o -name '*.png' -o -name '*.mp4' -o -name '*.log' -o -name '*.out' -o -name '*.err' \) -delete; \
	fi
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
memcheck:
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose ./build-dbg/bin/PIC -c test/test-large-cylinder.json
