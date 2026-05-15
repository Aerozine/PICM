cmake -B build-release -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-release

find test/CH* -type f -name '*.json' | sort | while read -r f; do
	echo "Running: $f"
	./build-release/bin/PIC "$f"
done
