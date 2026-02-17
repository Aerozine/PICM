# particle-in-cell-methods

To run , open a terminal 
```
cmake -G Ninja -B build
```
Now that the build sheme has been setup
To compile everything correctly
```
cmake --build build --config Release
```
And to execute it 
```
./build/bin/PIC -c <jsonfile>
```
Depedencies are handled into the CmakeList (fetched if not present)
- Nlohmann Json lib
