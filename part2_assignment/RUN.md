# Build and Run Instructions (Ubuntu)

This guide provides the complete instructions for compiling and running the `factory` and `belts` executables.

## 1. Prerequisites

First, install the necessary build tools (`g++`, `make`, `cmake`, and `git`) with a single command:

```bash
sudo apt update
sudo apt install build-essential cmake git
```

# 1. Clone your repository (if you haven't already)
git clone [https://github.com/kawaljeets22/factorio-assignment-kawaljeet.git]
cd factorio-assignment-final

# 2. Create a 'build' directory and move into it
mkdir build
cd build

# 3. Configure the project using CMake
#
#    !! -- IMPORTANT -- !!
#
#    This *first* command will take a very long time.
#    (Est. 10-30+ minutes, depending on your CPU and internet)
#
#    It is downloading and *compiling all of Google OR-Tools'
#    dependencies (like Abseil and Protobuf) from source.
#    This is normal. Please be patient and let it finish.
#    You will see a lot of text.
#
cmake .. -DCMAKE_BUILD_TYPE=Release

# 4. Compile your actual code
#    (This will be much faster, but may still take a minute)
cmake --build .

# Replace '../my_factory_test.json' with the path to your input file
cat ../my_factory_test.json | ./part2_assignment/factory/factory > factory_output.json

# Replace '../my_belts_test.json' with the path to your input file
cat ../my_belts_test.json | ./part2_assignment/belts/belts > belts_output.json
