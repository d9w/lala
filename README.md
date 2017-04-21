# Learning Aquatic Locomotion with Animats

This is the corresponding repository for the submission "Learning Aquatic
Locomotion with Animats" to the European Conference on Artificial Life, 2017.
All of the included materials are subject to change depending on reviewer
feedback, and rights are reserved to the authors until further notice following
the paper acceptance decision.

Videos of the three morphologies can be found as follows:
+[Worm](https://vimeo.com/214087988)
+[Quadropus](https://vimeo.com/214086348)
+[Stingray](https://vimeo.com/214084084)

[MecaCell](https://github.com/jdisset/MecaCell/) should be installed first in a
different directory. The commit used in this work for the submission was
[3435975](https://github.com/jdisset/MecaCell/commit/3435975c89064ceef1fc13935a1d0bc311cf9417),
on the dev branch. MecaCell should be compiled following its instructions and
the directory containing MecaCell should then be referenced in CMakeLists.txt,
line 7:

```cmake
set(MecaCellDir ~/MecaCell)
```

To compile the LALA source code, use CMake:

```bash
mkdir build
cd build
cmake ..
make
cd ..
```

If this has properly compiled, you can run the results of the evolution
presented in the paper:

```bash
./bin/console -f config/best_worm.json
```

You can also use the MecacellViewer to visualize the results:

```bash
./bin/viewer -f config/best_stingray.json
```

CMA-ES is provided using [julia](https://julialang.org). To run the evolution, you must
first create a ``shape_defaults.json`` file containing the defaults for the
desired morphology. CMA-ES is then launched as follows:

```bash
cp config/worm_defaults.json config/shape_defaults.json
julia runcmaes.jl 1234
```
