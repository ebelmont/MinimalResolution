# Building yoneda2

Minimal build steps to get `yoneda2` (and its prerequisite data-generating tools) running.
For the full project's build options see `README.md`.

```
cmake -B build -S .
cmake --build build --target yoneda2 mr_ex mr_mot motTab e2p -j4
/bin/cp -f build/yoneda2 build/mr_ex build/mr_mot build/motTab build/e2p .
```

- `mr_ex`, `mr_mot`, `motTab`, `e2p` generate the resolution/multiplication-table data files
  that `yoneda2` reads; build them alongside it.
- The last `cp` step is required: `cmake --build` only updates binaries under `build/`, not
  the repo-root copies actually invoked (e.g. `./yoneda2`) — re-run it after every rebuild.
- If configure fails with a missing-source error for the `yoneda` target, that target was
  removed (see git log on CMakeLists.txt) because `yoneda.cpp`, the deprecated prior
  implementation, is intentionally absent from this checkout and out of scope.

Minimal example for running:

./build/e2p 40  
./build/motTab 40
./build/mr_ex 40 32
./build/mr_mot 40 30
./yoneda2 40 30 4 6 6 9

which computes {4-6} * {6-9}. E.g. you can check this equals ./yoneda2 40 30 6 9 4 6

There is also a "table mode"
./yoneda2 40 30 4 6
which computes {4-6} times all generators in range.
