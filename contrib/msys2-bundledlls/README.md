# msys2-bundledlls
A bare bones Bash script to collect DLL files linked against an executable and copying them for bundling with the program.

It is an alternative to the [mingw-bundledlls](https://github.com/mpreisler/mingw-bundledlls) Python script that aims to do the same, but goes through all linked DLLs and removing ones found in the system using a blacklist. Using a blacklist for system libraries found in Windows is rather fragile, as it will need to be maintained to not cause false negatives.

Instead, this script works off a whitelist by filtering by DLLs that are found in the current MSYS2 environment (e.g. `/ucrt64/bin/`), which is a lot less error prone and will copy over any DLLs not provided by Windows. It's also a simple Bash script instead of written in Python and does not require having Python installed in the environment to run.

It is currently being used for [Principia](https://principia-web.se)'s Windows builds.

## Usage
Clone the repository or download the `bundledlls` script. If you trust me enough, you can download it directly using cURL like this:

```bash
curl https://raw.githubusercontent.com/rollerozxa/msys2-bundledlls/master/bundledlls > ../bundledlls
```

The first argument is the executable and the second argument is the directory in which the DLLs will be copied into. Usually this would be the same directory as the executable, but depending on your distribution pipeline you might want put them somewhere else for packaging.

For example, if you are building [Luanti in MSYS2](https://dev.luanti.org/compiling-on-windows-using-msys2/), you would want to run the following:

```bash
./bundledlls ../bin/luanti.exe ../bin/
```

## License
0BSD
