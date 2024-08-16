# Compilation of Elmer under Nix

## Compiling from upstream

The Nix flake exposes derivations (packages) that can be built using `nix build github:ElmerCSC/elmerfem#<version>`.

### List of included derivations

#### `default`
Barebones Elmer with OpenMP and MPI support.

#### `gui`
Elmer with GUI.

#### `full`
A derivation with more parallel computing features such as HYPRE.

All derivations have Elmer/Ice.

Specifying no version builds the `default` derivation.

After building, the binaries are accessible under `result/bin`.

## Compiling from a local repository

In the local Elmer repository run `nix build` or `nix build .#<version>`.

## Verbose output

By default `nix build` only shows one line of the build output at once.
To enable more verbose output, add the `-L` flag to the end of the build command.
