# Installing/compiling on NixOS

[Nix](https://nixos.org) is a declarative package manager.

## Elmer variants included

### `default`
Barebones Elmer with OpenMP and MPI support.

### `gui`
Elmer with GUI.

### `full`
A derivation with more parallel computing features such as HYPRE.

All derivations have Elmer/Ice.

Specifying no variant implies the `default` derivation.

## Nix overlay

The Nix flake exposes an overlay which contains the three Elmer variants under the following package names: `elmer`, `elmer-gui`, `elmer-full`.

## Binary cache

Add the following to your Nix config to use the Elmer binary cache.

```nix
nix.settings = {
  substituters = [ "https://elmerfem.cachix.org" ];
  trusted-public-keys = [ "elmerfem.cachix.org-1:nWIb5JzEzC2/W6qiuaC0urJRG+S7KvTn9WatX43gkHk=" ];
};
```

## Building of Elmer with Nix

### Compiling from upstream

`nix build github:ElmerCSC/elmerfem` for the default variant, or `nix build github:ElmerCSC/elmerfem#<variant>`.

After building, the binaries are accessible under `result/bin`.

### Compiling in a local repository

In the local Elmer repository run `nix build` or `nix build .#<variant>`.

### Verbose output

By default `nix build` only shows one line of the build output at once.
To enable more verbose output, add the `-L` flag to the end of the build command.

## For developers

### Running tests

#### All tests
`nix flake check -L`

#### Specific variant
`nix eval --raw .#checks.x86_64-linux.<variant>`

### Updating nixpkgs

- Edit the nixpkgs URL in `flake.nix`
  - Find the line `nixpkgs.url = "github:NixOS/nixpkgs/nixos-XX.YY";`.
  - Change `XX.YY` to the new nixpkgs release.
- Run `nix flake update` in the root of the repository.
- Check that everything works by running `nix flake check -L`.

### Adding new dependencies or build flags

- Find the variant you want to patch.
  - The derivation declarations are in the form `<variant> = {doCheck ? false}: { ...`.
- Add compile time dependencies to `nativeBuildInputs`.
- Add runtime dependencies to `buildInputs`.
- Add CMake flags to `cmakeFlags`.
