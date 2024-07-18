{
  description = "Elmer FEM";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

    flake-utils.url = "github:numtide/flake-utils";

    nix-filter.url = "github:numtide/nix-filter";

    mumps = {
      url = "github:mk3z/mumps";
      inputs.nixpkgs.follows = "nixpkgs";
    };

    csa = {
      url = "github:mk3z/csa-c";
      inputs.nixpkgs.follows = "nixpkgs";
    };

    mmg = {
      url = "github:mk3z/mmg/develop";
      inputs.nixpkgs.follows = "nixpkgs";
    };

    parmmg = {
      url = "github:mk3z/parmmg/develop";
      inputs = {
        nixpkgs.follows = "nixpkgs";
        mmg.follows = "mmg";
      };
    };
  };

  outputs = {
    self,
    nixpkgs,
    flake-utils,
    nix-filter,
    ...
  } @ inputs:
    flake-utils.lib.eachDefaultSystem (
      system: let
        pkgs = nixpkgs.legacyPackages.${system};
        mumps = inputs.mumps.packages.${system}.default;
        csa = inputs.csa.packages.${system}.default;
        mmg = inputs.mmg.packages.${system}.default;
        parmmg = inputs.parmmg.packages.${system}.default;

        basePkg = {
          name,
          nativeBuildInputs,
          buildInputs,
          cmakeFlags,
          doCheck,
          checkOptions,
          ...
        } @ inputs: let
          storepath = placeholder "out";
        in
          pkgs.stdenv.mkDerivation {
            inherit name doCheck storepath;

            pname = "${name}-devel";

            src = nix-filter {
              root = self;
              exclude = [
                (nix-filter.lib.matchExt "nix")
                "flake.lock"
                ".git"
                ".github"
                ".gitignore"
                ".gitmodules"
                ".travis.yml"
                ".vscode"
              ];
            };

            hardeningDisable = ["format"];

            nativeBuildInputs = with pkgs;
              [
                cmake
                gfortran
                pkg-config
                autoPatchelfHook
              ]
              ++ inputs.nativeBuildInputs;

            buildInputs = with pkgs;
              [
                mpi
                blas
                liblapack
                tbb
              ]
              ++ inputs.buildInputs;

            cmakeFlags =
              [
                "-DELMER_INSTALL_LIB_DIR=${storepath}/lib"
                "-DCMAKE_INSTALL_LIBDIR=lib"
                "-DCMAKE_INSTALL_INCLUDEDIR=include"
                "-DWITH_LUA:BOOL=TRUE"

                "-DWITH_OpenMP:BOOLEAN=TRUE"
                "-DWITH_MPI:BOOLEAN=TRUE"

                "-Wno-dev"
              ]
              ++ inputs.cmakeFlags;

            checkPhase = ''
              runHook preCheckPhase
              ctest -j $NIX_BUILD_CORES -L ${checkOptions}
              runHook postCheckPhase
            '';

            autoPatchelfIgnoreMissingDeps = ["libmpi_stubs.so"];

            preConfigure = ''
              patchShebangs ./
            '';

            nativeCheckInputs = with pkgs; [
              mpiCheckPhaseHook
              openssh
            ];
          };

        default = {
          doCheck ? false,
          checkOptions ? "-L quick",
        }:
          basePkg {
            inherit doCheck checkOptions;
            name = "elmer";
            nativeBuildInputs = [];
            buildInputs = [];
            cmakeFlags = [];
          };

        gui = {
          doCheck ? false,
          checkOptions ? "-L quick",
        }:
          basePkg {
            inherit doCheck checkOptions;
            name = "elmer-gui";

            nativeBuildInputs = [pkgs.libsForQt5.wrapQtAppsHook];

            buildInputs = with pkgs; [
              libsForQt5.qtbase
              libsForQt5.qtscript
              libsForQt5.qwt
              libGL
              libGLU
              opencascade-occt_7_6
              vtkWithQt5
            ];

            cmakeFlags = [
              "-DWITH_ELMERGUI:BOOLEAN=TRUE"
              "-DWITH_QT5:BOOLEAN=TRUE"
              "-DWITH_OCC:BOOLEAN=TRUE"
              "-DWITH_VTK:BOOLEAN=TRUE"
              "-DCMAKE_OpenGL_GL_PREFERENCE=GLVND"
            ];
          };

        ice = {
          doCheck ? false,
          checkOptions ? ''-L fast -E "(Hydro_Coupled)|(Hydro_SedOnly)|(Proj_South)"'',
        }:
          basePkg {
            inherit doCheck checkOptions;
            name = "elmer-ice";

            nativeBuildInputs = [];

            buildInputs = with pkgs;
              [
                hdf5-mpi
                hypre
                nn
                scalapack
              ]
              ++ [
                csa
                mumps
              ];

            cmakeFlags = [
              "-DWITH_ElmerIce:BOOL=TRUE"

              "-DWITH_NETCDF:BOOL=TRUE"
              "-DNETCDF_LIBRARY=${pkgs.netcdf-mpi}/lib/libnetcdf.so"
              "-DNETCDFF_LIBRARY=${pkgs.netcdffortran}/lib/libnetcdff.so"
              "-DNETCDF_INCLUDE_DIR=${pkgs.netcdf-mpi}/include"
              "-DCMAKE_Fortran_FLAGS=-I${pkgs.netcdffortran}/include"

              "-DWITH_Hypre:BOOL=TRUE"

              "-DWITH_Mumps:BOOL=TRUE"

              "-DWITH_ScatteredDataInterpolator:BOOL=TRUE"
              "-DCSA_LIBRARY=${csa}/lib/libcsa.a"
              "-DCSA_INCLUDE_DIR=${csa}/include"
              "-DNN_INCLUDE_DIR=${pkgs.nn}/include"
              "-DNN_LIBRARY=${pkgs.nn}/lib/libnn.a"

              "-DWITH_MMG:BOOL=TRUE"
              "-DMMG_INCLUDE_DIR=${mmg}/include"
              "-DMMG_LIBRARY=${mmg}/lib/libmmg.so"

              "-DWITH_PARMMG:BOOL=TRUE"
              "-DPARMMG_INCLUDE_DIR=${parmmg}/include"
              "-DPARMMG_LIBRARY=${parmmg}/lib/libparmmg.so"

              "-DWITH_GridDataReader:BOOL=TRUE"

              "-DWITH_Trilinos:BOOL=FALSE"
            ];
          };
      in {
        checks = {
          default = default {doCheck = true;};
          gui = gui {doCheck = true;};
          ice = ice {doCheck = true;};
        };

        packages = {
          default = default {};
          gui = gui {};
          ice = ice {};
        };
      }
    );
}
