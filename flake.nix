{
  description = "Elmer FEM";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.05";

    flake-utils.url = "github:numtide/flake-utils";

    nix-filter.url = "github:numtide/nix-filter";

    mumps = {
      url = "github:mk3z/mumps";
      inputs.nixpkgs.follows = "nixpkgs";
    };

    hypre = {
      url = "github:mk3z/hypre";
      inputs.nixpkgs.follows = "nixpkgs";
    };

    csa = {
      url = "github:mk3z/csa-c";
      inputs.nixpkgs.follows = "nixpkgs";
    };

    nn = {
      url = "github:mk3z/nn-c";
      inputs.nixpkgs.follows = "nixpkgs";
    };

    mmg = {
      url = "github:mk3z/mmg/develop";
      inputs.nixpkgs.follows = "nixpkgs";
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
        hypre = inputs.hypre.packages.${system}.default;
        csa = inputs.csa.packages.${system}.default;
        nn = inputs.nn.packages.${system}.default;
        mmg = inputs.mmg.packages.${system}.default;

        basePkg = {
          name,
          nativeBuildInputs,
          buildInputs,
          cmakeFlags,
          doCheck,
          checkPhase ? ''
            runHook preCheckPhase
            ctest -j $NIX_BUILD_CORES -L fast
            runHook postCheckPhase
          '',
        } @ inputs: let
          storepath = placeholder "out";
        in
          pkgs.stdenv.mkDerivation {
            inherit name doCheck checkPhase storepath;

            pname = "elmerfem";

            src = nix-filter {
              root = self;
              exclude = [
                (nix-filter.lib.matchExt "nix")
                "flake.lock"
                ".git"
                ".gitignore"
                ".gitmodules"
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

            autoPatchelfIgnoreMissingDeps = ["libmpi_stubs.so"];

            preConfigure = ''
              patchShebangs ./
            '';

            nativeCheckInputs = with pkgs; [
              mpiCheckPhaseHook
              openssh
            ];
          };

        nogui = {doCheck ? false}:
          basePkg {
            inherit doCheck;
            name = "elmerfem";
            nativeBuildInputs = [];
            buildInputs = [];
            cmakeFlags = [];
          };

        gui = {doCheck ? false}:
          basePkg {
            inherit doCheck;
            name = "elmerfem-gui";

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

        ice = {doCheck ? false}:
          basePkg {
            inherit doCheck;
            name = "elmerice";

            nativeBuildInputs = [];

            buildInputs = with pkgs;
              [
                hdf5-mpi
                scalapack
              ]
              ++ [
                csa
                hypre
                mumps
                nn
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
              "-DNN_INCLUDE_DIR=${nn}/include"
              "-DNN_LIBRARY=${nn}/lib/libnn.a"

              "-DWITH_MMG:BOOL=TRUE"
              "-DMMG_INCLUDE_DIR=${mmg}/include"
              "-DMMG_LIBRARY=${mmg}/lib/libmmg.so"

              "-DWITH_GridDataReader:BOOL=TRUE"

              "-DWITH_Trilinos:BOOL=FALSE"
            ];

            checkPhase = ''
              runHook preCheckPhase
              ctest -j $NIX_BUILD_CORES -L fast -E "(Hydro_Coupled)|(Hydro_SedOnly)|(Proj_South)"
              runHook postCheckPhase
            '';
          };
      in {
        checks = {
          nogui = nogui {doCheck = true;};
          gui = gui {doCheck = true;};
          ice = ice {doCheck = true;};
        };

        packages = {
          default = nogui {};
          nogui = nogui {};
          gui = gui {};
          ice = ice {};
        };
      }
    );
}
