INCLUDE_DIRECTORIES("${CMAKE_BINARY_DIR}/fem/src")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/hutiter/include")
INCLUDE(AddModules)

ADD_ELMER_MODULES(SKIP contrib)

ADD_SUBDIRECTORY(contrib)

ADD_DEPENDENCIES(EliminatePeriodic EliminateDirichlet)
ADD_DEPENDENCIES(FourierLoss MagnetoDynamics)
TARGET_LINK_LIBRARIES(EliminatePeriodic EliminateDirichlet)
