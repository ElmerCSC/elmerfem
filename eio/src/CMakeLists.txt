SET(COMMON_SRCS 
	eio-config.h 
	EIOConstraintAgent.cpp
	EIODualMeshAgent.cpp 
	EIOGeometryAgent.cpp 
	EIOMeshAgent.cpp 
	EIOModelDataAgent.cpp 
	EIOModelManager.cpp 
	EIOPartWriter.cpp 
	EIOSolverAgent.cpp)

INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/eio/include")
INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}")

ADD_LIBRARY(eioc eio_api_c.cpp ${COMMON_SRCS})
ADD_LIBRARY(eiof eio_api_f.cpp ${COMMON_SRCS})

