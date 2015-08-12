##
#  CMake script for the step-1 tutorial program:
##

# Set the name of the project and target:
SET(TARGET "Stokes")

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
IF(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

# Declare all source files the target consists of:
SET(TARGET_SRC
  ${TARGET}.cc
  StokesProblem.cc
  ExactSolutionEx1.cc
  ExactSolutionEx2.cc
  ExactSolutionEx3.cc
  ExactSolutionEx4.cc
  SchurComplement.cc
  RightHandSide.cc
  Parameters.cc
  # You can specify additional files here!
  )


ENABLE_TESTING()
ADD_EXECUTABLE(stokes_test 
  stokes_test.cc
  StokesProblem.cc
  ExactSolutionEx1.cc
  ExactSolutionEx2.cc
  ExactSolutionEx3.cc
  ExactSolutionEx4.cc
  SchurComplement.cc
  RightHandSide.cc
  Parameters.cc
  )
ADD_TEST(NAME stokes_test COMMAND stokes_test)

# Usually, you will not need to modify anything beyond this point...

FIND_PACKAGE(deal.II 8.3.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_SETUP_TARGET(stokes_test)
DEAL_II_INVOKE_AUTOPILOT()
