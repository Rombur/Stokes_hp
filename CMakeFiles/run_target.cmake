SET(ENV{PATH} "/u/aghesmati/Documents/Spring2014/Implementation_Stokes:/w/aghesmati/dealII_install/lib:$ENV{PATH}")
EXECUTE_PROCESS(COMMAND Stokes
  RESULT_VARIABLE _return_value
  )
IF(NOT "${_return_value}" STREQUAL 0)
  MESSAGE(SEND_ERROR "
Program terminated with exit code: ${_return_value}")
ENDIF()
