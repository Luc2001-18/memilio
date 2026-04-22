# Install script for directory: /home/luc/Project_YAM/memilio/cpp

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libmemilio.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libmemilio.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libmemilio.so"
         RPATH "/home/luc/Project_YAM/memilio/cpp/examples/lib:/home/luc/Project_YAM/memilio/cpp/examples/bin")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE SHARED_LIBRARY FILES "/home/luc/Project_YAM/memilio/cpp/examples/lib/libmemilio.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libmemilio.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libmemilio.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libmemilio.so"
         OLD_RPATH "/opt/ohpc/pub/spack/develop/opt/spack/linux-rocky8-zen2/gcc-12.2.0/hdf5-1.14.0-ufs5mrlgiarhgsxtu4wp73lp5f7ciiod/lib:/home/luc/Project_YAM/memilio/cpp/examples/lib:"
         NEW_RPATH "/home/luc/Project_YAM/memilio/cpp/examples/lib:/home/luc/Project_YAM/memilio/cpp/examples/bin")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/libmemilio.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "/home/luc/Project_YAM/memilio/cpp/memilio" FILES_MATCHING REGEX "/memilio\\/[^/]*\\/[^/]*\\.h$")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "/home/luc/Project_YAM/memilio/cpp/examples/memilio" FILES_MATCHING REGEX "/memilio\\/[^/]*\\/[^/]*\\.h$")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/cmake/memilio" TYPE FILE FILES
    "/home/luc/Project_YAM/memilio/cpp/examples/memilio-config-version.cmake"
    "/home/luc/Project_YAM/memilio/cpp/examples/memilio-config.cmake"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/luc/Project_YAM/memilio/cpp/examples/memilio/ad/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/memilio/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/abm/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/d_abm/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/ode_secir/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/ode_secirts/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/ode_secirvvs/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/lct_secir/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/lct_secir_2_diseases/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/glct_secir/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/ide_secir/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/ide_seir/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/ode_seir/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/ode_seir_vector/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/ode_seirv/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/ode_seirdb/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/ode_seair/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/ode_sir/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/sde_sir/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/sde_sirs/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/sde_seirvv/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/graph_abm/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/smm/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/hybrid/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/models/ode_mseirs4/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/examples/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/tests/cmake_install.cmake")
  include("/home/luc/Project_YAM/memilio/cpp/examples/sbml_model_generation/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/luc/Project_YAM/memilio/cpp/examples/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
