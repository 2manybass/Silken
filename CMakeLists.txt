cmake_minimum_required(VERSION 3.10)
project(Silken)

# Set standard C++ version
set(CMAKE_CXX_STANDARD 17)

# 1. Find the required external libraries
find_package(PkgConfig REQUIRED)
pkg_check_modules(GTK3 REQUIRED gtk+-3.0)
pkg_check_modules(SNDFILE REQUIRED sndfile)
find_library(FFTW3_LIB fftw3)
pkg_check_modules(GTKMM REQUIRED gtkmm-3.0)

# 2. Add your source files
add_executable(Silken
    main.cpp
    AnalysisEngine.cpp
    SpectralView.cpp
    NoteView.cpp
)

# 3. Include and Link everything
target_include_directories(Silken PRIVATE ${GTK3_INCLUDE_DIRS})
target_include_directories(Silken PRIVATE ${GTKMM_INCLUDE_DIRS})
target_link_libraries(Silken 
    ${GTK3_LIBRARIES} 
    ${SNDFILE_LIBRARIES} 
    ${FFTW3_LIB}
    ${GTKMM_LIBRARIES}
    m # link the math library
)
