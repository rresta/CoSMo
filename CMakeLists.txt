cmake_minimum_required(VERSION 3.7)
project(thesis_code)

set(CMAKE_CXX_STANDARD 14)

#include_directories("C:\\Users\\Rossella\\seqan\\include")
#include_directories(include)

#set (SEQA N_INCLUDE_PATH "C:\\Users\\Rossella\\seqan\\include\\seqan\\")
#INCLUDE_DIRECTORIES(${SEQAN_INCLUDE_DIRS})
#ADD_DEFINITIONS(${SEQAN_DEFINITIONS})
#C:\Users\Rossella\seqan\include
if(WIN32 OR WIN64)
    include_directories(..\\seqan\\include)
    include_directories(.\\include)
else()
    include_directories(../seqan/include)
    include_directories(./include)
endif()

set(SOURCE_FILES main.cpp)
add_executable(thesis_code ${SOURCE_FILES})

#TARGET_LINK_LIBRARIES(seqan ${SEQAN_LIBRARIES})