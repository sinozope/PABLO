#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# Add a target to generate the examples
option(BUILD_EXAMPLES "Create and install the examples" OFF)

include_directories("${PROJECT_SOURCE_DIR}/src")

if(BUILD_EXAMPLES)
	# List of examples
	set(EXAMPLE_LIST "")
	list(APPEND EXAMPLE_LIST "test1b")
	list(APPEND EXAMPLE_LIST "test3")
	list(APPEND EXAMPLE_LIST "test4")
	list(APPEND EXAMPLE_LIST "test5")
	list(APPEND EXAMPLE_LIST "test7")
	list(APPEND EXAMPLE_LIST "test13")
	list(APPEND EXAMPLE_LIST "test14")
	list(APPEND EXAMPLE_LIST "test15")
	list(APPEND EXAMPLE_LIST "test16")
	list(APPEND EXAMPLE_LIST "test17")
	list(APPEND EXAMPLE_LIST "test116")
	list(APPEND EXAMPLE_LIST "test120")
	list(APPEND EXAMPLE_LIST "bubbles")
	list(APPEND EXAMPLE_LIST "bubbles3D")

	#Rules to build the examples
	foreach(EXAMPLE_NAME IN LISTS EXAMPLE_LIST)
		set(EXAMPLE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/${EXAMPLE_NAME}.cpp")
		list(APPEND EXAMPLES_SOURCES ${EXAMPLE_SOURCES})

		add_executable(${EXAMPLE_NAME} "${EXAMPLE_SOURCES}")
		target_link_libraries(${EXAMPLE_NAME} ${PABLO_LIBRARY})
		install(TARGETS ${EXAMPLE_NAME} DESTINATION doc/${CMAKE_PROJECT_NAME}-${PABLO_VERSION}/examples)
	endforeach()

	add_custom_target(examples DEPENDS ${EXAMPLE_LIST})
	add_custom_target(clean-examples COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

	set (EXAMPLES_SOURCES ${EXAMPLES_SOURCES} PARENT_SCOPE)
endif()

