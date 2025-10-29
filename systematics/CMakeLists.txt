add_library(systematics ElectronSystematics.cpp)
target_link_libraries(systematics ${ROOT_LIBRARIES} util MAT MAT-MINERvA)
install(TARGETS systematics DESTINATION lib)
