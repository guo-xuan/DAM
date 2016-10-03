################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/AssociationEvaluation.cpp \
../src/Config.cpp \
../src/DAM.cpp \
../src/DependentModule.cpp \
../src/GwasData.cpp \
../src/HashTable.cpp \
../src/HyperGroup.cpp \
../src/IndependentModule.cpp \
../src/MCMC.cpp \
../src/MapTable.cpp \
../src/Variant.cpp 

OBJS += \
./src/AssociationEvaluation.o \
./src/Config.o \
./src/DAM.o \
./src/DependentModule.o \
./src/GwasData.o \
./src/HashTable.o \
./src/HyperGroup.o \
./src/IndependentModule.o \
./src/MCMC.o \
./src/MapTable.o \
./src/Variant.o 

CPP_DEPS += \
./src/AssociationEvaluation.d \
./src/Config.d \
./src/DAM.d \
./src/DependentModule.d \
./src/GwasData.d \
./src/HashTable.d \
./src/HyperGroup.d \
./src/IndependentModule.d \
./src/MCMC.d \
./src/MapTable.d \
./src/Variant.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -fopenmp -std=c++11 -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


