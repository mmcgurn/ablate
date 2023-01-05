target_sources(ablateLibrary
        PRIVATE
        simpleFormula.cpp
        functionWrapper.cpp
        functionPointer.cpp
        constantValue.cpp
        parsedSeries.cpp
        linearTable.cpp
        formula.cpp
        formulaBase.cpp

        PUBLIC
        simpleFormula.hpp
        functionWrapper.hpp
        functionPointer.hpp
        mathFunction.hpp
        functionFactory.hpp
        constantValue.hpp
        parsedSeries.hpp
        linearTable.hpp
        formula.hpp
        formulaBase.hpp
        )

add_subdirectory(geom)