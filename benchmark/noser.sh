#!/bin/bash
nosetests --processes=16 --process-timeout=60 \
    OperatorTest:OperatorTest \
    ChandraTest:ChandraTest \
#    ChandraTest:AnalyticalSolutionTest \
#    ChandraTest_2:ChandraTest_2 \
#    CylinderTest:CylinderTest \


    #SolutionSuite \
