#!/bin/bash
nosetests \
    ChandraTest:ChandraTest \
    ChandraTest:AnalyticalSolutionTest \
    ChandraTest_2:ChandraTest_2 \
    CylinderTest:CylinderTest \
    OperatorTest:OperatorTest \
    --processes=16 --process-timeout=60
