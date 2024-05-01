#!/bin/bash

FILE_TARGETS='target_models.csv'

# Run all simulations
for i in $(cat $FILE_TARGETS); do
    echo
    echo "Processing ${i/,/}";
    echo
    bash autorun_neuronsim.sh  ${i/,/};
    bash autorun_analysis.sh  ${i/,/};
done

## Run all analyses
#for i in $(cat $FILE_TARGETS); do
#    echo
#    echo "Processing ${i/,/}";
#    echo
#    bash autorun_analysis.sh  ${i/,/};
#done

