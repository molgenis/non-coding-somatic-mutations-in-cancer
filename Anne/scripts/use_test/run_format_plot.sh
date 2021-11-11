#!/bin/bash

PATH_GLOBAL=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/
PATH_PLOT=/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/FORMAT/


while IFS= read -r line; do
  echo "tester: $line"
  python3 ${PATH_PLOT}make_plots_FORMAT.py ${line} ${PATH_PLOT}
done < "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/files.txt"