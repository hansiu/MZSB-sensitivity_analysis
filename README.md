# MZSB-sensitivity_analysis

University project for a class. The task was to perform sensitivity analysis of any biological model.

I do not provide the final report and the analysis of the results here (especially as they are written in polish), just this scramble of the code that I made. Maybe someone will want to use pysb library with biomodels or SALib in the future and maybe this scramble could help them.

## Contents

In model.py file you can see how I used pysb to run the model of mouse Iron metabolism from biomodels and SALib to analyze the sensitivity (I used Morris, as it was possible with this data and Sobol to learn another method).

In analyze.py you can see how I obtained the plots for visualization of results in the report.
