# CM124-EM

python3 em-phasing.py test_data.txt output.txt window_size

Recommended to run on window size less than 20 (accuracy increases with window size but program runtime increases exponentially)

To test the accuracy of the output of em-phasing.py, run the following command:
Rscript calculate switch accuracy.R [estimated haplotypes file] [true haplotypes file]