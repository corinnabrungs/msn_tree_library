# Well plate visualization

If you used plates for your analysis you can visualize the compound detection depending on the
acquisition method, for example positive and negative ionization mode. An example for a 96 and 384
well
plate is shown below that were acquired in positive and negative ion mode.

![Plate96](pictures\plate_example96.png)
![Plate384](pictures\plate_example384.png)

You need two notebooks for the visualization.

The first one [method_comparison](notebooks/comparison_positive_negative.ipynb) combines the
spectral libraries (.mgf files) and highlights the detection (positive, negative, both, missing).
Accordingly, you have to provide the corresponding .mgf files as well as the metadata sheet. The
variable lib is used to detect the start of the unique_sample_id within the created usi. The mapping
name can be adjusted to the method, we used polarity in our case. Just be careful, which dataframe
is used left or right for the mapping. After running the complete notebook, you create a new file (
stored under your outfile path and name) that you need for the second script.

The second one [well_plate_visualization](notebooks/piechart_dataframe_wellplates.ipynb) is
creating the visualization. For the filename you use the file created before, you adjust the plate
size (so far only 96 and 384 is possible), and the name of the groups depends on the mapping you
used from the previous script. The outfile path and name can be changed in the code.



