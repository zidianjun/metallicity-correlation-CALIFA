# metallicity-correlation-CALIFA
Metallicity correlations in CALIFA nearby galaxies

The majority of the codes is galaxy_metallicity.py, possessing classes named after Galaxy and GalaxyFigure. The most time-consuming function is MCMC()
. Typically it takes 1h for one galaxy. Thus, multi-process is highly recommended for looping all the 100 galaxies.

config.py has the essential parameters.

Manual:

1. Run "python main.py" to get the "total_chain_*galaxy_name*.txt" files in ./output/output and the figures in ./figures.
2. Run "python concatenate.py" to get the "total_chain_*diag*.csv" files in ./output.
3. Run "python correlation_length.py" to get the "correlation_length.csv" file in ./output.
4. Run the other python codes to get the figures in the publication.
  4.1 Run "python scatter.py" to get Fig. 1 and Fig. 11.
  4.2 Use the functions met_map() and met_fluc of GalaxyFigure to get Fig. 2.
  4.3 Use the functions met_fluc_corr() and mcmc_plot() of GalaxyFigure to get Fig. 3 and Fig. 4.
  4.4 Use the function corner.corner() in MCMC.py (defaulted as commented, use it by uncommenting) to get Fig. 5.
  4.5 Use the function thres_func() in main.py to get the "*galaxy_name*_thres.csv" file.
