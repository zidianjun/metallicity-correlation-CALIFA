# metallicity-correlation-CALIFA
Metallicity correlations in CALIFA nearby galaxies

The majority of the codes is galaxy_metallicity.py, possessing classes named after Galaxy and GalaxyFigure. The most time-consuming function is MCMC()
. Multi-processing is highly recommended for looping all the 100 galaxies. Typically it takes 400 core\*hour - if you set a 50-core multi-thread process, it will take about 8h (depends on the device).

config.py has the essential parameters.

Manual:

1. Run "python main.py" to get the "total_chain_*galaxy_name*.txt" files in ./output/output and the figures in ./figures.

2. Run "python concatenate.py" to get the "total_chain_*diag*.csv" files in ./output.

3. Run "python correlation_length.py" to get the "correlation_length.csv" file in ./output.

4. Run the other python codes to get the figures in the publication.

    4.1. Run "python fig1.py" to get Fig. 1.
  
    4.2. Use the functions example() to get Fig. 2.
  
    4.3. Use the functions met_fluc_corr() and mcmc_plot() of GalaxyFigure to get Fig. 3 and Fig. 4.
  
    4.4. Use the function corner.corner() in MCMC.py (defaulted as commented, use it by uncommenting) to get Fig. 5.
  
    4.5. Use the function thres_func() in main.py to get the "*galaxy_name*\_thres.csv" files. Run "python thres.py" to get Fig. 6.
    
    4.6. Run "python test.py", "python corner.py", "python density.py", "python vel_disp.py", "python morph.py", and "python corr_scale.py" to get Fig. 7 - 12.
    
    4.7. Set AGN_criterion = 'Kauffmann' and redo Step 1 - 3 to get the "correlation_length_Ka03.csv" file in ./output. Run "python lc_vs_lc.py --suffix='Ka03' " to get Fig. A1.
    
    4.8. Set decomposition = 'decom' and redo Step 1 - 3 to get the "correlation_length_MA17.csv" file in ./output. Run "python lc_vs_lc.py --suffix='MA17' " to get Fig. B1. Run "python decom.py" to get Fig. B2.
     
    4.9. Set adp_bin = True and redo Step 1 - 3 to get the "correlation_length_adp.csv" file in ./output. Run "python lc_vs_lc.py --suffix='adp' " to get Fig. C1.
