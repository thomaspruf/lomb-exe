# lomb-exe

The Lomb-Scargle periodogram is a highly sensitive tool in the search for rhythms in both evenly or unevenly sampled data.  Researchers may want to use it without having to install and learn a scripting language such as R, Julia, or Python. In all of these languages, Lomb-Scargle periodograms modules are indeed available.
Here, I provide executables for Mac, Linux, and Windows designed to compute the Lomb-Scargle periodogram for data stored in a comma separated file (.csv, as can be created e.g. in Microsoft Excel©.) Change the name of "lomb Mac" or "lomb Unix" to lomb.

The program is called at a minimum with
lomb data.csv      (use ./lomb data.csv on Mac or Linux)
where “data.csv” is an example file name in the same directory. Analysed will be the columns named t and y for measurement times and signal. Other columns will be ignored. To use other names, such as „times” and “activity” use switches (with double minus):
lomb mydata.csv --time=times --signal=activity
Main results will be printed:

Peak at     : 0.0416527824058    24.008	 Location of periodogram peak as frequency/period
Peak height : 0.4002705098900587	
P-value     : 9.892074577395284e-134       Probability that peak occurred by chance.
Sign. Level : 0.0174			 Significance level at  alpha
alpha       : 0.05			 Significance threshold
N out       : 600			 Number of periodogram points.


To add a downloadable plot in a browser window and to save the data on disk use switches –p and –s (single minus), respectively, as in:
lomb data.csv –p –s
This shows a plot of the periodogram power as a function of frequency. A peak that exceeds the significance level is probably not occurring by chance. The periodogram is always shown in the standardized form, i.e. powers range from 0 to 1.
The above command also creates the .csv file lomb.results, containing the main results followed by N out frequencies, followed by N out powers. If desired, this can be read into another program to plot power as a function of frequency.

To choose the density of frequencies to inspect, so as to not miss a narrow periodogram peak, use the oversampling factor switch (ofac). An ofac of 5, which is commonly used, means 5 times more frequencies are scanned:
lomb data.csv --ofac=5

To limit the range of frequencies to be inspected, use:
lomb data.csv  --from=1.0 --to=2.5

Using another significance threshold:
lomb data.csv --ofac=5 --from=1.0 --to=2.5 –p –s --alpha=0.01

Getting help:
lomb -h

Getting the version:
lomb -v

An example call using all options:
lomb  myfile.txt --time=timestamp --signal=tb -p -s --from=0.1 --to=1.0 
--ofac=4 --alpha=0.01 -h

Note
For a description of the properties of the Lomb-Scargle Periodogram, its computation and comparison with other methods see Ruf, T. (1999). Program lomb uses the algorithm given by Press et al (1994). The Lomb-Scargle periodogram was originally proposed by Lomb N.R. (1976) and further extended by Scargle J.D. (1982). An improved method for assessing the statistical significance of candidate periodicities by Baluev (2008), based on extreme value theory, is used. This implementation uses code modified from the Astropy.timeseries (© Copyright 2011–2021, The Astropy Developers.BSD3) Python package (VanderPlas et al. 2012, 2015).

Author
Thomas Ruf thomas.ruf@vetmeduni.ac.at based on code by Press et al (1994).

References
Baluev, R. V. (2008). Assessing the statistical significance of periodogram peaks. Monthly Notices of the Royal Astronomical Society, 385(3), 1279-1285.
Lomb N.R. (1976) Least-squares frequency analysis of unequally spaced data. Astrophysics and Space Science 39:447–462
Press W.H., Teukolsky S.A., Vetterling S.T., Flannery, B.P. (1994) Numerical recipes in C: the art of scientific computing.2nd edition. Cambridge University Press, Cambridge, 994pp.
Ruf, T. (1999) The Lomb-Scargle Periodogram in Biological Rhythm Research: Analysis of Incomplete and Unequally Spaced Time-Series. Biological Rhythm Research 30: 178–201.
Scargle J.D. (1982) Studies in astronomical time series. II. Statistical aspects of spectral analysis of unevenly spaced data. The Astrophysical Journal 302: 757–763.
VanderPlas, J., Connolly, A. Ivezic, Z. & Gray, A. (2012) Introduction to astroML: Machine learning for astrophysics. Proceedings of the Conference on Intelligent Data Understanding
VanderPlas, J. & Ivezic, Z. (2015) Periodograms for Multiband Astronomical Time Series.The Astrophysical Journal 812.1:18


