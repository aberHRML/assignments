# assignments 1.0.2

* add a fix for an error caused by a breaking change in [`tidygraph`](https://tidygraph.data-imaginist.com/index.html) v1.3.0.

# assignments 1.0.1

* The default ppm threshold has been reduced to 4.

* The importance of the ^37^Cl adduct has been increased in the default negative mode adducts.

* The ^13^C2 isotope has been removed from the default isotopes.

* The default retention time difference limit for relationships has been changed to 2 seconds for RP-LC-HRMS and NP-LC-HRMS.

* The absolute values of correlation coefficients are now used to calculate average component weights.

* Where components contain a feature represented by more than one adduct and isotope combination, only the node with the highest AIS is now retained.

# assignments 1.0.0

* Added a `NEWS.md` file to track changes to the package.

* The `Assignment` S4 class now inherits from the `AssignmentParameters` S4 class.

* The molecular formula generation is now handled by [`mzAnnotation::ipMF()`](https://aberhrml.github.io/mzAnnotation/reference/ipMF.html).

* Improved molecular formula selection routine based on the Seven Golden Rules from [Kind et al. 2007](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-105).

* The adduct and isotope assignment routine now conducted over multiple iterations.

* Graphical components are now selected using an improved plausibility score.

* Graphical components are now only retained if they contain at least one non-isotopic assignment.

* The individual assignment step methods (`calcCorrelations()`, `calcRelationships()`, `addIsoAssign()`, `transformAssign()`) are now exported.

* Added the `availableTechniques()` function to return the supported analytical techniques.

* Numerous documentation improvements.

* Added a usage introduction vignette.

* The package documentation is now available at <https://aberhrml.github.io/assignments/>
