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
