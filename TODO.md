## To Do

1. ~~Adjust the PEVmean and CDmean algorithms such that the genetic variance is equal to V_u * #.polymorphic markers, not total markers~~
2. Calculate expected response to selection and compare to realized
3. Explore allele frequency
    1. How often do alleles at low frequency drift to higher frequency
    2. Can we track markers / QTL that move into intermediate allele frequency?
4. Filtering minor alleles more stringently and observing the results in the scenario of not updating the training population.
    1. Filter on 0.03 minimum MAF (original) and 0.10 MAF
    2. Plot the site frequency spectra of markers at each cycle for each of the two MAF filtering strategies
5. Try to diagnose why updating the training population in some way is better than not doing anything at all. This represents the greatest difference among results, so it would be worth investigating.
6. Try to calculate expected marker effects. Presumably this would be a function of the QTL effect and LD between the QTL and the marker. If successful, find the differences between the observed and expected.
