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


## Suggestions from Genomewide Selection Group
1. ~~More closely match the breeding cycles to seasons.~~ This is already implemented in the experimental design. If crosses are made in the Y1 Fall, F1s are grown in the Y2 Winter GH, F2s are grown in Y2 Summer, and F3s in the Y2 Fall. Genotyping is done as the F3s grow, and predictions are made on the F3s, which can then be used as parents. While the next cycle breeding lines are being inbred, the previous cycle F4 lines can be grown in the Y3 Winter nursery, and the F5s can be evaluated in the field in the Y3 Summer. Data from these trials could be available for the following fall for updated predictions. Therefore, TP updating could occur during each cycle.
2. Measure diversity among parents to maintain V_g
3. Use more than 700 snps in the long term
4. Use procedures to optimize for accuracy and diversity
5. Repeat cycle 0 to cycle 1 to look for parent overlap
6. Different TP sizes for different optimization schemes
7. Justify the random mating procedure for crossing. A: This is recurrent selection, so random mating would not be beyond the realm of possibility, especially if the selection procedure was being used to improve the population as a whole, from which elite inbreds would be derived.
