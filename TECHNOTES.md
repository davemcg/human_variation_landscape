## Notes on how ExAC calculated their gene-level Z-scores.

### Sources
http://www.nature.com/nature/journal/v536/n7616/extref/nature19057-s1.pdf
- ExAC paper supplement
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4222185/
- Daly paper used as basis for ExAC Z-score calculation



At the most basic level, this is a test of observed vs expected number of variants at the gene-level. 
But like everything, it is VERY complicated to do properly. You should try to control for the following stuff:

1. Sequencing depth
  * ExAC uses a median depth >50 at the exon level
  * For exons with <50 median, they calculated a regression using synonymous variants against depth
2. Divergence of region from primates
  * Hellmann, I. et al. Why do human diversity levels vary at a megabase scale? Genome Research 15, 1222- 1231 (2005).
3. Type of mutation
  * Synonymous, misssense, truncating
  * I am using missense-only right now
  * But should use synonymous-only as a control of sorts?
4. Filtering. ExAC: "For this paper, we focus on the canonical transcript as defined by Ensembl v75 for each protein-coding gene and drop all transcripts that do not begin with a methionine, end with a stop codon, or whose length are not divisible by three."
5. Outliers
  * Should ID weird windowns and just toss them
  * I def have some windows with HUGE >100 variants

This is how Samocha/Daly describe their approach:

We wanted to create an accurate model of de novo mutation for each gene. In order to do so, we extended a previous sequence context-based model of de novo mutation to derive gene-specific probabilities of mutation for each of the following mutation types: synonymous, missense, nonsense, essential splice site, and frameshift3. In brief, the local sequence context was used to determine the probability of each base in the coding region mutating to each other possible base and then determine the coding impact of each possible mutation. These probabilities of mutation were summed across genes to create a per-gene probability of mutation for the aforementioned mutation types (see Supplementary Note for more details). Here, we applied the method to exons and immediately flanking essential splice sites, but note that the framework is applicable to non-genic sequences. While fitting the expected rates of mutation to observed data, we added a term for local primate divergence across 1 Mb (to capture additional unmeasured sources of regional mutational variability) and another for the average depth of sequence of each nucleotide (to capture inefficiency of variant discovery at lower sequencing depths); both terms significantly improved the fit of the model to observed data (details in Supplementary Note). We also investigated a regional replication timing term22, but found no evidence for it significantly improving the model (Supplementary Note).


