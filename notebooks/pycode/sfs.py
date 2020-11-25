
def ffs(mt, hl):
    # Calculate a histogram of alternate allele counts.
    n_sites, n_samples = mt.count()
    
    # The number of alleles in biallelic sites is: 2 * n_samples
    n_alleles = 2 * n_samples
    ac_hist = mt.aggregate_entries(hl.agg.filter(mt.alleles.length() == 2,
                  # Calculate a histogram of alternate allele counts.
                  hl.agg.hist(mt.variant_qc.AC[1], 1, n_alleles, n_alleles)))
    
    # Get the site counts from the histogram.
    site_counts = ac_hist.bin_freq
    
    # Total sites evaluated is the sum of the site counts.
    n_sites = sum(site_counts)
    
    # Site count divided by total sites.
    site_freq = [sc / n_sites for sc in site_counts]

    # Get the allele counts.
    allele_counts = [int(ac) for ac in ac_hist.bin_edges]
    
    # Calculate the folded frequency spectrum.
    folded_freq = []
    nn = len(site_freq)
    for ac, sf in zip(allele_counts, site_freq):
        ff = sf + site_freq[nn - ac]
        # The FFS is only defined for values corresponding to:
        if ac < nn / 2:
            folded_freq.append(ff)
        else:
            break
    
    # Make a Hail table with the allele counts and site frequencies.
    table = []
    for ac, fc in zip(allele_counts, folded_freq):
        row = {'ac': ac, 'ff': fc}
        table.append(row)

    ht_ffs = hl.Table.parallelize(hl.literal(table, 'array<struct{ac:int,ff:float32}>'))
    
    return ht_ffs