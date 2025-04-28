# Example data usage
## For MSanalyst test
* `example.mgf` and `example_quant.csv` are used for `mn.py` and  `reanalysing.py`.

* `319.mgf` is used for `ms2search.py`.

* `exampleDB.mgf` and `exampleDB.xlsx` are used for `customizing.py`

## For figures in manuscript
### Figure 2
* Data: `std.mgf` and `std_quant.csv`
* Library searching parameters: 
  1. --pepmass_match_tolerance 10
  2. --library_matching_method modified_cosine
  3. --library_matching_similarity 0.7 (≥ 0.7)
  4. --is_library_matching_similarity 0.7 (≥ 0.7)
  5. --peak_percentage_threshold 0.7 (≥ 0.7)
  6. --is_library_matching_peaks 5 (≥ 5)
  7. --library_matching_peaks 5 (≥ 5)

## Figure 3
* Data: `KutzOsmac.mgf` and `KutzOsmac_quant.csv`
* Library searching parameters: 
  1. --pepmass_match_tolerance 5
  2. --library_matching_method modified_cosine
  3. --library_matching_similarity 0.7 (≥ 0.7)
  4. --is_library_matching_similarity 0.7 (≥ 0.7)
  5. --peak_percentage_threshold 0.7 (≥ 0.7)
  6. --is_library_matching_peaks 5 (≥ 5)
  7. --library_matching_peaks 5 (≥ 5)

* Networking parameters: 
  1. Network 1:
     1. --self_clustering_method neutral_loss
     2. --self_clustering_similarity 0.7
     3. --self_clustering_peaks 3
  2. Network 2:
     4. --self_clustering_method symmetric_chi_squared
     5. --self_clustering_similarity 0.7
     6. --self_clustering_peaks 4
     
    


