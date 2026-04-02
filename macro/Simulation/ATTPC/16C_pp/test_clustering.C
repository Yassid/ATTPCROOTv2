/// @file test_clustering.C
/// @brief Clustering spacing scan for UKF track fitting.
///
/// Scans different minimum cluster spacing values by calling run_ukf_digi
/// for each configuration. Due to ROOT 6.26 interpreter limitations with
/// complex C++ (unique_ptr moves), this scan is best run via the shell:
///
///   source build/config.sh
///   for s in 1.0 2.0 3.0 5.0 7.0 10.0 15.0; do
///     echo "=== Spacing $s mm ==="
///     root -b -q "run_ukf_digi.C(-1, 1.0, $s)" 2>&1 | grep -E "Fitted|Failed|Mom|Avg"
///   done
///
/// Results from 16C(p,p) at 20 deg, 200 events, 50 proton tracks:
///
///   | Spacing | Fitted | Bias(%) | RMS(%) | Clusters |
///   |---------|--------|---------|--------|----------|
///   | 1.0 mm  |   50   |  -1.56  |  1.00  |    66    |
///   | 2.0 mm  |   50   |  -1.58  |  0.96  |    63    |  <-- best RMS
///   | 3.0 mm  |   49   |  -1.58  |  1.04  |    60    |
///   | 5.0 mm  |   50   |  -1.67  |  1.05  |    53    |
///   | 7.0 mm  |   50   |  -1.73  |  1.22  |    50    |
///   | 10.0 mm |   50   |  -1.70  |  1.14  |    46    |
///   | 15.0 mm |   50   |  -1.78  |  1.35  |    43    |
///
/// Conclusions:
///   - Bias is stable (~-1.6%) across all spacings → energy loss model, not clustering
///   - Optimal spacing is ~2 mm for this track type (best RMS with 100% convergence)
///   - Diminishing returns above ~5 mm spacing
///   - 100% convergence for all spacings tested (except 3 mm: 98%)

void test_clustering()
{
   std::cout << "This macro documents clustering scan results." << std::endl;
   std::cout << "To run the scan, use the shell loop in the file header." << std::endl;
   std::cout << "Or run individual spacings:" << std::endl;
   std::cout << "  root -b -q 'run_ukf_digi.C(-1, 1.0, 2.0)'" << std::endl;
}
