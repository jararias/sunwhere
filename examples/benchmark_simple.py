"""
Simple benchmark: sunwhere performance test
Tests PSA and NREL algorithms with 1 and 100 sites
"""

import time
import numpy as np
import pandas as pd
import sunwhere


def simple_benchmark():
    """Quick performance test for sunwhere"""
    
    print("\n" + "="*70)
    print("SUNWHERE PERFORMANCE BENCHMARK")
    print("="*70)
    
    # Test data: 1 year hourly
    times = pd.date_range('2024-01-01', periods=8760, freq='h', tz='UTC')
    print(f"\nDataset: {len(times)} timestamps (1 year, hourly)")
    
    results = []
    
    for n_sites in [1, 100]:
        print(f"\n{'-'*70}")
        print(f"Testing with {n_sites} site(s)")
        print(f"{'-'*70}")
        
        # Generate locations
        if n_sites == 1:
            lats = 40.4
            lons = -3.7
        else:
            np.random.seed(42)
            lats = -90 + 180 * np.random.random(n_sites)
            lons = -180 + 360 * np.random.random(n_sites)
        
        for algorithm in ['psa', 'nrel']:
            print(f"\nAlgorithm: {algorithm.upper()}")
            
            # Warm-up run
            _ = sunwhere.sites(times, lats, lons, algorithm=algorithm, engine='numexpr')
            
            # Timed runs
            times_list = []
            for run in range(5):
                start = time.perf_counter()
                result = sunwhere.sites(times, lats, lons, algorithm=algorithm, engine='numexpr')
                _ = result.sza.values  # Force computation
                end = time.perf_counter()
                elapsed = (end - start) * 1000  # Convert to ms
                times_list.append(elapsed)
                print(f"  Run {run+1}: {elapsed:.2f} ms")
            
            mean_time = np.mean(times_list)
            std_time = np.std(times_list)
            print(f"  Average: {mean_time:.2f} ± {std_time:.2f} ms")
            
            results.append({
                'sites': n_sites,
                'algorithm': algorithm,
                'mean_ms': mean_time,
                'std_ms': std_time
            })
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70 + "\n")
    
    print(f"{'Sites':<10} {'Algorithm':<12} {'Mean (ms)':<15} {'Std (ms)':<15}")
    print("-" * 70)
    for r in results:
        print(f"{r['sites']:<10} {r['algorithm'].upper():<12} {r['mean_ms']:<15.2f} {r['std_ms']:<15.2f}")
    
    # Performance metrics
    print("\n" + "="*70)
    print("PERFORMANCE METRICS")
    print("="*70 + "\n")
    
    for algorithm in ['psa', 'nrel']:
        one_site = [r for r in results if r['sites'] == 1 and r['algorithm'] == algorithm][0]
        hundred_sites = [r for r in results if r['sites'] == 100 and r['algorithm'] == algorithm][0]
        
        speedup = hundred_sites['mean_ms'] / (one_site['mean_ms'] * 100)
        efficiency = speedup * 100
        
        print(f"{algorithm.upper()} Algorithm:")
        print(f"  1 site:     {one_site['mean_ms']:.2f} ms")
        print(f"  100 sites:  {hundred_sites['mean_ms']:.2f} ms")
        print(f"  Expected (100x1 site): {one_site['mean_ms']*100:.2f} ms")
        print(f"  Speedup: {1/speedup:.2f}x faster than sequential")
        print(f"  Efficiency: {efficiency:.1f}% (ideal vectorization = 100%)")
        print()
    
    # Throughput
    print("THROUGHPUT (calculations per second):")
    print("-" * 70)
    for r in results:
        total_calcs = len(times) * r['sites']
        throughput = total_calcs / (r['mean_ms'] / 1000)
        print(f"{r['algorithm'].upper():<12} {r['sites']:>3} sites: {throughput:>12,.0f} calcs/sec")
    
    print("\n" + "="*70)
    print("BENCHMARK COMPLETE")
    print("="*70 + "\n")


if __name__ == "__main__":
    simple_benchmark()
