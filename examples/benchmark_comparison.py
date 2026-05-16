"""
Benchmark comparison: sunwhere vs pvlib vs solposx

Compares performance for solar zenith angle calculation across:
- 1 site and 100 sites
- PSA and NREL algorithms
- NumExpr engine (when available)

Usage:
    python benchmark_comparison.py

Requirements:
    - sunwhere (required)
    - pvlib (optional): pip install pvlib
    - solposx (optional): pip install solposx
    - tabulate: pip install tabulate
"""

import time
import warnings
import numpy as np
import pandas as pd

try:
    from tabulate import tabulate
    HAS_TABULATE = True
except ImportError:
    HAS_TABULATE = False
    print("⚠️  tabulate not installed. Install with: pip install tabulate")
    print("   Results will be shown in plain text format.\n")

# Import libraries
import sunwhere

try:
    import pvlib.solarposition as pvsol
    HAS_PVLIB = True
except ImportError:
    HAS_PVLIB = False
    print("⚠️  pvlib not installed. Install with: pip install pvlib")

try:
    from solposx.solarposition import psa as solposx_psa
    from solposx.solarposition import spa as solposx_spa
    HAS_SOLPOSX = True
except ImportError:
    HAS_SOLPOSX = False
    print("⚠️  solposx not installed. Install with: pip install solposx")

# Suppress warnings
warnings.filterwarnings('ignore')


def benchmark_function(func, n_runs=10):
    """Measure execution time of a function"""
    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        func()
        end = time.perf_counter()
        times.append(end - start)
    return np.mean(times), np.std(times)


def run_benchmarks():
    """Run comprehensive benchmark comparison"""
    
    # Setup test data
    print("\n" + "="*70)
    print("BENCHMARK: sunwhere vs pvlib vs solposx")
    print("="*70)
    
    # # Time range: 1 year, hourly
    # times = pd.date_range('2024-01-01', periods=8760, freq='h', tz='UTC')
    # print(f"\nTime range: {len(times)} timestamps (1 year hourly)")
    
    # Time range: 1 year, minutely
    times = pd.date_range('2024-01-01', periods=60*24*365, freq='1min', tz='UTC')
    print(f"\nTime range: {len(times)} timestamps (1 year minutely)")

    # Test configurations
    n_sites_list = [1, 10]
    algorithms = ['psa', 'nrel']
    
    results = []
    
    for n_sites in n_sites_list:
        print(f"\n{'─'*70}")
        print(f"TESTING WITH {n_sites} SITE(S)")
        print(f"{'─'*70}")
        
        # Generate random locations
        if n_sites == 1:
            lats = 40.4  # Madrid
            lons = -3.7
        else:
            np.random.seed(42)
            lats = -90 + 180 * np.random.random(n_sites)
            lons = -180 + 360 * np.random.random(n_sites)
        
        for algorithm in algorithms:
            print(f"\n  Algorithm: {algorithm.upper()}")
            print(f"  {'-'*66}")
            
            # SUNWHERE with numexpr
            def sunwhere_calc():
                result = sunwhere.sites(
                    times, lats, lons, 
                    algorithm=algorithm, 
                    engine='numexpr',
                    refraction=False
                )
                _ = result.sza.values  # Force computation
            
            try:
                mean_time, std_time = benchmark_function(sunwhere_calc, n_runs=5)
                print(f"  ✓ sunwhere ({algorithm}, numexpr): {mean_time*1000:.2f} ± {std_time*1000:.2f} ms")
                results.append({
                    'Library': 'sunwhere',
                    'Sites': n_sites,
                    'Algorithm': algorithm,
                    'Engine': 'numexpr',
                    'Time (ms)': f"{mean_time*1000:.2f}",
                    'Std (ms)': f"{std_time*1000:.2f}"
                })
            except Exception as e:
                print(f"  ✗ sunwhere ({algorithm}, numexpr): ERROR - {e}")
            
            # SUNWHERE with numpy (only PSA for comparison)
            if algorithm == 'psa':
                def sunwhere_calc_numpy():
                    result = sunwhere.sites(
                        times, lats, lons, 
                        algorithm=algorithm, 
                        engine='numpy',
                        refraction=False
                    )
                    _ = result.sza.values  # Force computation
                
                try:
                    mean_time, std_time = benchmark_function(sunwhere_calc_numpy, n_runs=5)
                    print(f"  ✓ sunwhere ({algorithm}, numpy):   {mean_time*1000:.2f} ± {std_time*1000:.2f} ms")
                    results.append({
                        'Library': 'sunwhere',
                        'Sites': n_sites,
                        'Algorithm': algorithm,
                        'Engine': 'numpy',
                        'Time (ms)': f"{mean_time*1000:.2f}",
                        'Std (ms)': f"{std_time*1000:.2f}"
                    })
                except Exception as e:
                    print(f"  ✗ sunwhere ({algorithm}, numpy): ERROR - {e}")
            
            # PVLIB
            if HAS_PVLIB:
                # Map algorithm names
                pvlib_method_map = {
                    'psa': 'ephemeris',  # pvlib doesn't have PSA, using ephemeris
                    'nrel': 'nrel_numba'
                }
                pvlib_method = pvlib_method_map.get(algorithm)
                
                if pvlib_method:
                    if n_sites == 1:
                        def pvlib_calc():
                            result = pvsol.get_solarposition(
                                times, lats, lons,
                                method=pvlib_method
                            )
                            _ = result['zenith'].values
                        
                        try:
                            mean_time, std_time = benchmark_function(pvlib_calc, n_runs=5)
                            print(f"  ✓ pvlib ({pvlib_method}):          {mean_time*1000:.2f} ± {std_time*1000:.2f} ms")
                            results.append({
                                'Library': 'pvlib',
                                'Sites': n_sites,
                                'Algorithm': pvlib_method,
                                'Engine': 'numba' if 'numba' in pvlib_method else 'python',
                                'Time (ms)': f"{mean_time*1000:.2f}",
                                'Std (ms)': f"{std_time*1000:.2f}"
                            })
                        except Exception as e:
                            print(f"  ✗ pvlib ({pvlib_method}): ERROR - {e}")
                    else:
                        # pvlib requires iteration over multiple sites
                        def pvlib_calc_multi():
                            results_list = []
                            for i in range(n_sites):
                                lat_i = lats[i] if hasattr(lats, '__getitem__') else lats
                                lon_i = lons[i] if hasattr(lons, '__getitem__') else lons
                                result = pvsol.get_solarposition(
                                    times, lat_i, lon_i,
                                    method=pvlib_method
                                )
                                results_list.append(result['zenith'].values)
                            _ = np.array(results_list)
                        
                        try:
                            mean_time, std_time = benchmark_function(pvlib_calc_multi, n_runs=5)
                            print(f"  ✓ pvlib ({pvlib_method}):          {mean_time*1000:.2f} ± {std_time*1000:.2f} ms (iterating over sites)")
                            results.append({
                                'Library': 'pvlib',
                                'Sites': n_sites,
                                'Algorithm': pvlib_method,
                                'Engine': 'numba' if 'numba' in pvlib_method else 'python',
                                'Time (ms)': f"{mean_time*1000:.2f}",
                                'Std (ms)': f"{std_time*1000:.2f}"
                            })
                        except Exception as e:
                            print(f"  ✗ pvlib ({pvlib_method}): ERROR - {e}")
            
            # SOLPOSX
            if HAS_SOLPOSX:
                # Map algorithm names for solposx
                solposx_func_map = {
                    'psa': ('psa', solposx_psa),
                    'nrel': ('spa', solposx_spa)
                }
                
                if algorithm in solposx_func_map:
                    solposx_name, solposx_func = solposx_func_map[algorithm]
                    
                    if n_sites == 1:
                        def solposx_calc():
                            result = solposx_func(times, lats, lons)
                            _ = result['zenith']
                        
                        try:
                            mean_time, std_time = benchmark_function(solposx_calc, n_runs=5)
                            print(f"  ✓ solposx ({solposx_name}):            {mean_time*1000:.2f} ± {std_time*1000:.2f} ms")
                            results.append({
                                'Library': 'solposx',
                                'Sites': n_sites,
                                'Algorithm': solposx_name,
                                'Engine': 'python',
                                'Time (ms)': f"{mean_time*1000:.2f}",
                                'Std (ms)': f"{std_time*1000:.2f}"
                            })
                        except Exception as e:
                            print(f"  ✗ solposx ({solposx_name}): ERROR - {e}")
                    else:
                        # solposx requires iteration over multiple sites
                        def solposx_calc_multi():
                            results_list = []
                            for i in range(n_sites):
                                lat_i = lats[i] if hasattr(lats, '__getitem__') else lats
                                lon_i = lons[i] if hasattr(lons, '__getitem__') else lons
                                result = solposx_func(times, lat_i, lon_i)
                                results_list.append(result['zenith'])
                            _ = np.array(results_list)
                        
                        try:
                            mean_time, std_time = benchmark_function(solposx_calc_multi, n_runs=5)
                            print(f"  ✓ solposx ({solposx_name}):            {mean_time*1000:.2f} ± {std_time*1000:.2f} ms (iterating over sites)")
                            results.append({
                                'Library': 'solposx',
                                'Sites': n_sites,
                                'Algorithm': solposx_name,
                                'Engine': 'python',
                                'Time (ms)': f"{mean_time*1000:.2f}",
                                'Std (ms)': f"{std_time*1000:.2f}"
                            })
                        except Exception as e:
                            print(f"  ✗ solposx ({solposx_name}): ERROR - {e}")
    
    # Summary table
    print("\n" + "="*70)
    print("SUMMARY TABLE")
    print("="*70 + "\n")
    
    if results:
        df = pd.DataFrame(results)
        
        if HAS_TABULATE:
            print(tabulate(df, headers='keys', tablefmt='grid', showindex=False))
        else:
            # Fallback to simple formatting
            print(df.to_string(index=False))
        
        # Performance comparison
        print("\n" + "="*70)
        print("PERFORMANCE COMPARISON (relative to sunwhere)")
        print("="*70 + "\n")
        
        for n_sites in n_sites_list:
            print(f"\n{n_sites} site(s):")
            print("-" * 50)
            
            sunwhere_times = df[(df['Library'] == 'sunwhere') & (df['Sites'] == n_sites)]
            other_libs = df[(df['Sites'] == n_sites) & (df['Library'] != 'sunwhere')]
            
            # Map sunwhere algorithms to other libraries' algorithms
            algo_mapping = {
                'psa': ['ephemeris', 'default'],  # pvlib ephemeris, solposx default
                'nrel': ['nrel_numba']  # pvlib nrel_numba
            }
            
            for algo in algorithms:
                sw_time = sunwhere_times[sunwhere_times['Algorithm'] == algo]['Time (ms)'].values
                if len(sw_time) > 0:
                    sw_time = float(sw_time[0])
                    print(f"\n  {algo.upper()} algorithm:")
                    print(f"    sunwhere: {sw_time:.2f} ms (baseline)")
                    
                    # Get corresponding algorithms for other libraries
                    matching_algos = algo_mapping.get(algo, [])
                    
                    # Compare with other libraries for this algorithm
                    for _, row in other_libs.iterrows():
                        # Only compare if the algorithm matches
                        if row['Algorithm'] in matching_algos:
                            other_time = float(row['Time (ms)'])
                            speedup = other_time / sw_time
                            
                            lib_label = row['Library']
                        if n_sites > 1 and row['Library'] in ['pvlib', 'solposx']:
                            
                            if speedup > 1:
                                print(f"    {lib_label}: {other_time:.2f} ms ({speedup:.2f}x SLOWER)")
                            else:
                                print(f"    {lib_label}: {other_time:.2f} ms ({1/speedup:.2f}x FASTER)")
    
    print("\n" + "="*70)
    print("BENCHMARK COMPLETE")
    print("="*70 + "\n")


if __name__ == "__main__":
    run_benchmarks()
