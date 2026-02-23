#!/usr/bin/env python3
"""
Generate calculate_activity MATLAB files for salt mixtures.
Reads mixture_sources_list.xlsx and fits polynomials to water activity data.
"""

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import re
import os

def sanitize_filename(name):
    """Convert mixture name to valid filename."""
    # Remove special characters and replace spaces/+ with underscores
    name = re.sub(r'[^\w\s+]', '', name)
    name = re.sub(r'[\s+]+', '_', name)
    return name

def parse_mixture_name(mixture_str):
    """Extract component names from mixture string."""
    # Example: "NaCl + LiCl" -> ["NaCl", "LiCl"]
    components = [c.strip() for c in mixture_str.split('+')]
    return components

def fit_polynomial(mf_data, aw_data, max_degree=4):
    """
    Fit polynomial to water activity vs mass fraction data.
    Returns best fit coefficients and degree.
    """
    best_coeffs = None
    best_degree = 1
    best_rmse = float('inf')
    
    # Try polynomials from degree 1 to max_degree
    for degree in range(1, max_degree + 1):
        try:
            # Fit polynomial: aw = a0 + a1*mf + a2*mf^2 + ...
            coeffs = np.polyfit(mf_data, aw_data, degree)
            
            # Calculate RMSE
            aw_pred = np.polyval(coeffs, mf_data)
            rmse = np.sqrt(np.mean((aw_data - aw_pred)**2))
            
            # Keep if better
            if rmse < best_rmse:
                best_rmse = rmse
                best_coeffs = coeffs
                best_degree = degree
                
        except Exception as e:
            print(f"  Warning: Could not fit degree {degree}: {e}")
            continue
    
    return best_coeffs, best_degree, best_rmse

def extract_mixture_data(df, mixture_col, mixture_name):
    """Extract data for a specific mixture from the dataframe."""
    col_idx = df.columns.get_loc(mixture_col)
    
    # Get the data block (typically 8 columns wide)
    mixture_data = df.iloc[:, col_idx:col_idx+8]
    
    # Find the header row (contains 'aw')
    header_row_idx = None
    for i in range(min(10, len(mixture_data))):
        row_vals = mixture_data.iloc[i].values
        if any('aw' in str(val).lower() or 'a_w' in str(val).lower() for val in row_vals):
            header_row_idx = i
            break
    
    if header_row_idx is None:
        return None, None, None
    
    # Get headers
    headers = [str(h).strip() for h in mixture_data.iloc[header_row_idx].values]
    
    # Find mass fraction and aw columns
    mf_cols = [i for i, h in enumerate(headers) if 'mass fraction' in h.lower()]
    aw_col = [i for i, h in enumerate(headers) if h.lower() in ['aw', 'a_w']]
    
    if not aw_col:
        return None, None, None
    
    aw_col_idx = aw_col[0]
    
    # Extract data rows
    data_rows = []
    for i in range(header_row_idx + 1, len(mixture_data)):
        row = mixture_data.iloc[i]
        
        # Check if aw value exists and is numeric
        if pd.notna(row.iloc[aw_col_idx]):
            try:
                aw_val = float(row.iloc[aw_col_idx])
                
                # Extract mass fractions
                mf_vals = []
                for mf_idx in mf_cols:
                    if pd.notna(row.iloc[mf_idx]):
                        mf_vals.append(float(row.iloc[mf_idx]))
                    else:
                        mf_vals.append(0.0)
                
                if len(mf_vals) == len(mf_cols):
                    data_rows.append(mf_vals + [aw_val])
                    
            except (ValueError, TypeError):
                break
        else:
            break
    
    if not data_rows:
        return None, None, None
    
    # Convert to arrays
    data_array = np.array(data_rows)
    
    # Get mass fraction columns and aw
    mf_data = data_array[:, :-1]  # All but last column
    aw_data = data_array[:, -1]    # Last column
    
    # Get component names from mixture name instead of headers (headers are incorrect)
    component_names = parse_mixture_name(mixture_name)
    
    # Make sure we have the right number of components
    if len(component_names) != mf_data.shape[1]:
        # Fall back to header extraction if mismatch
        component_names = []
        for mf_idx in mf_cols:
            header = headers[mf_idx]
            # Extract component name from "mass fraction (ComponentName)"
            match = re.search(r'\(([^)]+)\)', header)
            if match:
                component_names.append(match.group(1).strip())
            else:
                component_names.append(f"Component{mf_idx}")
    
    return mf_data, aw_data, component_names

def fit_bivariate_polynomial(mf1, mf2, aw, max_degree=3):
    """
    Fit a bivariate polynomial to water activity data.
    aw = sum_{i,j} c_{ij} * mf1^i * mf2^j
    """
    from sklearn.preprocessing import PolynomialFeatures
    from sklearn.linear_model import LinearRegression
    
    best_model = None
    best_rmse = float('inf')
    best_degree = 1
    
    for degree in range(1, max_degree + 1):
        try:
            # Create polynomial features
            poly = PolynomialFeatures(degree=degree, include_bias=True)
            X = np.column_stack([mf1, mf2])
            X_poly = poly.fit_transform(X)
            
            # Fit linear regression
            model = LinearRegression(fit_intercept=False)
            model.fit(X_poly, aw)
            
            # Calculate RMSE
            aw_pred = model.predict(X_poly)
            rmse = np.sqrt(np.mean((aw - aw_pred)**2))
            
            if rmse < best_rmse:
                best_rmse = rmse
                best_model = (model, poly)
                best_degree = degree
        except Exception as e:
            continue
    
    return best_model, best_degree, best_rmse

def generate_matlab_file(mixture_name, components, mf_data, aw_data, output_dir):
    """Generate MATLAB calculate_activity file for a mixture."""
    
    # Sanitize filename
    filename = f"calculate_activity_{sanitize_filename(mixture_name)}.m"
    filepath = os.path.join(output_dir, filename)
    
    # Determine number of components
    n_components = mf_data.shape[1]
    
    # Fit polynomials for each component's mass fraction vs water activity
    # We'll fit aw as a function of mass fractions
    # For 2-component mixtures, we can use multiple approaches
    
    if n_components == 2:
        # For binary mixtures, try both univariate and bivariate fits
        
        # Univariate fits
        fits = []
        for i in range(n_components):
            coeffs, degree, rmse = fit_polynomial(mf_data[:, i], aw_data, max_degree=4)
            fits.append({
                'component': components[i],
                'coeffs': coeffs,
                'degree': degree,
                'rmse': rmse,
                'type': 'univariate',
                'component_idx': i
            })
        
        # Bivariate fit
        try:
            from sklearn.preprocessing import PolynomialFeatures
            from sklearn.linear_model import LinearRegression
            
            bivariate_model, bi_degree, bi_rmse = fit_bivariate_polynomial(
                mf_data[:, 0], mf_data[:, 1], aw_data, max_degree=3
            )
            
            if bivariate_model is not None:
                fits.append({
                    'component': 'both',
                    'model': bivariate_model,
                    'degree': bi_degree,
                    'rmse': bi_rmse,
                    'type': 'bivariate'
                })
        except ImportError:
            pass  # sklearn not available, skip bivariate fit
        except Exception as e:
            print(f"  Warning: Bivariate fit failed: {e}")
        
        # Use the fit with lower RMSE as primary
        primary_fit = min(fits, key=lambda x: x['rmse'])
        
        if primary_fit['type'] == 'univariate':
            primary_idx = primary_fit['component_idx']
        else:
            primary_idx = None
        
        # Get data range
        mf_min = [np.min(mf_data[:, i]) for i in range(n_components)]
        mf_max = [np.max(mf_data[:, i]) for i in range(n_components)]
        aw_min = np.min(aw_data)
        aw_max = np.max(aw_data)
        
        # Generate MATLAB function
        func_name = f"calculate_activity_{sanitize_filename(mixture_name)}"
        
        with open(filepath, 'w') as f:
            # Function header
            f.write(f"function aw = {func_name}(mf1, mf2)\n")
            f.write(f"% Calculate water activity for {mixture_name}\n")
            f.write(f"% Component 1: {components[0]}\n")
            f.write(f"% Component 2: {components[1]}\n")
            f.write(f"%\n")
            f.write(f"% Inputs:\n")
            f.write(f"%   mf1 - mass fraction of {components[0]}\n")
            f.write(f"%   mf2 - mass fraction of {components[1]}\n")
            f.write(f"% Output:\n")
            f.write(f"%   aw - water activity\n")
            f.write(f"%\n")
            f.write(f"% Data range:\n")
            f.write(f"%   {components[0]}: {mf_min[0]:.4f} to {mf_max[0]:.4f}\n")
            f.write(f"%   {components[1]}: {mf_min[1]:.4f} to {mf_max[1]:.4f}\n")
            f.write(f"%   aw: {aw_min:.4f} to {aw_max:.4f}\n")
            f.write(f"\n")
            
            # Input validation
            f.write(f"% Validate inputs\n")
            f.write(f"if mf1 < {mf_min[0]:.6f} || mf1 > {mf_max[0]:.6f}\n")
            f.write(f"    warning('{components[0]} mass fraction outside calibrated range');\n")
            f.write(f"end\n")
            f.write(f"if mf2 < {mf_min[1]:.6f} || mf2 > {mf_max[1]:.6f}\n")
            f.write(f"    warning('{components[1]} mass fraction outside calibrated range');\n")
            f.write(f"end\n")
            f.write(f"\n")
            
            if primary_fit['type'] == 'bivariate':
                # Bivariate polynomial fit
                model, poly = primary_fit['model']
                f.write(f"% Bivariate polynomial fit\n")
                f.write(f"% Degree: {primary_fit['degree']}, RMSE: {primary_fit['rmse']:.6f}\n")
                f.write(f"\n")
                
                # Get feature names and coefficients
                feature_names = poly.get_feature_names_out(['mf1', 'mf2'])
                coeffs = model.coef_
                
                f.write(f"% Calculate water activity\n")
                f.write(f"aw = 0")
                
                for i, (fname, coeff) in enumerate(zip(feature_names, coeffs)):
                    if abs(coeff) > 1e-15:  # Skip near-zero coefficients
                        # Convert feature name to MATLAB expression
                        matlab_expr = fname.replace('mf1', 'mf1').replace('mf2', 'mf2')
                        matlab_expr = matlab_expr.replace('^', '.^').replace(' ', '.*')
                        
                        if matlab_expr == '1':
                            f.write(f" + {coeff:.10e}")
                        else:
                            f.write(f" + {coeff:.10e}*{matlab_expr}")
                
                f.write(f";\n")
                
            else:
                # Univariate polynomial fit
                f.write(f"% Polynomial fit using {primary_fit['component']} mass fraction\n")
                f.write(f"% Degree: {primary_fit['degree']}, RMSE: {primary_fit['rmse']:.6f}\n")
                
                # Write coefficients
                coeffs = primary_fit['coeffs']
                for i, coeff in enumerate(coeffs):
                    f.write(f"A_{len(coeffs)-1-i} = {coeff:.10e};\n")
                
                f.write(f"\n")
                f.write(f"% Calculate water activity\n")
                
                # Use primary component for calculation
                if primary_idx == 0:
                    f.write(f"mf = mf1;  % Using {components[0]}\n")
                else:
                    f.write(f"mf = mf2;  % Using {components[1]}\n")
                
                # Build polynomial expression
                f.write(f"aw = A_0")
                for i in range(1, len(coeffs)):
                    f.write(f" + A_{i}*mf^{i}")
                f.write(f";\n")
            
            f.write(f"\n")
            f.write(f"end\n")
        
        print(f"  Generated: {filename}")
        if primary_fit['type'] == 'bivariate':
            print(f"    Fit type: bivariate (degree {primary_fit['degree']}, RMSE: {primary_fit['rmse']:.6f})")
        else:
            print(f"    Fit type: univariate - {primary_fit['component']} (degree {primary_fit['degree']}, RMSE: {primary_fit['rmse']:.6f})")
        
        return filepath, primary_fit['rmse']
    
    else:
        print(f"  Warning: {n_components}-component mixture not yet supported")
        return None, None

def main():
    """Main function to process all mixtures."""
    
    # Read the Excel file
    excel_file = 'mixture_sources_list.xlsx'
    print(f"Reading {excel_file}...")
    df = pd.read_excel(excel_file)
    
    # Create output directory
    output_dir = '/Users/clara/Downloads/RE_ ML for AWH Discussion/calculate_activity_mixtures'
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}\n")
    
    # Find all mixture columns
    mixture_cols = [col for col in df.columns if col.startswith('Paper') and ':' in col]
    
    print(f"Found {len(mixture_cols)} mixtures\n")
    print("="*80)
    
    # Process each mixture
    results = []
    
    for i, mixture_col in enumerate(mixture_cols, 1):
        # Extract mixture name
        match = re.search(r':\s*(.+)', mixture_col)
        if not match:
            continue
        
        mixture_name = match.group(1).strip()
        print(f"\n{i}. Processing: {mixture_name}")
        
        # Extract data
        mf_data, aw_data, components = extract_mixture_data(df, mixture_col, mixture_name)
        
        if mf_data is None:
            print(f"  Error: Could not extract data")
            continue
        
        print(f"  Components: {', '.join(components)}")
        print(f"  Data points: {len(aw_data)}")
        
        # Generate MATLAB file
        filepath, rmse = generate_matlab_file(mixture_name, components, mf_data, aw_data, output_dir)
        
        if filepath:
            results.append({
                'mixture': mixture_name,
                'components': components,
                'n_points': len(aw_data),
                'rmse': rmse,
                'file': os.path.basename(filepath)
            })
    
    # Print summary
    print("\n" + "="*80)
    print("\nSUMMARY")
    print("="*80)
    print(f"Successfully processed {len(results)} mixtures\n")
    
    for result in results:
        print(f"{result['mixture']}")
        print(f"  File: {result['file']}")
        print(f"  Components: {', '.join(result['components'])}")
        print(f"  Data points: {result['n_points']}")
        print(f"  RMSE: {result['rmse']:.6f}")
        print()
    
    print(f"All files saved to: {os.path.abspath(output_dir)}")

if __name__ == "__main__":
    main()
