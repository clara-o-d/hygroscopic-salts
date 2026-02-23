#!/usr/bin/env python3
"""
Generate calculate_activity_temperature MATLAB files for salts with temperature dependence.
Reads temperature_sources_list.xlsx and fits polynomials to water activity data.
"""

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import re
import os

def sanitize_filename(name):
    """Convert salt name to valid filename."""
    # Remove special characters and replace spaces with underscores
    name = re.sub(r'[^\w\s]', '', name)
    name = re.sub(r'\s+', '_', name)
    return name

def fit_bivariate_polynomial(mf, temp, aw, max_degree=3):
    """
    Fit a bivariate polynomial to water activity data.
    aw = sum_{i,j} c_{ij} * mf^i * temp^j
    """
    try:
        from sklearn.preprocessing import PolynomialFeatures
        from sklearn.linear_model import LinearRegression
    except ImportError:
        return None, None, float('inf')
    
    best_model = None
    best_rmse = float('inf')
    best_degree = 1
    
    for degree in range(1, max_degree + 1):
        try:
            # Create polynomial features
            poly = PolynomialFeatures(degree=degree, include_bias=True)
            X = np.column_stack([mf, temp])
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

def extract_salt_data(df, col_idx, salt_name):
    """Extract temperature-dependent data for a specific salt."""
    
    # Headers are in row 5 (index 5)
    # Data starts at row 6 (index 6)
    
    # Column structure: molality, temp, osmotic_coeff, aw, mass_fraction
    mf_col = col_idx + 4  # mass fraction column
    temp_col = col_idx + 1  # temperature column
    aw_col = col_idx + 3  # water activity column
    
    data_rows = []
    for i in range(6, len(df)):
        row = df.iloc[i]
        
        # Check if we have data
        if pd.notna(row.iloc[aw_col]) and pd.notna(row.iloc[temp_col]) and pd.notna(row.iloc[mf_col]):
            try:
                aw = float(row.iloc[aw_col])
                temp = float(row.iloc[temp_col])
                mf = float(row.iloc[mf_col])
                
                # Validate data
                if 0 <= aw <= 1 and -50 <= temp <= 150 and 0 <= mf <= 1:
                    data_rows.append([mf, temp, aw])
            except (ValueError, TypeError):
                break
        else:
            # Check if we've reached the end of this salt's data
            if len(data_rows) > 0:
                break
    
    if not data_rows:
        return None, None, None
    
    # Convert to arrays
    data_array = np.array(data_rows)
    mf_data = data_array[:, 0]
    temp_data = data_array[:, 1]
    aw_data = data_array[:, 2]
    
    return mf_data, temp_data, aw_data

def generate_matlab_file(salt_name, mf_data, temp_data, aw_data, output_dir):
    """Generate MATLAB calculate_activity_temperature file for a salt."""
    
    # Sanitize filename
    filename = f"calculate_activity_temperature_{sanitize_filename(salt_name)}.m"
    filepath = os.path.join(output_dir, filename)
    
    # Fit bivariate polynomial
    model_data, degree, rmse = fit_bivariate_polynomial(mf_data, temp_data, aw_data, max_degree=4)
    
    if model_data is None:
        print(f"  Error: Could not fit polynomial for {salt_name}")
        return None, None
    
    model, poly = model_data
    
    # Get data ranges
    mf_min = np.min(mf_data)
    mf_max = np.max(mf_data)
    temp_min = np.min(temp_data)
    temp_max = np.max(temp_data)
    aw_min = np.min(aw_data)
    aw_max = np.max(aw_data)
    
    # Generate MATLAB function
    func_name = f"calculate_activity_temperature_{sanitize_filename(salt_name)}"
    
    with open(filepath, 'w') as f:
        # Function header
        f.write(f"function aw = {func_name}(mf, T)\n")
        f.write(f"% Calculate water activity for {salt_name} as a function of mass fraction and temperature\n")
        f.write(f"%\n")
        f.write(f"% Inputs:\n")
        f.write(f"%   mf - mass fraction of {salt_name}\n")
        f.write(f"%   T  - temperature (°C)\n")
        f.write(f"% Output:\n")
        f.write(f"%   aw - water activity\n")
        f.write(f"%\n")
        f.write(f"% Data range:\n")
        f.write(f"%   Mass fraction: {mf_min:.4f} to {mf_max:.4f}\n")
        f.write(f"%   Temperature: {temp_min:.1f}°C to {temp_max:.1f}°C\n")
        f.write(f"%   Water activity: {aw_min:.4f} to {aw_max:.4f}\n")
        f.write(f"%\n")
        f.write(f"% Fit quality:\n")
        f.write(f"%   Polynomial degree: {degree}\n")
        f.write(f"%   RMSE: {rmse:.6f}\n")
        f.write(f"%   Number of data points: {len(aw_data)}\n")
        f.write(f"\n")
        
        # Input validation
        f.write(f"% Validate inputs\n")
        f.write(f"if mf < {mf_min:.6f} || mf > {mf_max:.6f}\n")
        f.write(f"    warning('{salt_name} mass fraction outside calibrated range (%.4f to %.4f)', {mf_min:.6f}, {mf_max:.6f});\n")
        f.write(f"end\n")
        f.write(f"if T < {temp_min:.2f} || T > {temp_max:.2f}\n")
        f.write(f"    warning('Temperature outside calibrated range (%.1f°C to %.1f°C)', {temp_min:.2f}, {temp_max:.2f});\n")
        f.write(f"end\n")
        f.write(f"\n")
        
        # Bivariate polynomial
        f.write(f"% Bivariate polynomial fit: aw = f(mf, T)\n")
        f.write(f"% Polynomial degree: {degree}, RMSE: {rmse:.6f}\n")
        f.write(f"\n")
        
        # Get feature names and coefficients
        feature_names = poly.get_feature_names_out(['mf', 'T'])
        coeffs = model.coef_
        
        f.write(f"% Calculate water activity\n")
        f.write(f"aw = 0")
        
        for i, (fname, coeff) in enumerate(zip(feature_names, coeffs)):
            if abs(coeff) > 1e-15:  # Skip near-zero coefficients
                # Convert feature name to MATLAB expression
                matlab_expr = fname.replace('mf', 'mf').replace('T', 'T')
                matlab_expr = matlab_expr.replace('^', '.^').replace(' ', '.*')
                
                if matlab_expr == '1':
                    f.write(f" + {coeff:.12e}")
                else:
                    f.write(f" + {coeff:.12e}*{matlab_expr}")
        
        f.write(f";\n")
        f.write(f"\n")
        
        # Ensure physical bounds
        f.write(f"% Ensure water activity is between 0 and 1\n")
        f.write(f"aw = max(0, min(1, aw));\n")
        f.write(f"\n")
        f.write(f"end\n")
    
    print(f"  Generated: {filename}")
    print(f"    Degree: {degree}, RMSE: {rmse:.6f}, Data points: {len(aw_data)}")
    
    return filepath, rmse

def main():
    """Main function to process all salts with temperature data."""
    
    # Read the Excel file
    excel_file = 'temperature_sources_list.xlsx'
    print(f"Reading {excel_file}...")
    df = pd.read_excel(excel_file)
    
    # Create output directory
    output_dir = '/Users/clara/Downloads/RE_ ML for AWH Discussion/calculate_activity_temperature'
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}\n")
    
    # Identify salt columns (those with MW as column name)
    salt_info = []
    for i, col in enumerate(df.columns):
        if isinstance(col, (int, float)) and col > 10:  # Molecular weights
            # Get salt name from row 2
            salt_name_raw = df.iloc[2, i]
            if pd.notna(salt_name_raw):
                # Extract salt name (e.g., "Paper 1: LiCl" -> "LiCl")
                match = re.search(r':\s*(.+)', str(salt_name_raw))
                if match:
                    salt_name = match.group(1).strip()
                    salt_info.append((i, col, salt_name))
    
    print(f"Found {len(salt_info)} salts with temperature data\n")
    print("="*80)
    
    # Process each salt
    results = []
    
    for col_idx, mw, salt_name in salt_info:
        print(f"\nProcessing: {salt_name} (MW = {mw})")
        
        # Extract data
        mf_data, temp_data, aw_data = extract_salt_data(df, col_idx, salt_name)
        
        if mf_data is None:
            print(f"  Error: Could not extract data")
            continue
        
        print(f"  Data points: {len(aw_data)}")
        print(f"  Mass fraction range: {np.min(mf_data):.4f} to {np.max(mf_data):.4f}")
        print(f"  Temperature range: {np.min(temp_data):.1f}°C to {np.max(temp_data):.1f}°C")
        print(f"  Water activity range: {np.min(aw_data):.4f} to {np.max(aw_data):.4f}")
        
        # Generate MATLAB file
        filepath, rmse = generate_matlab_file(salt_name, mf_data, temp_data, aw_data, output_dir)
        
        if filepath:
            results.append({
                'salt': salt_name,
                'n_points': len(aw_data),
                'rmse': rmse,
                'mf_range': (np.min(mf_data), np.max(mf_data)),
                'temp_range': (np.min(temp_data), np.max(temp_data)),
                'file': os.path.basename(filepath)
            })
    
    # Print summary
    print("\n" + "="*80)
    print("\nSUMMARY")
    print("="*80)
    print(f"Successfully processed {len(results)} salts\n")
    
    for result in results:
        print(f"{result['salt']}")
        print(f"  File: {result['file']}")
        print(f"  Data points: {result['n_points']}")
        print(f"  Mass fraction: {result['mf_range'][0]:.4f} to {result['mf_range'][1]:.4f}")
        print(f"  Temperature: {result['temp_range'][0]:.1f}°C to {result['temp_range'][1]:.1f}°C")
        print(f"  RMSE: {result['rmse']:.6f}")
        print()
    
    print(f"All files saved to: {os.path.abspath(output_dir)}")

if __name__ == "__main__":
    main()
