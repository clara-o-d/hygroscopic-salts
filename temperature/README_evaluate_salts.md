# Atacama Desert Salt Evaluation for AWH

## Overview

The `evaluate_all_salts_atacama.m` script evaluates **all salts with temperature-dependent Pitzer parameters** to find the best candidates for Atmospheric Water Harvesting (AWH) in Atacama Desert conditions.

## Key Metric: Water Uptake Swing

**Water Uptake Swing** = (Max water uptake) - (Min water uptake) over the daily cycle

Where water uptake is defined as:
```
Water Uptake = kg H₂O / (kg H₂O + kg salt)
```

**Why this matters for AWH:**
- **High swing** = Large difference between night absorption and day release
- **More harvestable water** per kg of salt per day
- **Better AWH system efficiency**

## What the Script Does

### 1. **Loads All Salts**
- Reads `pitzer_binary.csv` for all cation-anion pairs
- Filters to only salts with temperature-dependent parameters
- Includes: LiCl, LiBr, NaCl, KCl, MgCl₂, CaCl₂, and many more

### 2. **Evaluates Atacama Performance**
For each salt:
- Calculates water uptake at all 82 Atacama (RH, T) points
- Uses temperature-dependent Pitzer model
- Accounts for different stoichiometries (1:1, 2:1, etc.)

### 3. **Calculates Swing**
- Finds minimum water uptake (typically during hot, dry day)
- Finds maximum water uptake (typically during cool, humid night)
- Computes difference: `swing = max - min`

### 4. **Ranks and Visualizes**
- Sorts salts by swing (best first)
- Generates 4 comparison plots
- Saves results to CSV

## Output Files

### Figures (saved to `figures/temperature/`)

1. **`Atacama_Salt_Comparison_Swing.png`**
   - Horizontal bar chart of top 15 salts
   - Ranked by water uptake swing
   - Shows which salts have best daily variation

2. **`Atacama_Salt_Swing_vs_MW.png`**
   - Scatter plot: Swing vs Molecular Weight
   - Color indicates data coverage
   - Top 5 salts labeled
   - Shows tradeoff between performance and weight

3. **`Atacama_Salt_Range_Comparison.png`**
   - Shows full range for each salt
   - Blue dot = minimum (day conditions)
   - Red dot = maximum (night conditions)
   - Line = daily range
   - Easy visual comparison

4. **`Atacama_Salt_Coverage.png`**
   - Bar chart showing % of Atacama points within valid range
   - Some salts don't work at very low RH (below deliquescence)
   - Higher coverage = more reliable across conditions

### Data File

**`atacama_salt_performance.csv`**
- Complete results table
- Columns: salt_name, cation, anion, water_uptake_min, water_uptake_max, swing, n_valid_points, MW, nu
- Can be imported into Excel/Python for further analysis

## How to Run

```matlab
cd '/Users/clara/Downloads/RE_ ML for AWH Discussion/temperature'
evaluate_all_salts_atacama
```

**Runtime:** ~2-5 minutes depending on number of salts

## Interpreting Results

### Console Output

The script prints:
```
Rank  Salt               Swing   Min     Max     Points  MW      Nu
----  ----------------  ------  ------  ------  ------  ------  ---
 1    Li+ → Cl-         0.4321  0.1234  0.5555   78/82   42.4   2
 2    Li+ → Br-         0.4102  0.1456  0.5558   76/82   86.9   2
...
```

**Columns:**
- **Rank**: Performance ranking (1 = best)
- **Salt**: Cation-anion pair
- **Swing**: Daily water uptake range (higher = better)
- **Min**: Minimum uptake (day conditions)
- **Max**: Maximum uptake (night conditions)
- **Points**: Valid Atacama points / total (higher = more reliable)
- **MW**: Molecular weight (g/mol)
- **Nu**: Number of ions from dissociation

### Example Interpretation

If LiCl shows:
- Swing = 0.45
- Min = 0.10
- Max = 0.55

**This means:**
- During cold, humid night: 55% water, 45% salt by mass
- During hot, dry day: 10% water, 90% salt by mass
- **Daily harvest potential**: 45% of total mass can be water
- If you have 1 kg of system (salt + water), you can harvest 0.45 kg (450 mL) of water per day!

## Expected Top Performers

Based on literature and typical properties:

### Excellent (Swing > 0.40)
- **LiCl**: Industry standard, very hygroscopic
- **LiBr**: Similar to LiCl, slightly heavier
- **CaCl₂**: Good performance, cheaper than Li salts

### Good (Swing 0.30-0.40)
- **MgCl₂**: Moderate performance, widely available
- **ZnCl₂**: Good capacity, concerns about toxicity

### Moderate (Swing 0.20-0.30)
- **NaCl, KCl**: Common salts, lower performance
- Most other chlorides and bromides

### Poor (Swing < 0.20)
- Sulfates (high deliquescence RH)
- Salts with DRH > 40% (much of Atacama is drier)

## Practical Considerations

The script finds the **thermodynamically optimal** salts, but real-world selection also considers:

### Cost
- **Cheap**: NaCl, KCl, CaCl₂, MgCl₂
- **Moderate**: LiBr
- **Expensive**: LiCl

### Availability
- **Abundant**: Na, K, Ca, Mg salts
- **Limited**: Li salts

### Environmental/Health
- **Safe**: NaCl, KCl, CaCl₂, MgCl₂
- **Concerns**: ZnCl₂ (toxicity), LiCl (corrosive)

### Kinetics
- Script shows **equilibrium** performance
- Real systems need fast absorption/desorption
- Some salts are kinetically slow

## Advanced Analysis

### Swing vs Coverage Tradeoff

Some salts show:
- High swing BUT low coverage = only works in narrow RH range
- Lower swing BUT high coverage = works reliably across conditions

**For Atacama:** High coverage (>80%) is important because conditions vary

### Temperature Dependence

The script uses **temperature-dependent Pitzer parameters**, so it properly accounts for:
- How salt properties change with temperature
- Different behavior at night (cold) vs day (hot)
- More accurate than isothermal models

### Molecular Weight Impact

Lighter salts (low MW) are preferable for mobile systems:
- LiCl (42.4 g/mol) vs CaCl₂ (111 g/mol)
- But LiCl is expensive!
- Plot 2 helps visualize this tradeoff

## Validation

Results should be consistent with:
- Published AWH literature (LiCl typically ranks #1)
- Experimental water uptake data
- Commercial AWH system choices

**If results look wrong:**
- Check solver convergence (look at coverage %)
- Verify Pitzer parameters loaded correctly
- Ensure temperature units are consistent (°C vs K)

## Next Steps

1. **Run the script** to see actual rankings
2. **Compare top 3-5 salts** in detail
3. **Consider cost/availability** for your application
4. **Experimental validation** with top candidates
5. **Kinetic testing** (script only shows equilibrium)

## Related Files

- `licl_temp_RH.m` - Detailed LiCl analysis with 3D plots
- `mgcl2_temp_RH.m` - Detailed MgCl₂ analysis with 3D plots
- `data/Atacama_RH.csv` - Atacama humidity data
- `data/Atacama_Temp.csv` - Atacama temperature data
- `data/parsed_thermodb/pitzer_binary.csv` - Pitzer parameter database

## References

**Atacama Desert:**
- One of driest places on Earth
- High solar irradiance (good for regeneration)
- Large daily T swings (good for AWH cycles)
- Challenging test case

**AWH with Hygroscopic Salts:**
- Kim et al. (2018). Nature Communications
- LaPotin et al. (2021). Joule
- Gao et al. (2021). Science Advances

---

**Questions?** Check the code comments or console output for detailed diagnostic information.
