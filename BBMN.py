import pandas as pd
import re
import csv

def process_spectra(target_mz, msp_file, csv_file, mz_tolerance=0.02, intensity_threshold=10000, 
                    output_csv='output_csv.csv', output_msp='output_msp.msp'):
    df = pd.read_csv(csv_file)
    
    matched_compounds = set()
    
    current_compound = None
    current_segment = []
    msp_segments = {}
    found_mz_values = set()
    
    with open(msp_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            if not line:
                if current_compound and current_segment:
                    if len(found_mz_values) == len(target_mz):   
                        matched_compounds.add(current_compound)
                        msp_segments[current_compound] = '\n'.join(current_segment) + '\n\n'
                
                current_segment = []
                found_mz_values = set()
                continue
                
            current_segment.append(line)
            
            if line.startswith('Name:'):
                match = re.search(r'\((.*?)\)', line)
                if match:
                    current_compound = match.group(1)
            
            elif re.match(r'^\d+\.?\d*\s+\d+\.?\d*$', line):
                values = line.split()
                mz = float(values[0])
                intensity = float(values[1])
                
                for target in target_mz:
                    if abs(mz - target) <= mz_tolerance and intensity > intensity_threshold:
                        found_mz_values.add(target)
    
    if current_compound and current_segment and len(found_mz_values) == len(target_mz):
        matched_compounds.add(current_compound)
        msp_segments[current_compound] = '\n'.join(current_segment) + '\n\n'
    
    matched_rows = df[df['Compound'].isin(matched_compounds)]
    
    if not matched_rows.empty:
        matched_rows = matched_rows.fillna('')
        
        with open(output_csv, 'w', newline='') as f:
            f.write(',,,,,,,,,,"Normalised abundance","Raw abundance",,,,,,,,,,\n')
            f.write(',,,,,,,,,,"Condition 1","Condition 1",,,,,,,,,,\n')
            
            writer = csv.writer(f, quoting=csv.QUOTE_ALL)
            writer.writerow(matched_rows.columns)
            writer.writerows(matched_rows.values)
            
        

        with open(output_msp, 'w') as f:
            for compound in matched_rows['Compound']:
                if compound in msp_segments:
                    f.write(msp_segments[compound])
        
        
        for compound in matched_compounds:
            print(compound)
    else:
        print("not matched!")


target_mz = [180.0671] 
msp_file = r'input_msp.msp'
csv_file = r'intput_csv.csv'
mz_tolerance = 0.0005   
intensity_threshold = 6000 

process_spectra(
    target_mz, 
    msp_file, 
    csv_file, 
    mz_tolerance=mz_tolerance,
    intensity_threshold=intensity_threshold,
    output_csv='output_csv.csv',
    output_msp='output_msp.msp'
)