import numpy as np

def process_file(file_prefix):
    try:
        # Load the data
        input_file = f'T{file_prefix}.dat'
        output_file_tf = f'TF{file_prefix}.dat'
        output_file_f = f'F{file_prefix}.dat'
        output_file_pf = f'PF{file_prefix}.dat'
        output_file_z = f'Z{file_prefix}.dat'

        # Read input data
        data = np.loadtxt(input_file)
        print(f"Read data from {input_file}")

        # Process to TF file
        col1 = data[:, 0] * 13.6056980659 - 14.2377   #Change according to fermi energy
        col6 = data[:, 5] * 1e-14 #tau
        col11 = data[:, 4] * 1e6
        tf_data = np.column_stack((col1, data[:, 1], data[:, 2], data[:, 3], data[:, 4], col6, data[:, 6], data[:, 7], data[:, 8], data[:, 9], col11))
        np.savetxt(output_file_tf, tf_data)
        print(f"Processed TF data and saved to {output_file_tf}")

        # Process to F file
        col3 = 2.44e-8 * col6 * data[:, 1]
        col4 = col6 * 0.01
        col5 = col11
        f_data = np.column_stack((col1, data[:, 1], col3, col4, col5))
        np.savetxt(output_file_f, f_data)
        print(f"Processed F data and saved to {output_file_f}")

        # Process to PF file
        col2_pf = col5**2 * col4 * 1e-6
        pf_data = np.column_stack((col1, col2_pf))
        np.savetxt(output_file_pf, pf_data)
        print(f"Processed PF data and saved to {output_file_pf}")

        # Process to Z file
        denominators = {
            '300':   4.543753,  #replace with respective Kl values
            '500':  2.726252,
            '700':  1.947323,
            '900':  1.514584,
            '1100': 1.239205,
            '1300': 1.048558
        }
        denom = denominators[file_prefix]
        col3_z = (col5**2 * col4 * data[:, 1] * 1e-10) / (col3 + denom)
        z_data = np.column_stack((col1, data[:, 1], col3_z))
        np.savetxt(output_file_z, z_data)
        print(f"Processed Z data and saved to {output_file_z}")

    except Exception as e:
        print(f"An error occurred while processing {file_prefix}: {e}")

# List of file prefixes to process
file_prefixes = ['300', '500', '700', '900', '1100', '1300']

# Process each file
for prefix in file_prefixes:
    process_file(prefix)

