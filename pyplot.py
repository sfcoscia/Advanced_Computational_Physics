from scipy.io import FortranFile
import numpy as np

# Open the Fortran binary file
with FortranFile('./200/fort.18', 'r') as f:
    # Initialize an empty list to store the data
    data = []
    
    # Loop over the file to read the records
    while True:
        try:
            # Read a record of 6 double precision values (np.float64)
            record = f.read_record(dtype=np.float64)
            
            # Store the record
            data.append(record)
        except EOFError:
            # Break when the end of the file is reached
            break

# Convert the list of records into a single NumPy array
data = np.array(data)

# Print the first 5 records for inspection
print("First 5 records of the data:\n", data[:5])
