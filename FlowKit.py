# This is a Python script for automation of flow cytometry data analysis with flowkit

# Read in libraries
import bokeh
from bokeh.plotting import show
import matplotlib.pyplot as plt

import flowkit as fk
import glob

bokeh.io.output_notebook()

_ = plt.ioff()

fk.__version__

fcs_path = 'BM_Data_trial/Raw_files.fcs/BM_+9G4_Ig589g_001_024.fcs'
sample = fk.Sample(fcs_path)

#Get metadata
sample.get_metadata()
# Sample channels
sample.channels
# Parameter names
sample.pnn_labels
# Parameter specific names
sample.pns_labels

# Get detectors
detectors = [sample.pnn_labels[i] for i in sample.fluoro_indices]

comp_files = glob.glob('BM_Data_trial/Raw_files.fcs/Compensation*.fcs')
comp_samples = [fk.Sample(fcs_file) for fcs_file in comp_files]



# Plot histogram
p = sample.plot_histogram('FSC-H', source='raw')
show(p)

# Contour plot
f = sample.plot_contour('FSC-H', 'SSC-H', source='raw')
plt.show()

# Contour plot with settings
x_min = y_min = 0
x_max = y_max = 250
f = sample.plot_contour(
  'FSC-H', 'SSC-H', source='raw', x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max
)
plt.show()

# Add in points
f = sample.plot_contour('FSC-H', 'SSC-H', source='raw', plot_events=True, fill=True)
plt.show()

# Interactive scatter plot
p = sample.plot_scatter(
  'FSC-H', 'SSC-H',
  source='raw', y_min=0., y_max=130, x_min=0., x_max=280, color_density=True, bin_width=8
)
show(p)



### Transformation ###
# Default values in flowjo
biex_xform = fk.transforms.WSPBiexTransform('biex', max_value=262144.000029, positive=4.418540, width=-10, negative=0)
sample.apply_transform(biex_xform)

