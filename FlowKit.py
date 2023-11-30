# This is a Python script for automation of flow cytometry data analysis with flowkit

# Read in libraries
import bokeh
import pandas as pd
from bokeh.plotting import show
import matplotlib.pyplot as plt
import os

import flowkit as fk
import glob

bokeh.io.output_notebook()

_ = plt.ioff()

fk.__version__

fcs_path = 'BM_Data_trial/Raw_files.fcs/BM_+9G4_Ig589g_001_024.fcs'

sample = fk.Sample(fcs_path)

sample_path = 'BM_Data_trial/Raw_files.fcs'
wsp_path = 'BM_Data_trial/20210923_BM_trial.wsp'
workspace = fk.Workspace(wsp_path, fcs_samples=sample_path)

workspace.get_sample_ids()

sample_id = 'BM_+9G4_Ig589g_001_024.fcs'
sample = workspace.get_sample(sample_id)
sample.channels
print(workspace.get_gate_hierarchy(sample_id, 'ascii'))

workspace.analyze_samples(sample_id=sample_id)
results = workspace.get_gating_results(sample_id)
results.report.head()


# Create a new session
full_sample_path = os.path.join(sample_path, sample_id)
session = fk.Session(fcs_samples=full_sample_path)
sample = session.get_sample(sample_id)

# Load compensation matrix
detectors = [sample.pnn_labels[i] for i in sample.fluoro_indices]
den_comp_mat = fk.Matrix('den_comp', sample.metadata['spill'], detectors)
session.add_comp_matrix(den_comp_mat)

# Data Transformation
# Note: these are the default values in FlowJo (and in FlowKit, so you don't have to remember these specific values)
biex_xform = fk.transforms.WSPBiexTransform('biex', max_value=262144.000029, positive=4.418540, width=-10, negative=0)
session.add_transform(biex_xform)

# Gate: Lymphocytes PolygonGate
FSC_dim = fk.Dimension('FSC-A', compensation_ref='uncompensated', transformation_ref=biex_xform.id, range_min=20000, range_max=100000)
SSC_dim = fk.Dimension('SSC-A', compensation_ref='uncompensated', transformation_ref=biex_xform.id, range_min=0, range_max=70000)
lymphocyte_gate = fk.gates.PolygonGate(
    'Lymphocytes',
    dimensions=[FSC_dim, SSC_dim],
    vertices=[[20000, 70000], [100000, 70000], [100000, 0], [20000, 0]]
)
session.add_gate(lymphocyte_gate, gate_path=('root',))

session.analyze_samples(sample_id=sample_id)
p = session.plot_gate(sample_id, "Lymphocytes")
show(p)


# remove gate 
# session.remove_gate('Lymphocytes', gate_path=('root',))
# Gate: Single cells PolygonGate
FSC_A_dim = fk.Dimension('FSC-A', 'uncompensated', transformation_ref=biex_xform.id)
FSC_H_dim = fk.Dimension('FSC-H', 'uncompensated', transformation_ref=biex_xform.id)

# Gate: Single Cells PolygonGate

# Gate: Live RectangleGate

# Gate: 9G4+ B cells RectangleGate

# Gate: B cells RectangleGate

# Gate: 9G4+ Mature_CD93-CD23+ RectangleGate

# Gate: 9G4+ T1_CD93+CD23- RectangleGate

# Gate: Mature_CD93-CD23+ RectangleGate

# Gate: T1_CD93+CD23- RectangleGate



print(g_strat.get_gate_hierarchy(output='ascii'))




# Plot transformed data
p = sample.plot_scatter(1, 4, source='raw')
show(p)






#Get metadata
sample.get_metadata()
# Sample channels
sample.channels
# Parameter names
sample.pnn_labels
# Parameter specific names
sample.pns_labels

# Data Transformation
# Note: these are the default values in FlowJo (and in FlowKit, so you don't have to remember these specific values)
biex_xform = fk.transforms.WSPBiexTransform('biex', max_value=262144.000029, positive=4.418540, width=-10, negative=0)
sample.apply_transform(biex_xform)
# Plot transformed data
p = sample.plot_scatter(9, 10, source='xform')
show(p)

# Apply compensation matrix from metadata 'spill'
sample.apply_compensation(sample.metadata['spill'])
# Plot compensated data
fig = sample.plot_scatter(1, 4, source='xform', x_axis_type="scientific", y_axis_type="scientific")
show(fig)


# Gating strategy
import numpy as np

fsc_data = sample.get_channel_events(1)
ssc_data = sample.get_channel_events(4)

# Create dimensions for FSC and SSC
center = (np.mean(fsc_data), np.mean(ssc_data))
var_fsc = np.var(fsc_data)
var_ssc = np.var(ssc_data)
cov = [[var_fsc, 0], [0, var_ssc]]
dist = 2  # containing about 95% of the data
dim_fsc = fk.Dimension('FSC-A', range_min=20000, range_max=100000)
dim_ssc = fk.Dimension('SSC-A', range_min=0, range_max=70000)

# Create a gate for lymphocytes based on FSC and SSC
lymph_gate = fk.gates.EllipsoidGate('Lymphocytes', dimensions=[dim_fsc, dim_ssc], coordinates=center, covariance_matrix=cov, distance_square=dist)

# Create a GatingStrategy instance
g_strat = fk.GatingStrategy()

# Add the lymphocyte gate to the gating strategy
g_strat.add_gate(lymph_gate, gate_path=('root',))

res = g_strat.gate_sample(sample)
res.report






# Apply gating: 9G4+ Mature_CD93-CD23+

dim_a = fk.Dimension('R670-A', range_max=50000)
dim_b = fk.Dimension('YG586-A', range_min=50000)
rect_top_left_gate = fk.gates.RectangleGate('top-left', dimensions=[dim_a, dim_b])
g_strat = fk.GatingStrategy()
g_strat.add_gate(rect_top_left_gate, gate_path=('root',))
res = g_strat.gate_sample(sample)
res.report

events_9G4 = sample.get_channel_events(8, 'raw')
events_9G4.min(), events_9G4.max()
fk.GatingStrategy()
# 9G4+ gate
gate_9g4 = fk.gates.ThresholdGate(2000, 9, 'above')  # Adjust threshold based on your data
# CD93- gate
gate_cd93 = fk.gates.ThresholdGate(1000, 10, 'below')  # Adjust threshold based on your data
# CD23+ gate
gate_cd23 = fk.gates.ThresholdGate(2000,14, 'above')  # Adjust threshold based on your data

