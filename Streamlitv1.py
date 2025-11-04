import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import statsmodels 
import matplotlib.cm as cm #colormap
from matplotlib.lines import Line2D
import matplotlib.colors as colors
import seaborn as sns
import seaborn.objects as so
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import matplotlib.ticker as tkr
from adjustText import adjust_text
import matplotlib.patches as mpl_patches
from matplotlib import gridspec
from sklearn import linear_model
import plotly.express as px
from sklearn import manifold
from sklearn.decomposition import PCA
from sklearn.metrics import euclidean_distances 
from sklearn import preprocessing
import umap
import streamlit as st
from IPython.utils import io

st.set_page_config(layout='wide')


mpl.rcParams['font.sans-serif'] = "Arial"
mpl.rcParams['font.weight'] = "normal"
#mpl.rc('font',**{'family':'sans-serif', 'sans-serif':['Helvetica']})
#rc('font',**{'family':'serif','serif':['Courier New'], 'weight':'bold'})
#rc('text', usetex=False)
plt.rcParams['font.size'] = 7
#plt.rcParams["font.weight"] = 'bold'
#plt.rcParams["axes.labelweight"] = 'bold'
mpl.rcParams['axes.linewidth'] = 1
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "Arial"
})

# set tick width
mpl.rcParams['xtick.major.size'] = 4
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['ytick.major.size'] = 4
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams["legend.frameon"] = False
mpl.rcParams["axes.axisbelow"] = True

taxa_filled_02 = pd.read_csv('GP15-17-OCE-02-3um.csv')
taxa_filled_3 = pd.read_csv('GP15-17-OCE-3-51um.csv')

sizefract = st.sidebar.radio('Size Fraction:',['0.2 - 3 µm', '3 - 51 µm'],index = 0)

if sizefract =='0.2 - 3 µm':
    taxa_domain_summed = pd.read_csv('Int_files_Streamlit/domain_02.csv')
    taxa_clade1_summed = pd.read_csv('Int_files_Streamlit/clade1_02.csv')
    taxa_clade2_summed = pd.read_csv('Int_files_Streamlit/clade2_02.csv')
    taxa_clade3_summed = pd.read_csv('Int_files_Streamlit/clade3_02.csv')
    taxa_clade4_summed = pd.read_csv('Int_files_Streamlit/clade4_02.csv')
    taxa_clade5_summed = pd.read_csv('Int_files_Streamlit/clade5_02.csv')
    taxa_genus_summed = pd.read_csv('Int_files_Streamlit/genus_02.csv')
    taxa_species_summed = pd.read_csv('Int_files_Streamlit/species_02.csv')
    taxa_kegg_summed = pd.read_csv('Int_files_Streamlit/kegg_02.csv')
    taxa_cog_summed = pd.read_csv('Int_files_Streamlit/cog_02.csv')
    taxa_protname_summed = pd.read_csv('Int_files_Streamlit/protname_02.csv')
    taxa_path_summed = pd.read_csv('Int_files_Streamlit/path_02.csv')
    taxa_EC_summed = pd.read_csv('Int_files_Streamlit/EC_02.csv')
    taxa_PFAM_summed = pd.read_csv('Int_files_Streamlit/PFAM_02.csv')
    taxa_module_summed = pd.read_csv('Int_files_Streamlit/module_02.csv')
#    taxa_GO_summed = pd.read_csv('Int_files_Streamlit/GO_02.csv')
    taxa_TC_summed = pd.read_csv('Int_files_Streamlit/TC_02.csv')
    taxa_filled = taxa_filled_02
    stn = pd.read_csv('GP15-17-OCE_stns.txt', sep='\t', encoding='latin-1')
    stn_names = pd.read_csv(r'Stn_Names_comma.csv')
    stn_ID = stn_names.Stn.astype(str)
    stn_ID = np.array(stn_ID)
    stn_keys = [1,3,4,5,6,8,9,10,12,14,15,16,18,18.3,19,21,23,25,27,29,31,33,35,37,39,1,3,6,8,10,12,14,16,18,20,22,25,27,29,32,33,35,37,38]
    stn_vals = ['GP15_02_1', 'GP15_02_3', 'GP15_02_4', 'GP15_02_5', 'GP15_02_6',
       'GP15_02_8', 'GP15_02_9', 'GP15_02_10', 'GP15_02_12', 'GP15_02_14',
       'GP15_02_15', 'GP15_02_16', 'GP15_02_18', 'GP15_02_18-3',
       'GP15_02_19', 'GP15_02_21', 'GP15_02_23', 'GP15_02_25',
       'GP15_02_27', 'GP15_02_29', 'GP15_02_31', 'GP15_02_33',
       'GP15_02_35', 'GP15_02_37', 'GP15_02_39', 'GP17-OCE_02_1',
       'GP17-OCE_02_3', 'GP17-OCE_02_6', 'GP17-OCE_02_8',
       'GP17-OCE_02_10', 'GP17-OCE_02_12', 'GP17-OCE_02_14',
       'GP17-OCE_02_16', 'GP17-OCE_02_18', 'GP17-OCE_02_20',
       'GP17-OCE_02_22', 'GP17-OCE_02_25', 'GP17-OCE_02_27',
       'GP17-OCE_02_29', 'GP17-OCE_02_32', 'GP17-OCE_02_33',
       'GP17-OCE_02_35', 'GP17-OCE_02_37', 'GP17-OCE_02_38']
    stn_dict = dict(zip(stn_keys, stn_vals))
    

if sizefract =='3 - 51 µm':
    taxa_domain_summed = pd.read_csv('Int_files_Streamlit/domain_3.csv')
    taxa_clade1_summed = pd.read_csv('Int_files_Streamlit/clade1_3.csv')
    taxa_clade2_summed = pd.read_csv('Int_files_Streamlit/clade2_3.csv')
    taxa_clade3_summed = pd.read_csv('Int_files_Streamlit/clade3_3.csv')
    taxa_clade4_summed = pd.read_csv('Int_files_Streamlit/clade4_3.csv')
    taxa_clade5_summed = pd.read_csv('Int_files_Streamlit/clade5_3.csv')
    taxa_genus_summed = pd.read_csv('Int_files_Streamlit/genus_3.csv')
    taxa_species_summed = pd.read_csv('Int_files_Streamlit/species_3.csv')
    taxa_kegg_summed = pd.read_csv('Int_files_Streamlit/kegg_3.csv')
    taxa_cog_summed = pd.read_csv('Int_files_Streamlit/cog_3.csv')
    taxa_protname_summed = pd.read_csv('Int_files_Streamlit/protname_3.csv')
    taxa_path_summed = pd.read_csv('Int_files_Streamlit/path_3.csv')
    taxa_EC_summed = pd.read_csv('Int_files_Streamlit/EC_3.csv')
    taxa_PFAM_summed = pd.read_csv('Int_files_Streamlit/PFAM_3.csv')
    taxa_module_summed = pd.read_csv('Int_files_Streamlit/module_3.csv')
#    taxa_GO_summed = pd.read_csv('Int_files_Streamlit/GO_3.csv')
    taxa_TC_summed = pd.read_csv('Int_files_Streamlit/TC_3.csv')
    taxa_filled = taxa_filled_3  
    stn = pd.read_csv('GP15-17-OCE_stns3.txt', sep='\t', encoding='latin-1')
    stn3 = pd.read_csv('stnlabel_to_ID_3.txt', sep='\t', encoding='latin-1')
    stn_ID = stn3['Stn ID'].astype(str)
    stn_ID = np.array(stn_ID)
    stn_keys=[3,4,5,6,8,9,10,12,14,15,16,18,18.3,19,21,23,25,27,29,31,33,35,37,39,1,3,6,8,10,12,14,16,18,20,22,25,27,32,33,35,37,38]
    stn_vals = ['GP15_3_3',
     'GP15_3_4',
     'GP15_3_5',
     'GP15_3_6',
     'GP15_3_8',
     'GP15_3_9',
     'GP15_3_10',
     'GP15_3_12',
     'GP15_3_14',
     'GP15_3_15',
     'GP15_3_16',
     'GP15_3_18',
     'GP15_3_18-3',
     'GP15_3_19',
     'GP15_3_21',
     'GP15_3_23',
     'GP15_3_25',
     'GP15_3_27',
     'GP15_3_29',
     'GP15_3_31',
     'GP15_3_33',
     'GP15_3_35',
     'GP15_3_37',
     'GP15_3_39',
     'GP17-OCE_3_1',
     'GP17-OCE_3_3',
     'GP17-OCE_3_6',
     'GP17-OCE_3_8',
     'GP17-OCE_3_10',
     'GP17-OCE_3_12',
     'GP17-OCE_3_14',
     'GP17-OCE_3_16',
     'GP17-OCE_3_18',
     'GP17-OCE_3_20',
     'GP17-OCE_3_22',
     'GP17-OCE_3_25',
     'GP17-OCE_3_27',
     'GP17-OCE_3_32',
     'GP17-OCE_3_33',
     'GP17-OCE_3_35',
     'GP17-OCE_3_37',
     'GP17-OCE_3_38']
    stn_dict3 = dict(zip(stn_keys, stn_vals))

################
optiontax = st.selectbox('Taxonomic Resolution:',['Domain', 'Supergroup', 'Phylum', 'Class','Order','Family','Genus','Species'],index = 0)

if optiontax =='Domain':
    data = taxa_domain_summed
if optiontax =='Supergroup':
    data = taxa_clade1_summed
if optiontax =='Phylum':
    data = taxa_clade2_summed
if optiontax =='Class':
    data = taxa_clade3_summed
if optiontax =='Order':
    data = taxa_clade4_summed
if optiontax =='Family':
    data = taxa_clade5_summed
if optiontax =='Genus':
    data = taxa_genus_summed
if optiontax =='Species':
    data = taxa_species_summed

df = pd.DataFrame(data)

df_sort = df.sort_values('sum')

# Create stacked bar plot
fig = px.bar(df, x='stn', y='sum', color='clade', barmode='stack',template="plotly_white", width=1500, height=800,labels = {'stn':'Station','sum':'Fractional Sum', 'clade':'Clade'}) #, color_discrete_map= clade4_dict
fig.update_traces(marker_line_width=0)
fig.update_xaxes(showline=True, linewidth=1, linecolor='black', mirror=True, tickmode = "array", tickvals = stn_vals, ticktext = stn_keys)
fig.update_yaxes(showline=True, linewidth=1, linecolor='black', mirror=True, range = [0,1])
fig.update_layout(
    font_family="Arial", font_size = 14, autosize=True
)
#fig.update_layout(legend=dict(
#    orientation="h",
#    entrywidth=70,
#    yanchor="bottom",
#    y=1.02,
#    xanchor="right",
#    x=1
#))

st.plotly_chart(fig,use_container_width=True)

#############

taxselect = st.selectbox('Taxon:',list(data.clade.unique()))


data1 = data[data['clade'].str.contains(taxselect)]

data2 = data1.groupby("stn").agg(
    lat = pd.NamedAgg(column="lat", aggfunc="min"),
    lon = pd.NamedAgg(column="lon", aggfunc="min"),
    summed = pd.NamedAgg(column="sum", aggfunc="sum"),
    stn = pd.NamedAgg(column="stn", aggfunc="min"),
    param_group = pd.NamedAgg(column="clade", aggfunc="sum")
)


fig = px.scatter_map(data2, lat = 'lat', lon = 'lon', color = 'summed',
                     hover_name="stn", size="summed",size_max = 30,color_continuous_scale="YlOrRd",opacity = 0.7)

fig.update_layout(height=1000, width = 700)
fig.update_layout(coloraxis_colorbar_title_text = '% per station')
fig.update_layout(
    autosize=False,
    hovermode='closest',
    map=dict(
        bearing=0,
        pitch=0,
        zoom=1,
        bounds=dict(
            west=-170,
            south=-70,
            east=-49,
            north=63),
        style="satellite"
    ))

fig.update_layout(
    font_family="Arial", font_size = 14, font_color = 'black'
    
)

st.plotly_chart(fig,use_container_width=False)

##############


markersinfo = '''
## Protein Markers
### Iron
| Protein Name | Function |
| :- | :- |
FecA	| Fe(III) dicitrate transport protein
FecB | Fe complex transport system SBP
FutA1/AfuA	| Fe(III) transport system SBP
IrpA/IcmP	| Putative Fe regulated protein
Ftn	| Ferritin - Fe storage protein
IsiB/NifF/FldA |	Flavodoxin - electron transport protein
EfeO	| Fe uptake system component
EfeU	| High affinity Fe transporter
Acr	| Fe complex outer membrane receptor
Sir | Ferredoxin-sulfite reductase
Ho1, PbsA1, HmuO |	Heme oxygenase
Fiu	| Catecholate siderophore receptor
SfuA	| Fe(III) transport system SBP
BhuR/HemR/TbpA	| Hemoglobin/transferrin/lactoferrin receptor

### Zinc
| Protein Name | Function |
| :- | :- |
FtsH	| Cell division metalloprotease
lraI	| Zn/Mn transport system SBP
COG0523 Cluster 1	| Putative Zn chaperone
ZnuA	| Zn transport system SBP
ZnuD	| Zn transport system ATP-binding protein

### Vitamin B12
| Protein Name | Function |
| :- | :- |
NrdJ	| Ribonucleotide reductase, class II (B12-dependent)
NrdA/Z	| Ribonucleotide reductase, class Ia (B12-independent)
NrdB	| Ribonucleotide reductase, class Ia (B12-independent)
BtuB	| B12 outer membrane transporter
MetH	| Methionine synthase (B12-dependent)
MetE	| Methionine synthase (B12-independent)
CUBN	| B12 outer membrane receptor
CobA	| B12 synthesis protein (cob(I)alamin adenosyltransferase)
CobS	| B12 synthesis protein (cobaltochelatase)

### Nitrate
| Protein Name | Function |
| :- | :- |
UrtA	| Urea transport system SBP
UrtD	| Urea transport system ATP-binding protein
UrtE	| Urea transport system ATP-binding protein
UreA	| Urease
UreC	| Urease
NtrP	| Nitrate/Nitrite transporter
NtcA	| Global nitrogen transcriptional regulator

### Phosphate
| Protein Name | Function |
| :- | :- |
PstB	| Phosphate transport system ATP-binding protein
PstS	| Phosphate transport system SBP
PhnD	| Phosphonate transport system SBP
SqdB	| Sulfolipid synthesis protein
PTC1	| Protein phosphatase
Acid Phosphatase	| Acid phosphatase
'''
st.markdown(markersinfo)

####

protselect = st.selectbox('Select a protein or enter your own:',['fecA', 'cobW','irpA','metE','metH','ureC','urtA','ftsH','btuB'],accept_new_options=True, index = 0)

data = taxa_protname_summed#.groupby(taxa_protname_summed['protname'])

#data = data.get_group('fecA')

data1 = data[data['protname'].str.contains(protselect)]

data2 = data1.groupby("stn").agg(
    lat = pd.NamedAgg(column="lat", aggfunc="min"),
    lon = pd.NamedAgg(column="lon", aggfunc="min"),
    summed = pd.NamedAgg(column="sum", aggfunc="sum"),
    stn = pd.NamedAgg(column="stn", aggfunc="min"),
    param_group = pd.NamedAgg(column="protname", aggfunc="sum")
)


fig = px.scatter_map(data2, lat = 'lat', lon = 'lon', color = 'summed',
                     hover_name="stn", size="summed",size_max = 30,color_continuous_scale="YlOrRd",opacity = 0.7)

#fig = px.scatter_mapbox(data, lat = 'lat', lon = 'lon', color = "sum",
#                     hover_name="stn", size="sum",size_max = 30,color_continuous_scale="YlOrRd", opacity = 1)
#fig.update_geos(fitbounds="locations")

#zoom 4, height 4000, width 1400 for print, size max 40

fig.update_layout(height=1000, width = 700)
fig.update_layout(coloraxis_colorbar_title_text = '% per station')
fig.update_layout(
    autosize=True,
    hovermode='closest',
    map=dict(
        bearing=0,
        pitch=0,
        zoom=1,
        bounds=dict(
        west=-170,
        south=-70,
        east=-49,
        north=63),
        style="satellite",
        
    ))
fig.update_layout(
    font_family="Arial", font_size = 14, font_color = 'black'
    
)

st.plotly_chart(fig,use_container_width=False)

##########################################


taxa_filled_small = taxa_filled[['Protein', 'KEGG_ko','EC','KEGG_TC','COG_category','PFAMs', 'KEGG_Pathway','KEGG_Module','GOs']] 

taxa_filled_small['protname'] = taxa_filled['Preferred_name']
taxa_filled_small['Description'] = taxa_filled['Description_y']
taxa_filled_small[['Domain', 'Supergroup', 'Phylum', 'Class','Order','Family','Genus','Species']] = taxa_filled[['domain','clade1','clade2','clade3','clade4','clade5','genus','species']]
taxa_filled_small = taxa_filled_small.fillna('Other')
taxa_filled_small[stn_ID] = taxa_filled[stn_ID]


toplot = taxa_filled_small[taxa_filled_small['protname'].str.contains(protselect)] #ko:K16087, ko:K16091, ko:K09815,ko:K02077,ko:K11959
maxval = toplot[stn_ID].max().max()
xlabels = stn_keys
fig = plt.figure(figsize = (8, 3))
fig.patch.set_alpha(0.0)
cmap = ListedColormap(['#01ff07','#fe01b1','#ff9408','#a9561e','#490648','#aa23ff','#13bbaf','#247afd','#ec2d01','#3f9b0b', 'gray','tan','orange'])

ax = plt.subplot()
ax.patch.set_alpha(0.0)
ax.set_title(str(protselect), color = "#e8e8e8", size = 9, loc = 'left')
ax = pd.plotting.parallel_coordinates(toplot, cols = stn_ID, class_column = optiontax, colormap = cmap, axvlines = False, linewidth = 0.5)#color=toplot.colours)
ax.set_ylim([0,maxval])
ax.spines['top'].set_color("orange")
ax.spines['bottom'].set_color("orange")
ax.spines['left'].set_color("orange")
ax.spines['right'].set_color("orange")
ax.tick_params(axis='x', colors="#e8e8e8", labelsize = 7, size = 0)
ax.tick_params(axis='y', colors="orange", labelsize = 7, size = 0)
ax.grid(False)
ax.vlines(range(0,len(xlabels)),ymin = 0, ymax = maxval,color = "#e8e8e8", zorder = 0, alpha = 0.5, linewidth = 0.3)
with io.capture_output() as captured:
    ax.set_xticklabels(xlabels)
ax.get_legend().remove()
fig.legend(loc='outside right upper', frameon = True, framealpha = 1, fontsize = 6, facecolor='k', labelcolor ='#e8e8e8', draggable = True, edgecolor = 'orange')
#plt.text(0.1,0.1,str(protselect),size = 15, weight ='bold')
ax.ticklabel_format(scilimits=(0,0), axis = 'y')

st.pyplot(fig, width = 'stretch')
#fig.savefig('lines.png', format="png", dpi=1000, transparent=True)
#st.image('lines.png', width = 'content')


##########################################



