import streamlit as st
import numpy as np
import altair as alt

import pipefactory as pf

#######################
# Page configuration
st.set_page_config(
    page_title="PipeFactory",
    page_icon="",
    layout="wide",
    initial_sidebar_state="collapsed")


alt.themes.enable("dark")

#######################
# CSS styling
st.markdown("""
<style>

[data-testid="block-container"] {
    padding-left: 2rem;
    padding-right: 2rem;
    padding-top: 1rem;
    padding-bottom: 0rem;
    margin-bottom: -7rem;
}

[data-testid="stVerticalBlock"] {
    padding-left: 0rem;
    padding-right: 0rem;
}

[data-testid="stMetric"] {
    background-color: #393939;
    text-align: center;
    padding: 15px 0;
}

[data-testid="stMetricLabel"] {
  display: flex;
  justify-content: center;
  align-items: center;
}

[data-testid="stMetricDeltaIcon-Up"] {
    position: relative;
    left: 38%;
    -webkit-transform: translateX(-50%);
    -ms-transform: translateX(-50%);
    transform: translateX(-50%);
}

[data-testid="stMetricDeltaIcon-Down"] {
    position: relative;
    left: 38%;
    -webkit-transform: translateX(-50%);
    -ms-transform: translateX(-50%);
    transform: translateX(-50%);
}

</style>
""", unsafe_allow_html=True)


with st.sidebar:

    st.title("PipeFactory")


st.markdown("# PipeFactory : Info & Geometry \n Here you will be able to define the geometry of a pipe. This will allow you export a mesh to run another finite element analysis, or run the analysis directly in this web app.")


# Further processing can be done with the campaign_name variable  
cols = st.columns([1, 2, 3])
with cols[0]:
   st.markdown("**Model Name:**")
with cols[1]:
   campaign_name = st.text_input("", key="campaign_name", label_visibility="collapsed")

cols = st.columns([1, 1, 4])
with cols[0]:
   st.markdown("**Number of Sections:**")
with cols[1]:
   num_sections = st.number_input("", key="num_sections", label_visibility="collapsed", min_value=1, max_value=5, step=1, format="%d")

cols = st.columns([1, 1, 4])
with cols[0]:
   st.markdown("**Element Type:**")
with cols[1]:
   element_type = st.selectbox('',('tri', 'quad', 'hex'),key="element_type", label_visibility="collapsed")


cols = st.columns([1, 1, 4])
with cols[0]:
   st.markdown("**Element Order:**")
with cols[1]:
   element_order = st.number_input("", key="element_order", label_visibility="collapsed", min_value=1, max_value=2, step=1, format="%d")

cols = st.columns([1, 1, 4])
with cols[0]:
   st.markdown("**Element Size:**")
with cols[1]:
   element_size = st.number_input("", key="element_size", label_visibility="collapsed", min_value=0.1, max_value=200.0, value = 30.0, step=0.01)

cols = st.columns(8)
export_types = ['.pipe', '.xml', '.inp', '.vtk']

exports = [ False ] * len(export_types)
exports[0] = True
   
with cols[0]:
   st.markdown('**Export Formats:**')     
for i, t in enumerate(export_types):
   with cols[i + 1]:
      exports[i] = st.checkbox(t, exports[i])



st.markdown("### Define Section Geometry \n Here you need to define how your pipe is made up, whether sections are straight or a bend and the parameters around each.")

type = [ None ] * num_sections
length = [ None ] * num_sections
roc = [ None ] * num_sections
angle = [ None ] * num_sections

cols = st.columns([0.7, 1.5, 0.5, 1.5, 0.5, 1.5,2])
with cols[2]:
   st.markdown("**Radius:**")
with cols[3]:
   radius = st.number_input("", key="pipe_radius", label_visibility="collapsed", min_value=0.0, max_value=200.0, step=0.01)

with cols[0]:
   st.markdown("**Thickness:**")
with cols[1]:
   thickness = st.number_input("", key="pipe_thickness", label_visibility="collapsed", min_value=0.0, max_value=10.0, step=0.01)



for k in range(int(num_sections)):
    cols = st.columns([0.7, 1.5, 0.5, 1.5, 0.5, 1.5,2])
    with cols[0]:
       st.markdown("**Section " + str(k) +"**")

    with cols[1]:
       type[k] = st.selectbox('',('Straight', 'Bend',),key="type_" + str(k), label_visibility="collapsed")

    if type[k] == "Straight":
       with cols[2]:
          st.markdown("**Length**")
       with cols[3]:
          length[k] = st.number_input("", key="length_straight_" + str(k), label_visibility="collapsed", min_value=0.1, max_value= 1e5, value=1000.)
    elif type[k] == "Bend":
       with cols[2]:
          st.markdown("**Angle**")
       with cols[3]:
          angle[k] = st.number_input("", key="angle_bend" + str(k), label_visibility="collapsed", min_value=-180.0, max_value= 180.0, value=0.0)
       with cols[4]:
        st.markdown("**Radius**")
       with cols[5]:
          roc[k] = st.number_input("", key="roc_bend_" + str(k), label_visibility="collapsed", min_value= 10.0, max_value= 10e5, value=1000.)

st.markdown("### Build Model \n It is time to build your pipe model. Here you can save you mesh into different formats so they can be used in other finite element codes or visualised prior to analysis.")

   
cols = st.columns([5, 1])
with cols[0]:
    submitted_model = st.button("Build Model", use_container_width=True)

mesh = None


if submitted_model:        
    st.write(type)
    section_list = []

    for i, t in enumerate(type):
        
        section_list.append({})
        
        if t.lower() == "straight":
            
            section_list[-1]['type'] = "straight"
            section_list[-1]['length'] = length[i]
            
            
        elif t.lower() == "bend":
            section_list[-1]['type'] = "bend"
            section_list[-1]['param'] = {}
            section_list[-1]['param']['radius'] = roc[i]
            section_list[-1]['param']['angle'] = angle[i]
    
    higher_order = False
    if element_order == 2:
            higher_order = True
    mesh = pf.Pipe(radius,thickness,section_list,(element_type,higher_order), element_size)
    mesh.export(campaign_name + ".vtk")



   

    
      
   





