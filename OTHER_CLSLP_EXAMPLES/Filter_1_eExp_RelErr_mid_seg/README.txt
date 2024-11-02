# Building and Running the CLSLP LPFIR Filter Synthesis Software
#
# -CMAKE is used for builds.  
# -Look at the BUILD_AND_RUN and GENERATE_PLOTS files for examples of the  
#       (LINUX) CMAKE commands for software builds.  
# -Command File GENERATE_PLOTS takes the raw output files creatd by CLSLP and
#       generates GNUPLOT graph files, plots the graphs, and renders these as 
#       PNG graphics files.
#

# Executes the BASH scripts
./BUILD_AND_RUN         # Compiles the C++ sources into an executable CLSLP
./GENERATE_PLOTS        # Executes the software CLSLP that generates an amplitde and phase graphs (png)


