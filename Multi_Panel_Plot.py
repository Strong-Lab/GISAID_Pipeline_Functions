import Plotting_Functions as plotf
import pandas as pd
import numpy as np
import os
import math
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from datetime import datetime
from datetime import date

def axes_line_plot(TS,
                   ax,
              graph_dates,
              protein=None,
              continent="Global",
              complex_name=None,
              color_by_variant=None,
              color_list=["firebrick","#1E399A","purple","forestgreen","lightcoral","chocolate","gold","darkcyan","black","#AAAAFF"],
              marker_list=["o","D","^","s"],
              axes_bounds="default",
              x_tick_properties=None,
              y_tick_properties=None,
              title_properties=None,
              legend=False, 
              legend_bbox_position=(1.02,0.9),
              fontsize_x_label=12,
              fontsize_y_label=12):
    """
    Creates a line plot from time series data

    Arguments
    ----------
    TS: A time series matrix.
    
    ax: A matplotlib axes object within which to draw the figure.

    Graph_dates: A list of dates that will be displayed on the x-axis tick labels. Can be in any format but must correspond to the dates in the TS matrix to properly represent data.
    
    Protein (string): the name of the protein for printing the title. 
    
    Complex_name: If applicable and specified, print the title of the complex (i.e. RdRP) instead of the protein name. Highly reccomended for aggregate analyses.
    
    color_by_variant (optional, default=None): Can pass a dictionary with variant names as keys and a dictionary of matplotlib line2D properties (such as marker size, marker color, etc.) to plot each unique variant with a desired style.
    
    color_list (list): If color_by_variant=None, a list of colors may be passed here. Plotting will cycle through the marker list, moving to the next color each time a new line is created from the data.
    
    marker_list (list): If color_by_variant=None, a list of markers may be passed here. Plotting will cycle through the marker list, moving to the next marker each time a new line is created from the data.
    
    figsize (tuple): Set the size of the figure in inchebs (width,height).
    
    axes_bounds: May be "default" or "RdRP". Setting the window to "RdRP" changes the width of the window to accomodate longer domain labels in the legend for the polymerase complex.
    
    center_title_to_figure (default=True): If false, the title is centered to the plotting window instead of the figure (the figure includes the plotting window and the legend)
    
    fontsize_x_label (int, default=12): Defines the font size of the x-axis label in points.
    
    fontsize_y_label (int, default=12): Defines the font size of the y-axis label in points.
    
    fontsize_x_ticks (int, default=10): Defines the font size of the x-axis tick labels in points.
     
    fontsize_y_ticks (int, default=12): Defines the font size of the y-axis tick labels in points.
    
    x_label_rotation (int, default=45): Degrees of counterclockwise rotation from horizontal position on axis (90=vertical labels, 0=horizontal labels).
    """
    import matplotlib.ticker as mtick
    
    ##### Argument Set-up #####
    #Fill default values of x_tick_properties if they are not specified
    x_tick_properties=plotf.fill_defaults(x_tick_properties,{'labelsize':10,'rotation':45})
               
    #Check and fill y-axis property kwargs
    y_tick_properties=plotf.fill_defaults(y_tick_properties,{'labelsize':12})
    
    title_properties=plotf.fill_defaults(title_properties,{'fontsize':26,'fontweight':'medium','pad':20})
    #####
    
    #Line Graph for Variant Prevalence Over Time
    x=graph_dates
    labels=list(TS.index)

    #Axes bounds relative to the figure vary based on size of legend. Special axes bounds used for the RdRP are saved here.
    #Axes bounds can be chosen with the string "axes_ratio"
    if axes_bounds=="default":
        ax_bound=[0.1,0.2,0.7,0.7]
    elif axes_bounds=="RdRP":
        ax_bound=[0.1,0.2,0.6,0.7]
    
    #Store the records of each line object created in the plot
    lines=[]
    #Standard: plot each variant with default color scheme unless color_by_variant is specified
    if color_by_variant==None:
        for i in range(0,len(TS),1):
            line_obj=ax.plot(x,TS.iloc[i,:],ls='-',marker=marker_list[i%len(marker_list)],ms=7,mew=1,mec="#FFFFFFEE",label=labels[i],color=color_list[i%len(color_list)])
            lines.append(line_obj[0])
    #If color_by_variant is specified, pass the kwargs corresponding to the current variant name to the plot function
    else:
        for i in range(0,len(TS),1):
            line_obj=ax.plot(x,TS.iloc[i,:],ls='-',label=labels[i],**color_by_variant[labels[i]])
            lines.append(line_obj[0])
            
    #Legend
    if legend==True:
        leg=ax.legend(bbox_to_anchor=legend_bbox_position,fontsize=14,loc="upper left",markerscale=1.5)
        leg.set_title("Mutation",prop={'weight':'medium','size':18})

    #X-tick properties (ticks and label placement)
    ax.tick_params(which="major",axis='x',width=1.25,**x_tick_properties)
    #X-tick text properties
    plt.setp(ax.get_xticklabels(), ha="right",rotation_mode="anchor")

    #Y-axis tick/label parameters
    ax.tick_params(which="major",axis='y',width=1.25,**y_tick_properties)

    #Adjust y-axis ticks
    ax.set_yticks(np.arange(0,1.1,0.1))
    ax.set_yticks(np.arange(0,1.1,0.02),minor=True);

    #Upper bound is the number of columns in the dataframe minus 1.0 (indexing starts at zero), plus 0.5 for space
    ax.set_xlim(left=-0.5,right=TS.shape[1]-0.5)
    ax.set_ylim(-0.02,1.01)

    #Y-axis tick labels as percentage
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1,decimals=0))
    
    #Set axes labels
    ax.set_ylabel("Prevalence of Mutation",fontsize=fontsize_y_label)
    ax.set_xlabel("Collection Date of Sample",fontsize=fontsize_x_label)

    #Define figure title
    if protein and continent=="Global":
        title="Most Common Mutations for {}: Global Prevalence".format(protein)
    elif protein and continent!="Global":
        title="Most Common Mutations for {}: {}".format(protein,continent)
    elif complex_name and continent=="Global":
        title="Most Common Mutations for {}: Global Prevalence".format(complex_name)
    elif complex_name and continent!="Global":
        title="Most Common Mutations for {}: {}".format(complex_name,continent)
    #Warn and print no title for invalid combinations
    elif (complex_name and protein) or (not complex_name and not protein):
        warnings.warn("Invalid combination of input for 'protein' and 'complex_name' arguments. Please specify input for either 'protein' or 'complex_name' (not both or neither).",category=UserWarning)
    
    #Center title to axis
    ax.set_title(title,**title_properties)
    
    #Create a grid
    ax.grid(which='major',axis='both',color="#DDDDDD",alpha=0.3)

    return lines

def axes_plot_n_seq(n_seq_path,ax,graph_dates,color_dict={'Europe':"firebrick",'North America':"#1E399A",'Oceania':"purple", 'Asia':"forestgreen",'South America':"chocolate","Africa":"darkcyan"},marker_list=["o","D","s"],remove_end=0,legend=False,legend_bbox_position=(1.02,0.9),log=False,figsize=(14,8),x_tick_properties=None,y_tick_properties=None,title_properties=None,line_kwargs=None,y_max=None,y_tick_freq=(500,100)):
    """
    Creates a plot for number of sequences by continent within the subplot it is called within.
    
    Arguments
    ----------
    
    y_max (int): May specify the upper bound of the y-limits using this argument.
    
    y_tick_freq (tuple): The first element of the tuple specifies the frequency of major tick marks with labels, and the second element specifes the frequency of minor tick marks. Default is (500,100).
    """
    ##### Argument Set-up #####
    #Fill default values of x_tick_properties if they are not specified
    x_tick_properties=plotf.fill_defaults(x_tick_properties,{'labelsize':10,'rotation':45})
               
    #Check and fill y-axis property kwargs
    y_tick_properties=plotf.fill_defaults(y_tick_properties,{'labelsize':12})
    
    #Check and fill title property kwargs
    title_properties=plotf.fill_defaults(title_properties,{'fontsize':26,'fontweight':'medium','pad':20})
    
    #Check for incompatible arguments log=True and user defined y-max
    if y_max and log==True:
        raise ValueError("Cannot specify argument y_max when log=True")
    #####           
       
    #Extract protein name from path: by default, protein name is after the last forward slash and before the first underscore
    protein=n_seq_path.split("/")[-1]
    protein=protein.split("_")[0]
    
    #Load total sequences by continent file
    if remove_end==0:
        n_seq=pd.read_csv(n_seq_path)
    elif remove_end>0:
        n_seq=pd.read_csv(n_seq_path).iloc[:,:-remove_end]
    elif remove_end<0:
        raise ValueError("Negative value specified for remove_end.")
    
    #Rename first column to "Continent" and make it the index
    newcols=list(n_seq.columns)
    newcols[0]="Continent"
    n_seq.columns=newcols
    n_seq.index=n_seq["Continent"]
    #Remove 'continent' column
    n_seq=n_seq.drop("Continent",axis=1)
    #Remove 'worldwide' row
    n_seq=n_seq.drop("Worldwide",axis=0)

    #Sort based on average number of sequences
    n_seq["avg"]=n_seq.mean(axis=1)
    n_seq=n_seq.sort_values(by="avg",ascending=False)
    n_seq=n_seq.drop("avg",axis=1)
    
    #Create line plot of total number of sequences by continent

    #Line Graph for Variant Prevalence Over Time
    x=graph_dates
    labels=list(n_seq.index)
    
    lines=[]
    #Create lines for each of the six continents
    for i in range(0,len(n_seq),1):
        line_obj=ax.plot(x,n_seq.iloc[i,:],ls='-',label=labels[i],color=color_dict[labels[i]],marker=marker_list[i%len(marker_list)],**line_kwargs)
        lines.append(line_obj[0])
    
    if legend==True:
        leg=ax.legend(bbox_to_anchor=legend_bbox_position,fontsize=14,loc="upper left",markerscale=2 if line_obj[0].get_markersize()<=6 else 1.5)
        leg.set_title("Continent",prop={'weight':'medium','size':18})
    
    #Set axes limits
    #Y limits: automatically calculate if not defined by the user, and do not calculate for log scale
    if log==False and not y_max:
        #Automatic calculation: if the number of sequences in any given week is greater than 7500, adjust the limits of the plot accordingly.
        max_value=n_seq.max(axis=1).max()
        if (max_value)>7500:
            #Increase the axis to the next 500 increment above the max; increment by another 500 if the max will be close to the new interval.
            if max_value%500>250:
                ymax=((max_value//500)+2)*500
            else:
                ymax=((max_value//500)+1)*500
        else:
            ymax=7500
    elif log==False and y_max:
        ymax=y_max

    elif log==True:
        max_value=n_seq.max(axis=1).max()
        #Automatic Y-axis max calculation for log-10 scale
        #The y-axis max value will be 10^(order_of_magnitude_of_max_value+1)
        order_of_magnitude=(math.log10(max_value))//1
        #for a max value of 5,000, this will be 10,000, for a max value of 20,000, this will be 100,000, and so on.
        ymax=10**(order_of_magnitude+1)

    #If log==False, set y-axis limits with the maximum value defined above. 
    if log==False:
        ax.set_ylim(-20,ymax)
        #X limits: number of columns minus one (zero-index) plus 0.5
        ax.set_xlim(left=-0.5,right=n_seq.shape[1]-0.5)
    
    #If log==True, set y-axis limits using a log scale.
    elif log==True:
        #Create a logarithmic scale for the y-axis
        plt.yscale(value="symlog",linthresh=1,subs=[2, 3, 4, 5, 6, 7, 8, 9])
        #Set y axis labels to standard notation (default is scientific with log scales)
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0f'))
        #Y-axes limits: set based on max value
        ax.set_ylim(bottom=-0.02,top=ymax)
    
    #X-tick properties (ticks and label placement)
    ax.tick_params(which="major",axis='x',width=1.25,**x_tick_properties)
    #X-tick text properties
    plt.setp(ax.get_xticklabels(), ha="right",rotation_mode="anchor")

    #Y-axis tick/label parameters
    ax.tick_params(which="major",axis='y',width=1.25,**y_tick_properties)
    #plt.setp(ax.get_yticklabels(),**y_tick_properties)

    #Adjust y-axis ticks (only use for linear graphs)
    if log==False:
        ax.set_yticks(np.arange(0,ymax+1,y_tick_freq[0]))
        ax.set_yticks(np.arange(0,ymax+1,y_tick_freq[1]),minor=True);
    
    #Set axes labels
    ax.set_ylabel("Number of Sequences Analyzed",fontsize=12)
    ax.set_xlabel("Collection Date of Sample",fontsize=12)

    #Figure Title
    ax.set_title(f"Weekly Number of Sequences Analyzed by Continent: {protein}",**title_properties)
    
    #Create a grid
    ax.grid(which='major',axis='both',color="#DDDDDD",alpha=0.3)

    return lines