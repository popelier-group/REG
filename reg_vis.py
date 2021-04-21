import numpy as np
import pandas as pd
from scipy.stats import norm
import matplotlib.pyplot as plt
from adjustText import adjust_text
from matplotlib.ticker import MultipleLocator  

def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)

def plot_violin(file_list=[], save = False, file_name='reg_violin.png'):
    labels = [i+1 for i in range(len(file_list))]
    fig=plt.figure()
    fig.suptitle(r'Distribution of $R$ in each segment')
    ax = fig.add_subplot(111)
    ax.violinplot(file_list, showmeans='True')
    set_axis_style(ax, labels)
    ax.set_xlabel('Segments')
 
    if save == True:
        fig.savefig(file_name, dpi=300)

    return ax


def generate_data_vis(file, file_list, n_term, save = False, file_name='data_vis.png', title = "REG data visualization" ):

    (mu, sigma) = norm.fit(file['R'])
    fig = plt.figure(constrained_layout=False, figsize=[10,8])
    fig.suptitle(title)       
    gs1 = fig.add_gridspec(nrows=6, ncols=6, top=0.9, bottom=0.65, hspace=0, wspace = 0.40)
    ax1 = fig.add_subplot(gs1[0:6,0:2]) # violin
    ax2 = fig.add_subplot(gs1[0:3,2:6]) # hist
    ax3 = fig.add_subplot(gs1[3:6,2:6]) #steam
    gs2 = fig.add_gridspec(nrows=6, ncols=6, top=0.55, bottom=0.05, hspace=0.05)
    ax4 = fig.add_subplot(gs2[:,:]) # data

##Violin_plot
    labels = [i+1 for i in range(len(file_list))]
    parts = ax1.violinplot(file_list, showmeans='True')
    set_axis_style(ax1, labels)
    for i in range(len(parts['bodies'])):
        if file.R.equals(file_list[i]) == True:
            parts['bodies'][i].set_color('green')
            parts['bodies'][i].set_facecolor('green')
            parts['bodies'][i].set_edgecolor('green')
    ax1.set_xlabel('Segments')
    ax1.grid(axis='y', alpha=0.75)
    ax1.set_title(r'Distribution of $R$')



##histogram plot
    ax2.hist(file['R'], bins=[0.1*i for i in range(-11,11,1)],
             density=0 , color = 'green', edgecolor='black', linewidth=1.2, alpha=0.4)
    ax2.grid(axis='y', alpha=0.75)
    ax2.grid(axis='x', alpha=0.75)
    ax2.set_xlim(-1,1)
    ax2.xaxis.set_ticks(np.arange(-1.0, 1.1, 0.1))
    ax2.xaxis.set_ticks_position('top')
    ax2.xaxis.set_label_position('top')
    ax2.set_ylabel('Count Number')
    ax2.set_xlabel(r'$R$')
    ax2.yaxis.set_ticks_position('right')
    #ax2.axvline(mu, ls='--', color='r')

##stem plot1
    ax3.set_xlim(-1, 1)
    ax3.stem(file['R'], file['REG'],'k', markerfmt=' ', use_line_collection ='True')
    ax3.grid(axis='x', alpha=0.75)
    ax3.set_ylabel('REG value')
    ax3.xaxis.set_ticks(np.arange(-1.0, 1.1, 0.1))
    ax3.yaxis.set_ticks([])
    ax3.set_xlabel(r'$R$')
    #ax3.axvline(mu, ls='--', color='r')
##Stem plot2
    pos = pd.DataFrame(columns=file.columns)
    neg = pd.DataFrame(columns=file.columns)    
    cond = file.REG < 0
    rows = file.loc[cond,:]
    neg = neg.append(rows, ignore_index=True)    
    cond = file.REG > 0
    rows = file.loc[cond,:]
    pos = pos.append(rows, ignore_index=True)        
    markerline, stemlines, baseline = ax4.stem(pos['R'], pos['REG'], 'b', use_line_collection ='True')
    markerline.set_markerfacecolor('b')
    markerline.set_markeredgecolor('b')    
    markerline2, stemlines2, baseline2 = ax4.stem(neg['R'], neg['REG'], 'r', use_line_collection ='True')
    markerline2.set_markerfacecolor('r')
    markerline2.set_markeredgecolor('r')    
    ax4.set_xlim(-1, 1)
    bbox_props = dict(boxstyle="round,pad=0.1", fc="w", ec="k", lw=0, alpha=0.8)
    cond = pos.R > 0
    temp = pos.loc[cond,:]
    n = temp.nlargest(n_term, 'REG').reset_index(drop=True)  
    text_p = [ax4.text(n['R'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props) for i in range(len(n))]
    adjust_text(text_p, arrowprops=dict(arrowstyle='->', color='black', lw=1.5)) 

    cond = neg.R < 0
    temp = neg.loc[cond,:]
    n = temp.nsmallest(n_term, 'REG').reset_index(drop=True)  
    text_n = [ax4.text(n['R'][i], n['REG'][i], n['TERM'][i], bbox=bbox_props) for i in range(len(n))]
    adjust_text(text_n, arrowprops=dict(arrowstyle='->', color='black', lw=1.5))

    ax4.set_title('Relevant IQA contributions')
    ax4.set_xlabel(r'$R$')
    ax4.set_ylabel('REG value')
    
    if save == True:
        fig.savefig(file_name, dpi=300)
    
    return

def plot_segment(coordinate, wfn_energy, critical_points, label=False, color=True, annotate=True,
                 title='REG segments', y_label='Energy', x_label='Coordinate', save =False, file_name='segments.png'):
        
    color_list=['blue', 'darkgreen', 'darkred', 'yellow', 'cyan', 'magenta', 'grey', 'salmon',
                'seagreen', 'aquamarine', 'lightgreen', 'silver', 'lime', 'indigo', 'indianred']

    fig = plt.figure(figsize=(12,5))
    graph = plt.subplot(111)
    graph.plot(coordinate,wfn_energy,'o', color='black')
    if label == True:
        for i in range(len(wfn_energy)):
            graph.annotate(str(i+1), [coordinate[i], wfn_energy[i]])

    
    graph.xaxis.set_major_locator(MultipleLocator(2*abs(coordinate[0]-coordinate[1])))
    graph.xaxis.set_minor_locator(MultipleLocator(abs(coordinate[0]-coordinate[1])))    
    graph.set_title(title, fontsize =16)  
    graph.set_xlim([graph.get_xlim()[0], graph.get_xlim()[1]])
    graph.set_ylabel(y_label, fontsize =14)
    graph.set_xlabel(x_label, fontsize=14)
    graph.axvline(coordinate[0], linestyle = 'dotted', color = 'black') #the first segment starts at the global minimum.
    
    for i in critical_points: #split the segments
        graph.axvline(coordinate[i], linestyle = 'dotted', color = 'black')
          
    if color == True: #Color_segments
        graph.axvspan(graph.get_xlim()[0], graph.get_xlim()[1], facecolor = color_list[len(critical_points)], alpha = 0.2)
        start = 0
        for i in range(len(critical_points)):
            if coordinate[start] > coordinate[critical_points[i]]:
                graph.axvspan(coordinate[critical_points[i]], coordinate[start], facecolor='white', alpha =1)
                graph.axvspan(coordinate[critical_points[i]], coordinate[start], facecolor=color_list[i], alpha =0.2)
            start=critical_points[i]            
        start = 0   
        for i in range(len(critical_points)):
            if coordinate[start] < coordinate[critical_points[i]]:
                graph.axvspan(coordinate[start], coordinate[critical_points[i]], facecolor='white', alpha =1)
                graph.axvspan(coordinate[start], coordinate[critical_points[i]], facecolor=color_list[i], alpha =0.2)
            start=critical_points[i]
                
    if annotate== True: #Write segment number:
        text = []
        start=0
        y_pos = graph.get_ylim()[1] 
        x_pos = coordinate[start]
        text.append(graph.text(x_pos,  y_pos, '1',fontsize=14))
        
        for i in range(len(critical_points)):
            start=critical_points[i]
            x_pos = coordinate[start]
            text.append(graph.text(x_pos, y_pos,  str(i+2), fontsize=14))
            start=critical_points[i]
        adjust_text(text)
    fig.autofmt_xdate()
    if save == True:
        fig.savefig(file_name, dpi=300)

    return
    
def pandas_REG_dataframe_to_table(dataframe, table_name, SAVE_FIG=True):
    if SAVE_FIG==True:
        dataframe['R'] = np.round(dataframe['R'], decimals=3)
        dataframe['REG'] = np.round(dataframe['REG'], decimals=2)
        fig, ax = plt.subplots()
        ax.axis('off')
        ax.axis('tight')
        t= ax.table(cellText=dataframe.values, colWidths = [0.4]*len(dataframe.columns),  colLabels=dataframe.columns,  cellLoc='center',loc='center')
        t.auto_set_font_size(False)
        t.set_fontsize(12)
        fig.savefig(table_name, dpi=300, bbox_inches="tight")

def create_term_dataframe(reg_dataframe, headers, i):
    temp = [reg_dataframe[0][i], reg_dataframe[1][i]]
    df = pd.DataFrame(temp).transpose()
    df.columns = ["REG", "R"]
    df.index = headers
    df = df.rename_axis('TERM').reset_index()
    return df

def filter_term_dataframe(prop_dataframe, original_prop_name, new_prop_name):
    col = ['TERM','REG','R']
    new_df = pd.DataFrame()
    for j in range(0, len(prop_dataframe['TERM'])):
        if original_prop_name in prop_dataframe['TERM'][j]:
             new_df =  new_df.append(prop_dataframe.iloc[j])
    new_df = new_df[col]
    new_df = new_df.sort_values('REG').reset_index(drop=True)
    new_df['TERM'] = new_df['TERM'].apply(lambda row:row.replace(original_prop_name + '-', new_prop_name +'('))
    new_df['TERM'] = new_df['TERM'].str.replace("_", ',')
    new_df['TERM'] = new_df['TERM'] + ')'
    return new_df